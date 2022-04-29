// static long double const PI = 3.14159265358979323846264338L ;
// #define PI 3.14159265358979323846

#define LDMAX(a, b) (((a) > (b)) ? (a) : (b))

/* log likelihood of log-normal AFT model */ 
double logL_AFT_logN(double *X, double *T, int *d, int n, double a0, double a1, double sig2)
{
	double logL=0, mu;
	int i;
	for (i=0; i<n; i++)
	{
		mu = T[i]-a0-a1*X[i];
		if (d[i]==1) {logL += (log(2*PI*sig2) + mu*mu/sig2)/(-2.0) ;}
		else if (d[i]==0) {logL += log(1-pnorm(mu/sqrt(sig2))); }
	}
	return(logL);
}

/* log likelihood of IV model for arbitrary censoring with bivariate normal errors and */
/* Dirichlet Process Mixtures (without the constant parts) */
/* censoring indicator d = 1 if left-censored (Y<L); 2 if interval-censored (L<Y<R); 
						   3 if right-censored (Y>R); 4 if event (L=Y=R) */

long double logL_IV_DPM_UKB(double *L, double *R, int *d, double *X, int n, double **G, int kG, double **U, int kU,
							double b1, double *a1,  double *a2, double *b2,
							double *mu1, double *sig1, double *mu2, double *sig2, double *rho)
{
	long double logL=0, e1, eL, eR, e3, log_fx, log_fy, SL, SR, a1G, a2U, b2U;
	int i, j;
	for (i=0; i<n; i++)
	{
		a1G=0; a2U=0; b2U=0;
		for (j=0; j<kG; j++)   a1G += a1[j]*G[i][j];
		for (j=0; j<kU; j++) { a2U += a2[j]*U[i][j];  b2U += b2[j]*U[i][j]; }
		e1 = X[i] - mu1[i] - a1G - a2U ;
		//e1 = X[i] - mu1[i] - a1G;
		//e2 = mu2[i] + b1*X[i] + b2U -T[i] + sqrt(sig2[i]/sig1[i])*rho[i]*e1 ;
		eL = mu2[i] + b1 * X[i] + b2U - L[i] + sqrt(sig2[i] / sig1[i]) * rho[i] * e1 ;
 		eR = mu2[i] + b1 * X[i] + b2U - R[i] + sqrt(sig2[i] / sig1[i]) * rho[i] * e1 ;
		//eL = mu2[i] + b1 * X[i] - L[i] + sqrt(sig2[i] / sig1[i]) * rho[i] * e1 ;
 		//eR = mu2[i] + b1 * X[i] - R[i] + sqrt(sig2[i] / sig1[i]) * rho[i] * e1 ;
		e3 = 1 - rho[i] * rho[i] ; 
		log_fx = (log(sig1[i]) + e1 * e1 / sig1[i]) / (-2.0) ;
		//log_fy = (log(e3) + log(sig2[i]) + eL * eL / e3 / sig2[i]) / (-2.0) ;
		log_fy = (log(e3) + log(sig2[i]) + eR * eR / e3 / sig2[i]) / (-2.0) ;
		SL = (pnorm(eL / sqrt(e3 * sig2[i]))) ;
		SR = (pnorm(eR / sqrt(e3 * sig2[i]))) ;

		logL += log_fx ;
		if (d[i]==4) logL += log_fy ; // Event
		//else if (d[i]==3) logL += log(SR) ; // Right-censored
		//else if (d[i]==2) logL += log(SL - SR) ; // Interval-censored
		//else if (d[i]==1) logL += log(1 - SL) ; // Left-censored
		else if (d[i]==3) logL += log(LDMAX(SR, 1e-320)) ; // Right-censored
		else if (d[i]==2) logL += log(LDMAX(SL - SR, 1e-320)) ; // Interval-censored
		else if (d[i]==1) logL += log(LDMAX(1.0 - SL,1e-320)) ; // Left-censored
		else printf("The censoring indicator is unknown!\n") ;
	}
	return(logL);
}

/* likelihood of IV model for one subject i with multiple instruments and multiple observed confounders (without the constant part in stage 1)*/ 
/* censoring indicator d = 1 if left-censored (Y<L); 2 if interval-censored (L<Y<R); 3 if right-censored (Y>R); 4 if event (L=Y=R) */
double L_IV_ind_UKB(double Li, double Ri, double di, double Xi, double *Gi, int kG, double *Ui, int kU, 
					double b1, double *a1,  double *a2, double *b2, double mu1, double sig1, double mu2, double sig2, double rho)
{
	double logL, e1, eL, eR, e3, log_fx, log_fy, SL, SR, a1G, a2U, b2U;
	int j;
	a1G=0; a2U=0; b2U=0;
	for (j=0; j<kG; j++)   a1G += a1[j]*Gi[j];
	for (j=0; j<kU; j++) { a2U += a2[j]*Ui[j];  b2U += b2[j]*Ui[j]; }
	e1 = Xi - mu1 - a1G - a2U ;
	//e1 = Xi - mu1 - a1G;
	eL = mu2 + b1 * Xi + b2U - Li + sqrt(sig2 / sig1) * rho * e1 ;
	eR = mu2 + b1 * Xi + b2U - Ri + sqrt(sig2 / sig1) * rho * e1 ;
	//eL = mu2 + b1 * Xi - Li + sqrt(sig2 / sig1) * rho * e1 ;
	//eR = mu2 + b1 * Xi - Ri + sqrt(sig2 / sig1) * rho * e1 ;
	e3 = 1 - rho * rho ; 
	log_fx = (log(2 * PI) + log(sig1) + e1 * e1 / sig1) / (-2.0) ;
	//log_fy = (log(2 * PI) + log(e3) + log(sig2) + eL * eL / e3 / sig2) / (-2.0) ;
	log_fy = (log(2 * PI) + log(e3) + log(sig2) + eR * eR / e3 / sig2) / (-2.0) ;
	SL = (pnorm(eL / sqrt(e3 * sig2))) ;
	SR = (pnorm(eR / sqrt(e3 * sig2))) ;

	logL = log_fx ; // (log(sig1) + e1*e1/sig1) / (-2.0) ;
	if (di==4) logL += log_fy ; // Event
	//else if (di==3) logL += log(SR) ; // Right-censored
	//else if (di==2) logL += log(SL - SR) ; // Interval-censored
	//else if (di==1) logL += log(1 - SL) ; // Left-censored
	else if (di==3) logL += log(LDMAX(SR, 1e-320)) ; // Right-censored
	else if (di==2) logL += log(LDMAX(SL - SR, 1e-320)) ; // Interval-censored
	else if (di==1) logL += log(LDMAX(1.0 - SL,1e-320)) ; // Left-censored
	else printf("The censoring indicator is unknown!\n") ;
	return(exp(logL));
}

/* log likelihood of IV-normal model with multiple instruments and multiple observed confounders (without the constant parts) */ 
double logL_IV_normal_UKB(double *L, double *R, int *d, double *X, int n, double **G, int kG, double **U, int kU, 
						  double a0, double *a1, double *a2, double sig1, double b0, double b1, double *b2, double sig2, double rho)
{
	// double logL, mu1, mu2, bot=(1-rho * rho) * sig2, a1G, a2U, b2U;
	double logL=0, e1, eL, eR, e3, log_fx, log_fy, SL, SR, a1G, a2U, b2U;
	int i, j,  nd1=sumint(d, n);

	// logL = (log(sig1) * n + (log(sig2) + log(1 - rho * rho)) * nd1) /(-2.0);
	for (i=0; i<n; i++)
	{
		a1G=0; a2U=0; b2U=0;
		for (j=0; j<kG; j++)   a1G += a1[j] * G[i][j];
		for (j=0; j<kU; j++) { a2U += a2[j] * U[i][j];  b2U += b2[j] * U[i][j]; }
		// mu1 = X[i] - a0 - a1G - a2U ;	
		// mu2 = b0 + b1 * X[i] + b2U - L[i] + sqrt(sig2 / sig1) * rho * mu1 ;
		// logL += mu1 * mu1 / sig1 / (-2.0) ;

		e1 = X[i] - a0 - a1G - a2U ;
		eL = b0 + b1 * X[i] + b2U - L[i] + sqrt(sig2 / sig1) * rho * e1 ;
 		eR = b0 + b1 * X[i] + b2U - R[i] + sqrt(sig2 / sig1) * rho * e1 ;
		//eL = b0 + b1 * X[i] - L[i] + sqrt(sig2 / sig1) * rho * e1 ;
 		//eR = b0 + b1 * X[i] - R[i] + sqrt(sig2 / sig1) * rho * e1 ;
		
		e3 = 1 - rho * rho ; 
		log_fx = (log(sig1) + e1 * e1 / sig1) / (-2.0) ;
		//log_fy = (log(e3) + log(sig2) + eL * eL / e3 / sig2) / (-2.0) ;
		log_fy = (log(e3) + log(sig2) + eR * eR / e3 / sig2) / (-2.0) ;
		SL = (pnorm(eL / sqrt(e3 * sig2))) ;
		SR = (pnorm(eR / sqrt(e3 * sig2))) ;

		logL += log_fx;
		if (d[i]==4) logL += log_fy ; // Event
		//else if (d[i]==3) logL += log(SR) ; // Right-censored
		//else if (d[i]==2) logL += log(SL - SR) ; // Interval-censored
		//else if (d[i]==1) logL += log(1 - SL) ; // Left-censored
		else if (d[i]==3) logL += log(LDMAX(SR, 1e-320)) ; // Right-censored
		else if (d[i]==2) logL += log(LDMAX(SL - SR, 1e-320)) ; // Interval-censored
		else if (d[i]==1) logL += log(LDMAX(1.0 - SL,1e-320)) ; // Left-censored		
		else printf("The censoring indicator is unknown!\n") ;
		// if (d[i]==1) logL += mu2*mu2/bot/(-2.0);
		// else if (d[i]==0) logL += log(pnorm(mu2/sqrt(bot))); 
	}
	return(logL);
}