
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


/* log likelihood of IV model with bivariate normal errors (without the constant parts) */ 
double logL_IV_normal(double *T, int *d, double *X, double *G, int n, double a0, double a1, double sig1, double b0, double b1, double sig2, double rho)
{
	double logL, mu1, mu2, bot=(1-rho*rho)*sig2;
	int i,  nd1=sumint(d, n);

	logL = (log(sig1)*n + (log(sig2)+log(1-rho*rho))*nd1) /(-2.0);
	for (i=0; i<n; i++)
	{
		mu1 = X[i]-a0-a1*G[i];	
		mu2 = b0+b1*X[i]-T[i]+sqrt(sig2/sig1)*rho*mu1;
		logL += mu1*mu1/sig1/(-2.0);
		if (d[i]==1) logL += mu2*mu2/bot/(-2.0);
		else if (d[i]==0) logL += log(pnorm(mu2/sqrt(bot))); 
	}
	return(logL);
}



/* log likelihood of IV model with bivariate normal errors and Dirichlet Process Mixtures (without the constant parts) */ 
double logL_IV_DPM_all(double *T, int *d, double *X, double *G, int n, double a1, double b1, 
					   double *mu1, double *sig1, double *mu2, double *sig2, double *rho)
{
	double logL=0, e1, e2, e3;
	int i;
	for (i=0; i<n; i++)
	{
		e1 = X[i]-mu1[i]-a1*G[i] ;
		e2 = mu2[i]+b1*X[i]-T[i]+sqrt(sig2[i]/sig1[i])*rho[i]*e1 ;
		e3 = 1-rho[i]*rho[i] ; 
		logL += (log(sig1[i]) + e1*e1/sig1[i]) / (-2.0) ;
		if (d[i]==1) logL += (log(e3) + log(sig2[i]) + e2*e2/e3/sig2[i]) / (-2.0) ;
		else if (d[i]==0) logL += log(pnorm(e2/sqrt(e3*sig2[i]))) ; 
	}
	return(logL);
}



/* likelihood of IV model for one subject i (without the constant part in stage 1)*/ 
double L_IV_ind(double Ti, int di, double Xi, double Gi, double a1, double b1, double mu1, double sig1, double mu2, double sig2, double rho)
{
	double logL, e1, e2, e3;
	e1 = Xi-mu1-a1*Gi ;
	e2 = mu2+b1*Xi-Ti+sqrt(sig2/sig1)*rho*e1 ;
	e3 = 1-rho*rho ; 
	logL = (log(sig1) + e1*e1/sig1) / (-2.0) ;
	if (di==1) logL += (log(2*PI) + log(e3) + log(sig2) + e2*e2/e3/sig2) / (-2.0) ;
	else if (di==0) logL += log(pnorm(e2/sqrt(e3*sig2))) ; 

	return(exp(logL));
}




/* log likelihood of IV model for arbitrary censoring with bivariate normal errors and Dirichlet Process Mixtures (without the constant parts) */ 
/* censoring indicator d = 1 if left-censored (Y<L); 2 if interval-censored (L<Y<R); 3 if right-censored (Y>R); 4 if event (L=Y=R) */
double logL_IV_AC_DPM_all(double *L, double *R, int *d, double *X, double *G, int n, double a1, double b1, 
					   double *mu1, double *sig1, double *mu2, double *sig2, double *rho)
{
	double logL=0, e1, e3, eL, eR, SL, SR;
	int i;
	for (i=0; i<n; i++)
	{
		e1 = X[i]-mu1[i]-a1*G[i] ;
		e3 = 1-rho[i]*rho[i] ;
		eL = mu2[i]+b1*X[i]-L[i]+sqrt(sig2[i]/sig1[i])*rho[i]*e1 ;
 		eR = mu2[i]+b1*X[i]-R[i]+sqrt(sig2[i]/sig1[i])*rho[i]*e1 ;
		SL = pnorm(eL/sqrt(e3*sig2[i]));
		SR = pnorm(eR/sqrt(e3*sig2[i]));

		logL += (log(sig1[i]) + e1*e1/sig1[i]) / (-2.0) ;
		if (d[i]==4) logL += (log(e3) + log(sig2[i]) + eL*eL/e3/sig2[i]) / (-2.0) ; // event 
		else if (d[i]==3) logL += log(SR) ;  // right-censored
		else if (d[i]==2) logL += log(SL - SR) ; // interval-censored
		else if (d[i]==1) logL += log(1 - SL);  // left-censored
		else printf("censoring indicator unknown\n");
	}
	return(logL);
}



/* likelihood of IV model for arbitrary censoring for one subject i (without the constant part in stage 1)*/ 
/* censoring indicator d = 1 if left-censored (Y<L); 2 if interval-censored (L<Y<R); 3 if right-censored (Y>R); 4 if event (L=Y=R) */ 
double L_IV_AC_ind(double L, double R, int d, double X, double G, 
				   double a1, double b1, double mu1, double sig1, double mu2, double sig2, double rho)
{
	double logL, e1, e3, eL, eR, SL, SR;
	e1 = X - mu1 - a1*G ;
	e3 = 1 - rho*rho ; 
	eL = mu2 + b1*X - L + sqrt(sig2/sig1)*rho*e1 ;
	eR = mu2 + b1*X - R + sqrt(sig2/sig1)*rho*e1 ;
	SL = pnorm(eL/sqrt(e3*sig2));
	SR = pnorm(eR/sqrt(e3*sig2));

	logL = (log(sig1) + e1*e1/sig1) / (-2.0) ;
	if (d==4) logL += (log(2*PI) + log(e3) + log(sig2) + eL*eL/e3/sig2) / (-2.0) ; // event 
	else if (d==3) logL += log(SR) ;  // right-censored
	else if (d==2) logL += log(SL - SR) ; // interval-censored
	else if (d==1) logL += log(1 - SL);  // left-censored
	else printf("censoring indicator unknown\n");
	return(exp(logL));
}


/* log likelihood of IV model for arbitrary censoring with bivariate normal errors (without the constant parts) */ 
/* censoring indicator d = 1 if left-censored (Y<L); 2 if interval-censored (L<Y<R); 3 if right-censored (Y>R); 4 if event (L=Y=R) */ 
double logL_IV_AC_normal(double *L, double *R, int *d, double *X, double *G, int n,
						 double a0, double a1, double sig1, double b0, double b1, double sig2, double rho)
{
	double logL, e1, e3, eL, eR, SL, SR;
	int i;
	logL = log(sig1) * n / (-2.0);
	for (i=0; i<n; i++)
	{
		e1 = X[i]-a0-a1*G[i];
		e3 = 1 - rho*rho;
		eL = b0 + b1*X[i] - L[i] + sqrt(sig2/sig1)*rho*e1 ;
		eR = b0 + b1*X[i] - R[i] + sqrt(sig2/sig1)*rho*e1 ;
		SL = pnorm(eL/sqrt(e3*sig2));
		SR = pnorm(eR/sqrt(e3*sig2));

		logL += e1*e1 / sig1 / (-2.0) ;
		if (d[i]==4) logL += (log(e3) + log(sig2) + eL*eL/e3/sig2) / (-2.0) ; // event 
		else if (d[i]==3) logL += log(SR) ;  // right-censored
		else if (d[i]==2) logL += log(SL - SR) ; // interval-censored
		else if (d[i]==1) logL += log(1 - SL);  // left-censored
		else printf("censoring indicator unknown\n");
	}
	return(logL);
}



