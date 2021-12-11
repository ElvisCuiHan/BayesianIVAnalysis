
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
double L_IV_ind(double Ti, double di, double Xi, double Gi, double a1, double b1, double mu1, double sig1, double mu2, double sig2, double rho)
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



/* log likelihood of IV-normal model with multiple instruments and multiple observed confounders (without the constant parts) */ 
double logL_IV_normal_cov(double *T, int *d, double *X, int n, double **G, int kG, double **U, int kU, 
						  double a0, double *a1, double *a2, double sig1, double b0, double b1, double *b2, double sig2, double rho)
{
	double logL, mu1, mu2, bot=(1-rho*rho)*sig2, a1G, a2U, b2U;
	int i, j,  nd1=sumint(d, n);

	logL = (log(sig1)*n + (log(sig2)+log(1-rho*rho))*nd1) /(-2.0);
	for (i=0; i<n; i++)
	{
		a1G=0; a2U=0; b2U=0;
		for (j=0; j<kG; j++)   a1G += a1[j]*G[i][j];
		for (j=0; j<kU; j++) { a2U += a2[j]*U[i][j];  b2U += b2[j]*U[i][j]; }
		mu1 = X[i] - a0 - a1G - a2U;	
		mu2 = b0 + b1*X[i] + b2U -T[i] + sqrt(sig2/sig1)*rho*mu1;
		logL += mu1*mu1/sig1/(-2.0);
		if (d[i]==1) logL += mu2*mu2/bot/(-2.0);
		else if (d[i]==0) logL += log(pnorm(mu2/sqrt(bot))); 
	}
	return(logL);
}



/* log likelihood of IV-DPM model with multiple instruments and multiple observed confounders (without the constant parts) */ 
double logL_IV_DPM_all_cov (double *T, int *d, double *X, int n, double **G, int kG, double **U, int kU,
							double b1, double *a1,  double *a2, double *b2,
							double *mu1, double *sig1, double *mu2, double *sig2, double *rho)
{
	double logL=0, e1, e2, e3, a1G, a2U, b2U;
	int i, j;
	for (i=0; i<n; i++)
	{
		a1G=0; a2U=0; b2U=0;
		for (j=0; j<kG; j++)   a1G += a1[j]*G[i][j];
		for (j=0; j<kU; j++) { a2U += a2[j]*U[i][j];  b2U += b2[j]*U[i][j]; }
		e1 = X[i] - mu1[i] - a1G - a2U ;
		e2 = mu2[i] + b1*X[i] + b2U -T[i] + sqrt(sig2[i]/sig1[i])*rho[i]*e1 ;
		e3 = 1-rho[i]*rho[i] ; 
		logL += (log(sig1[i]) + e1*e1/sig1[i]) / (-2.0) ;
		if (d[i]==1) logL += (log(e3) + log(sig2[i]) + e2*e2/e3/sig2[i]) / (-2.0) ;
		else if (d[i]==0) logL += log(pnorm(e2/sqrt(e3*sig2[i]))) ; 
	}
	return(logL);
}


/* likelihood of IV model for one subject i with multiple instruments and multiple observed confounders (without the constant part in stage 1)*/ 
double L_IV_ind_cov(double Ti, double di, double Xi, double *Gi, int kG, double *Ui, int kU, 
					double b1, double *a1,  double *a2, double *b2, double mu1, double sig1, double mu2, double sig2, double rho)
{
	double logL, e1, e2, e3, a1G, a2U, b2U;
	int j;
	a1G=0; a2U=0; b2U=0;
	for (j=0; j<kG; j++)   a1G += a1[j]*Gi[j];
	for (j=0; j<kU; j++) { a2U += a2[j]*Ui[j];  b2U += b2[j]*Ui[j]; }
	e1 = Xi - mu1 - a1G - a2U ;
	e2 = mu2 + b1*Xi + b2U -Ti + sqrt(sig2/sig1)*rho*e1 ;
	e3 = 1-rho*rho ; 
	logL = (log(sig1) + e1*e1/sig1) / (-2.0) ;
	if (di==1) logL += (log(2*PI) + log(e3) + log(sig2) + e2*e2/e3/sig2) / (-2.0) ;
	else if (di==0) logL += log(pnorm(e2/sqrt(e3*sig2))) ; 
	return(exp(logL));
}
