
/* quantiles of .025 and .975 from a sorted array */
void quantile95(double *X, int n, double q95[])
{
	int a1, a2;
	double p1=0.025, p2=0.975, h1, h2;

	h1 = (n-1)*p1+1;
	h2 = (n-1)*p2+1;
	a1 = (int)h1; 
	a2 = (int)h2;
	q95[0] = X[a1-1] + (X[a1]-X[a1-1])*(h1-a1);
	q95[1] = X[a2-1] + (X[a2]-X[a2-1])*(h2-a2);
}


/* quantiles of (.025, 0.50, .975) from a sorted array */
void quartiles(double *X, int n, double quar[])
{
	int a1, a2, a3;
	double p1=0.025, p2=0.5, p3=0.975, h1, h2, h3;

	h1 = (n-1)*p1+1;
	h2 = (n-1)*p2+1;
	h3 = (n-1)*p3+1;
	a1 = (int)h1; 
	a2 = (int)h2;
	a3 = (int)h3;
	quar[0] = X[a1-1] + (X[a1]-X[a1-1])*(h1-a1);
	quar[1] = X[a2-1] + (X[a2]-X[a2-1])*(h2-a2);
	quar[2] = X[a3-1] + (X[a3]-X[a3-1])*(h3-a3);
}



/* save individual values Xi by cluster indicator ci and cluster (distinct) values X_c */
void indvalues(double *Xi, int n, int *ci, double *X_c)
{
	int j;
	for (j=0; j<n; j++)
	{
		Xi[j]=X_c[ci[j]];
	}
}

/* save individual values Xi, Yi by cluster indicator ci and cluster (distinct) values X_c, Y_c, respectively */
void indvalues2(double *Xi, double *Yi, int n, int *ci, double *X_c, double *Y_c)
{
	int j, z;
	for (j=0; j<n; j++)
	{
		z=ci[j];
		Xi[j]=X_c[z];
		Yi[j]=Y_c[z];
	}
}


/* save individual values (mu1_i, sig1_i, mu2_i, sig2_i, rho_i) by cluster indicator ci and cluster (distinct) values (mu1_c, sig1_c, mu2_c, sig2_c, rho_c), respectively */
void indvalues5(double *mu1, double *sig1, double *mu2, double *sig2, double *rho, int n, int *ci, double *mu1_c, double *sig1_c, double *mu2_c, double *sig2_c, double *rho_c)
{
	int j, z;
	for (j=0; j<n; j++)
	{
		z=ci[j];
		mu1[j]=mu1_c[z];
		sig1[j]=sig1_c[z];
		mu2[j]=mu2_c[z];
		sig2[j]=sig2_c[z];
		rho[j]=rho_c[z];
	}
}

