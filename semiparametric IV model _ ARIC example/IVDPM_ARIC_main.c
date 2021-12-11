#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>

#include "misc4.c"
#include "misc4.h"
#include "RAN1.c"
#include "GASDEV.c"
#include "GetRandom.c"
#include "log_likelihood.c"
#include "IV_DPM_covar_multidays_tmp.h"

int main(void)
{
	long idum=-2014012901;

	/***** sample size *****/
	int n=12625, kG=1, kU=7;   //data dimensions: number of subjects, number of instruments, number of observed confounders
	int ch=1100; //length of chain
	int nsave=100; // number of samples in each result saving
	int m=10; //number of auxiliary parameters;
	int	i,j;
	FILE  *data;
	double  *X, *T, *tmp, **G, **U;
	int *d ;
	double prior_b1[2], *prior_a1_1, *prior_a1_2, *prior_a2_1, *prior_a2_2, *prior_b2_1, *prior_b2_2;
	double prior_c1[4], prior_c2[4], infop_c1[4], infop_c2[4], prior_v[3];

	/***** width of random walk *****/
	double width_b1=0.12;
	double width_a1[]={0.001};
	double width_a2[]={0.0064,0.0094,0.000084,0.00084,0.0036,0.010,0.0145};
	double width_b2[]={.81, 1.00, .0091, .093, .44, 1.00, 1.20};
	double width_c[]={0.005, 0.00062, 0.60, 2.2, 0.05};
	double width_v=0.3;

	/***** initial values *****/
	double init_b1, init_v;
	double *init_a1, *init_a2, *init_b2, **init_thetac;
	int init_ncl, *init_ci, *init_freqcl;

	X = (double *) malloc(n*sizeof(double));
	T = (double *) malloc(n*sizeof(double));
	d = (int *) malloc(n*sizeof(int));
	G = (double **) malloc(n*sizeof(double *));
	for (i=0; i<n; i++) G[i]=(double *) malloc(kG*sizeof(double));
	U = (double **) malloc(n*sizeof(double *));
	for (i=0; i<n; i++) U[i]=(double *) malloc(kU*sizeof(double));
	prior_a1_1 = (double *) malloc(kG*sizeof(double));
	prior_a1_2 = (double *) malloc(kG*sizeof(double));
	prior_a2_1 = (double *) malloc(kU*sizeof(double));
	prior_a2_2 = (double *) malloc(kU*sizeof(double));
	prior_b2_1 = (double *) malloc(kU*sizeof(double));
	prior_b2_2 = (double *) malloc(kU*sizeof(double));
	tmp = (double *) malloc((kG+kU)*sizeof(double)); // temporary space to scan data
	init_a1 = (double *) malloc(kG*sizeof(double));
	init_a2 = (double *) malloc(kU*sizeof(double));
	init_b2 = (double *) malloc(kU*sizeof(double));
	init_ci = (int *) malloc(n*sizeof(int));

	// read in initial values (including day 1)
data=fopen("temp_init_values.txt", "r");
	fscanf(data, "%lf\n", &init_b1);
	fscanf(data, "%lf\n", &init_a1[0]);
	fscanf(data, "%lf%lf%lf%lf%lf%lf%lf\n", &init_a2[0], &init_a2[1], &init_a2[2], &init_a2[3], &init_a2[4], &init_a2[5], &init_a2[6]);
	fscanf(data, "%lf%lf%lf%lf%lf%lf%lf\n", &init_b2[0], &init_b2[1], &init_b2[2], &init_b2[3], &init_b2[4], &init_b2[5], &init_b2[6]);
	fscanf(data, "%lf\n", &init_v);
	fscanf(data, "%d\n", &init_ncl);
	init_thetac = (double **) malloc(init_ncl*sizeof(double *));
	for (i=0; i<init_ncl; i++) init_thetac[i]=(double *) malloc(5*sizeof(double));
	init_freqcl = (int *) malloc(init_ncl*sizeof(int *));
	for (i=0; i<init_ncl; i++) 
	fscanf(data, "%lf%lf%lf%lf%lf%d\n", &init_thetac[i][0], &init_thetac[i][1], &init_thetac[i][2], &init_thetac[i][3], &init_thetac[i][4], &init_freqcl[i]);
	for (i=0; i<n; i++) fscanf(data, "%d", &init_ci[i]);
	fclose(data);


	/***** priors *****/
	prior_b1[0]=0; prior_b1[1]=1000000; // mean and variance of normal prior for b1;
	for (i=0; i<kG; i++) { prior_a1_1[i]=0;  prior_a1_2[i]=100000; } // mean and variance of normal prior for a1;
	for (i=0; i<kU; i++) { prior_a2_1[i]=0;  prior_a2_2[i]=100000; } // mean and variance of normal prior for a2;
	for (i=0; i<kU; i++) { prior_b2_1[i]=0;  prior_b2_2[i]=100000; } // mean and variance of normal prior for b2;
	prior_c1[0]=0; prior_c1[1]=0.0001; prior_c1[2]=0; prior_c1[3]=0.0001;
	prior_c2[0]=10000; prior_c2[1]=0.0001; prior_c2[2]=3000000; prior_c2[3]=0.0001;

	// slightly informative priors to generate a0, sig1, b0 and sig2 for the new clusters
	infop_c1[0]=1.667; infop_c1[1]=0.008; infop_c1[2]=67.391;  infop_c1[3]=19.543;  // means (from the normal IV method)
	infop_c2[0]=10; infop_c2[1]=5; infop_c2[2]=100; infop_c2[3]=20;  // SD
	// prior for v
	prior_v[0]=0.01; prior_v[1]=1.7; prior_v[2]=0.8;

	data=fopen("ARIC_surv_V2.dat", "r");
	for (i=0; i<n; i++)
    {
		fscanf(data, "%lf%d%lf%lf%lf%lf%lf%lf%lf%lf%lf\n", &T[i], &d[i], &X[i], &tmp[0], &tmp[1], &tmp[2], &tmp[3], &tmp[4], &tmp[5], &tmp[6], &tmp[7]);
		for (j=0;j<kG;j++) G[i][j]=tmp[j];
		for (j=0;j<kU;j++) U[i][j]=tmp[kG+j];
    } fclose(data);

/* print out data to make sure the data is read in
for (i=0; i<n; i++)
	{
	printf("%f	%d	%f%	", T[i],d[i],X[i]);
		for (j=0;j<kG;j++)	printf("%f	", G[i][j]);
		for (j=0;j<kU;j++)	printf("%f	", U[i][j]);
	printf("\n");
	}
*/

	
	// MCMC sampling
	IV_DPM_MH_covar(X, T, d, n, kG, kU, G, U, ch, prior_b1, prior_a1_1, prior_a1_2, prior_a2_1, prior_a2_2, prior_b2_1, prior_b2_2,
					prior_c1, prior_c2, infop_c1, infop_c2, prior_v, width_b1, width_a1, width_a2, width_b2, width_c, width_v,
					m, init_b1, init_a1, init_a2, init_b2, init_ncl, init_thetac, init_v, init_freqcl, init_ci, &idum, nsave);
					


  return 0;

}
