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
	long idum=-2014012901;   // ???????????????

	srand((unsigned)time(NULL));

	/***** sample size *****/
	int n=13936, kG=15, kU=5;   //data dimensions: number of subjects, number of instruments, number of observed confounders
	int ch=1100000; //length of chain
	int nsave=50; // number of samples in each result saving
	int m=10;      //number of auxiliary parameters;
	int i,j;
	FILE  *data;
	double  *X, *L, *R, *tmp, **G, **U;
	int *d ;
	double prior_b1[2], *prior_a1_1, *prior_a1_2, *prior_a2_1, *prior_a2_2, *prior_b2_1, *prior_b2_2;
	double prior_c1[4], prior_c2[4], infop_c1[4], infop_c2[4], prior_v[3];

	/***** width of random walk (0.0~0.1) *****/
	//double width_b1, *width_a1, *width_a2, *width_b2, *width_c, width_v;
	double width_b1, *width_a1, *width_a2, *width_b2,  width_v;
	width_a1 = (double *) malloc(kG*sizeof(double));
	width_a2 = (double *) malloc(kU*sizeof(double));
	width_b2 = (double *) malloc(kU*sizeof(double));
	//width_c  = (double *) malloc(kU*sizeof(double));  // ????????????????????????
	width_b1 = 0.0384;
	//for (i=0; i<kG; i++) width_a1[i] = 0.0064;
	//for (i=0; i<kU; i++) width_a2[i] = 0.0108;
	for (i=0; i<kG; i++) width_a1[i] = 0.0128 * 0.5;
	for (i=0; i<kU; i++) width_a2[i] = 0.0128 * 0.75;
	for (i=0; i<kU; i++) width_b2[i] = 0.140;
	//for (i=0; i<kU; i++) width_c[i]  = 0.1;          // ????????????????????????
	double width_c[]={0.024, 0.0072, 0.480, 0.88, 0.075};
	width_v=0.4;

	/***** initial values *****/
	double init_b1, init_v;
	double *init_a1, *init_a2, *init_b2, **init_thetac;
	int init_ncl, *init_ci, *init_freqcl;

	X = (double *) malloc(n*sizeof(double));
	L = (double *) malloc(n*sizeof(double));
	R = (double *) malloc(n*sizeof(double));
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
	char fname[256] = "temp_init_values_male.txt";
	//char fname[256] = "temp_init_values_snp.txt";
	data=fopen(fname, "r");
	fscanf(data, "%lf\n", &init_b1);
	for (i=0; i<kG; i++) fscanf(data, "%lf", &init_a1[i]);
	for (i=0; i<kU; i++) fscanf(data, "%lf", &init_a2[i]);
	for (i=0; i<kU; i++) fscanf(data, "%lf", &init_b2[i]);
	fscanf(data, "%lf\n", &init_v);
	fscanf(data, "%d\n", &init_ncl);
	init_thetac = (double **) malloc(init_ncl*sizeof(double *));
	for (i=0; i<init_ncl; i++) init_thetac[i]=(double *) malloc(5*sizeof(double));  // 5 --- ????????
	init_freqcl = (int *) malloc(init_ncl*sizeof(int *));
	for (i=0; i<init_ncl; i++) fscanf(data, "%lf%lf%lf%lf%lf%d\n", &init_thetac[i][0], &init_thetac[i][1], &init_thetac[i][2], &init_thetac[i][3], &init_thetac[i][4], &init_freqcl[i]);
	//for (i=0; i<n; i++) fscanf(data, "%d", &init_ci[i]);
	for (i=0; i<n; i++) init_ci[i] = 0;  // set all 9865 init_ci to 0.
	fclose(data);

	// incremental b1 with stepsize = 0.025 from -0.5 to 0.5
 	//if (init_b1 < -5.0 || init_b1 > 4.0) init_b1 = -5.0;
	//else init_b1 += 1.0;
/* 	data=fopen(fname, "w");
	fprintf(data, "%lf\n", init_b1);
 	for(i=0; i<kG; i++) fprintf(data, "%lf ", init_a1[i]); fprintf(data,"\n"); 
	for(i=0; i<kU; i++) fprintf(data, "%lf ", init_a2[i]); fprintf(data,"\n");
	for(i=0; i<kU; i++) fprintf(data, "%lf ", init_b2[i]); fprintf(data,"\n");
	fprintf(data, "%lf\n", init_v);
	fprintf(data, "%d\n", init_ncl);
	for (i=0; i<init_ncl; i++) fprintf(data, "%lf %lf %lf %lf %lf %d\n", init_thetac[i][0], init_thetac[i][1], init_thetac[i][2], init_thetac[i][3], init_thetac[i][4], init_freqcl[i]);
	for (i=0; i<n; i++) fprintf(data, "%d ", init_ci[i]);
	fclose(data); */

	/***** priors *****/
	const double prior_mean = 0.0;         // mean of normal prior
	const double prior_var = 1;    // variance of normal prior
	prior_b1[0] = prior_mean - 1; prior_b1[1] = prior_var;                         // mean and variance of normal prior for b1;
	for (i=0; i<kG; i++) { prior_a1_1[i] = prior_mean;  prior_a1_2[i] = prior_var; } // mean and variance of normal prior for a1;
	for (i=0; i<kU; i++) { prior_a2_1[i] = prior_mean;  prior_a2_2[i] = prior_var; } // mean and variance of normal prior for a2;
	for (i=0; i<kU; i++) { prior_b2_1[i] = prior_mean;  prior_b2_2[i] = prior_var; } // mean and variance of normal prior for b2;
	//prior_c1[0]=0.0;     prior_c1[1]=0.0001; prior_c1[2]=0.0;       prior_c1[3]=0.0001;
	//prior_c2[0]=10000.0; prior_c2[1]=0.0001; prior_c2[2]=30000.0; prior_c2[3]=0.0001; 

	/*// Female
	prior_c1[0]=4.9701;  prior_c1[1]=0.1249; prior_c1[2]=4.875;       prior_c1[3]=.143;
	prior_c2[0]=0.0149; prior_c2[1]=1.0; prior_c2[2]=1; prior_c2[3]=1; 
	// slightly informative priors to generate a0, sig1, b0 and sig2 for the new clusters
	infop_c1[0]=4.9701;   infop_c1[1]=0.1249;  infop_c1[2]=4.875;    infop_c1[3]=.143;  // means (from the normal IV method)
	infop_c2[0]=0.0149;    infop_c2[1]=1.0;    infop_c2[2]=1;     infop_c2[3]=1;    // SD
	// prior for v
	prior_v[0]=0.01;     prior_v[1]=4.8;     prior_v[2]=2.4;
	*/

	// Male
	prior_c1[0]=4.9701;  prior_c1[1]=0.1249; prior_c1[2]=4.296;       prior_c1[3]=.0643;
	prior_c2[0]=0.0149;  prior_c2[1]=1.0; 	 prior_c2[2]=1; 		  prior_c2[3]=1; 
	// slightly informative priors to generate a0, sig1, b0 and sig2 for the new clusters
	infop_c1[0]=4.9701;   infop_c1[1]=0.1249;  infop_c1[2]=4.296;    infop_c1[3]=.0643;  // means (from the normal IV method)
	infop_c2[0]=0.0149;    infop_c2[1]=1.0;    infop_c2[2]=1;     infop_c2[3]=1;    // SD
	// prior for v
	prior_v[0]=0.01;     prior_v[1]=4.8;     prior_v[2]=2.4;

	data=fopen("0330_DPMIV_test_data_eth1_male_subsetsnp_2clust.dat", "r"); // ????????????????
	char idxnum[256];
	for (i=0; i<n; i++)
    {
		//fscanf(data, "%s", idxnum);
		//printf("%s\n", idxnum);
		for (j=0; j<4+kG+kU; j++) 
		{
			fscanf(data, "%lf", tmp+j);  // ?????????????????
		}
		//fscanf(data, "%lf%lf%d%lf%lf%lf%lf%lf%lf%lf\n", \
		//    &L[i], &R[i], &d[i], &X[i], \
		//	&tmp[0], &tmp[1], &tmp[2], &tmp[3], &tmp[4], &tmp[5]);
		L[i] = tmp[0];  R[i] = tmp[1];  d[i] = tmp[2];  X[i] = tmp[3];
		for (j=0;j<kG;j++) G[i][j]=tmp[4+j];      // tmp[j];     ???????????
		for (j=0;j<kU;j++) U[i][j]=tmp[4+kG+j];   // tmp[kG+j];
		//for (j=0;j<kU;j++) U[i][j]=tmp[j];
		//for (j=0;j<kG;j++) G[i][j]=tmp[kU+j];
    } fclose(data);
	init_a1[0] = -0.00120427914499857;

/* print out data to make sure the data is read in
for (i=0; i<n; i++)
	{
	printf("%f	%d	%f%	", T[i],d[i],X[i]);
		for (j=0;j<kG;j++)	printf("%f	", G[i][j]);
		for (j=0;j<kU;j++)	printf("%f	", U[i][j]);
	printf("\n");
	}
*/

    printf("starting.....\n");
	time_t timep;
	struct tm *p;
	char tstamp[256] = {0};
	time(&timep);
	p = localtime(&timep);
	sprintf(tstamp, "%02d%02d-%02d%02d%02d_ch=%dK_b1=%.3f", 1+p->tm_mon, p->tm_mday, p->tm_hour, p->tm_min, p->tm_sec, ch/1000, init_b1);
	// MCMC sampling
	IV_DPM_MH_covar_tmp(X, L, R, d, n, kG, kU, G, U, ch, prior_b1, prior_a1_1, prior_a1_2, prior_a2_1, prior_a2_2, prior_b2_1, prior_b2_2,
					prior_c1, prior_c2, infop_c1, infop_c2, prior_v, width_b1, width_a1, width_a2, width_b2, width_c, width_v,
					m, init_b1, init_a1, init_a2, init_b2, init_ncl, init_thetac, init_v, init_freqcl, init_ci, &idum, nsave, tstamp);
					


  return 0;

}
