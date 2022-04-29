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
//#include "IV_AC_DPM_MH.h"
#include "IV_AC_DPM_MH_temp.h"

int main(void)
{
	long idum=-2014021610;
	/***** sample size and width of random walk for (a1, b1, mu1, sig1, mu2, sig2, rho, v) *****/
	int sim=1;  // number of simulations
	int n=500;   //sample size (300, 500)
	//double wid[]={0.065, 0.31, 0.08, 0.035, 0.25, 0.30, 0.28, 0.4}; //width of random walk for n=300
	double wid[]={0.05, 0.22, 0.045, 0.023, 0.20, 0.25, 0.22, 0.3}; //width of random walk for n=500
	int i, j;
	double  *X, *G, *U, *W, *e1, *e2, *e3, *Y, *L, *R, *C1, *C2, *eC1, *eC2;
	int *d ;
	/***** simulation parameter settings *****/
	double a0t=0.5, a1t=sqrt(0.1*0.05), a2t=sqrt(0.1*0.4), e1sd=sqrt(0.1*0.4), e2sd=sqrt(0.1*0.15);
	double b1t=0; // b1t={0, -0.5, -1};
	double b0t=5-a0t*b1t, b2t=-sqrt(0.5*0.1), e3sd=sqrt(0.5*0.9);
	double sig1t=0.095, sig2t=0.5; //only used for vague informative priors;

	/***** priors *****/
	double prior_a1[]={0, 10000}, prior_b1[]={0, 10000};
	double prior_c1[]={0, 0.001, 0, 0.001}, prior_c2[]={10000, 0.001, 10000, 0.001};
	double infop_c1[]={a0t, sig1t, b0t, sig2t}, infop_c2[]={10, 5, 10, 5}; //vague informative priors for auxiliary parameters
	double prior_v[]={0.01, 3.3, 0.8}; // (v_low, v_up, v_w) in prior of v: proportional to ((v_up-v)/(v_up-v_low))^v_w  (Conley et al. 2008)
	//{0.01, 5, 0.8} for n=100; {0.01, 3.3, 0.8} for n=300; {0.01, 2.9, 0.8} for n=500; {0.01, 2.6, 0.8} for n=800; {0.01, 2.2, 0.8} for n=2000;

	FILE  *b1_result;
	double res_b1[6],  res_v[7], res_acp[9];

	X = (double *) malloc(n*sizeof(double));
	U = (double *) malloc(n*sizeof(double));
	G = (double *) malloc(n*sizeof(double));
	W = (double *) malloc(n*sizeof(double));
	e1 = (double *) malloc(n*sizeof(double));
	e2 = (double *) malloc(n*sizeof(double));
	e3 = (double *) malloc(n*sizeof(double));
	Y = (double *) malloc(n*sizeof(double));
	L = (double *) malloc(n*sizeof(double));
	R = (double *) malloc(n*sizeof(double));
	C1 = (double *) malloc(n*sizeof(double));
	C2 = (double *) malloc(n*sizeof(double));
	eC1 = (double *) malloc(n*sizeof(double));
	eC2 = (double *) malloc(n*sizeof(double));
	d = (int *) malloc(n*sizeof(int));
	
	b1_result=fopen("b1_results.txt", "a");
	fprintf(b1_result, "i	b1	b1_med	b1_se	b1_CIL	b1_CIR	b1_cover	v	v_med v_se	v_CIL	v_CIR	k	k_se\n");
	fclose(b1_result);

	for (j=0; j<sim; j++)
	{
		for(i=0; i<n; i++)
		{
			G[i] = gasdev(&idum);
			U[i] = gasdev(&idum);
			e1[i] = gasdev(&idum)*e1sd;
			e2[i] = gasdev(&idum)*e2sd;
			W[i] = a0t + a1t*G[i] + a2t*U[i] + e1[i]; 
			X[i] = W[i] + e2[i];

			e3[i] = gasdev(&idum)*e3sd;
			Y[i] = b0t + b1t*W[i] + b2t*U[i] + e3[i];
			eC1[i] = gasdev(&idum)*e3sd;
			C1[i] = b0t + b1t*W[i] + b2t*U[i] + eC1[i];
			eC2[i] = gasdev(&idum)*e3sd;
			C2[i] = b0t + b1t*W[i] + b2t*U[i] + eC2[i];
			if (eC1[i]<eC2[i]) {L[i]=C1[i]; R[i]=C2[i];}
			else {L[i]=C2[i]; R[i]=C1[i];}

			if (Y[i]<L[i]) d[i]=1;
			else if (Y[i]<R[i]) d[i]=2;
			else d[i]=3;
		}

		IV_AC_DPM_MH(X, L, R, d, G, n, 5500, 500, 5, prior_a1, prior_b1, prior_c1, prior_c2, infop_c1, infop_c2, prior_v, wid, 10, res_b1, res_v, res_acp, &idum); 
		   
		if(res_b1[3]>b1t || res_b1[4]<b1t) res_b1[5]=0; else res_b1[5]=1;
		
		printf("%d %f %f %f %f %f %f %f %f %f %f %f %f %f\n", j+1, res_b1[0], res_b1[1], res_b1[2], res_b1[3], res_b1[4], res_b1[5], res_v[0], res_v[1], res_v[2], res_v[3], res_v[4], res_v[5], res_v[6]);
		printf("acceptance rate %f %f %f %f %f %f %f %f %f \n", res_acp[0], res_acp[1], res_acp[2], res_acp[3], res_acp[4], res_acp[5], res_acp[6], res_acp[7], res_acp[8]);
		
		b1_result=fopen("b1_results.txt", "a");
		fprintf(b1_result, "%d %f %f %f %f %f %f %f %f %f %f %f %f %f\n", j+1, res_b1[0], res_b1[1], res_b1[2], res_b1[3], res_b1[4], res_b1[5], res_v[0], res_v[1], res_v[2], res_v[3], res_v[4], res_v[5], res_v[6]);
		fclose(b1_result);

		

	}
	
	free(X); free(G); free(W); free(U); free(e1); free(e2); free(e3); 
	free(Y); free(L); free(R); free(d); free(C1); free(C2); free(eC1); free(eC2);
	
	return 0;

}
