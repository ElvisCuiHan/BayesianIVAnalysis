#include <math.h>
#include "gamma_function.c"
#include "ln_gamma_function.c"

/********* Function for IV linear regression for right-censored survival outcome with Dirichlet Process Mixtures *********/
/********* with multiple instruments and multiple observed confounders *********/
/*  Revision date: 2014-01-28
	Use Algorithm 8 in Neal (2000)   
	Use prior for strength parameter in DPM proposed by Conley et al. (2008)
	Note: This version saves the posterior samples of b1 (causal effect of true covariate W on outcome Y), 
	v (strength parameter of Dirichlet prior) and ncl (numebr of clusters) into a text file
	X = an array of univariate covariate (that might be subject to measurement error and/or unobserved confounding)
	T = an array of observed survival time
	d = an array of censoring indicator (1=event, 0=right-censored)
	n = sample size
	kG = number of instrumental variables
	kU = number of observed confounders
	G = a matrix (double array) of instrumental variables (dimension = n*kG)
	U = a matrix (double array) of observed confounders (dimension = n*kU)
	ch = number of iterations (length of the chain)
	prior_b1 = an array of the mean and variance for the normal prior for b1, e.g. {0, 10000}
	prior_a1_1 = an array of the means for the normal prior for a1, e.g. {0, 0, ..., 0} (length = kG)
	prior_a1_2 = an array of the variance for the normal prior for a1, e.g. {10000, 10000, ..., 10000} (length = kG)
	prior_a2_1 = an array of the means for the normal prior for a2, e.g. {0, 0, ..., 0} (length = kU)
	prior_a2_2 = an array of the variance for the normal prior for a2, e.g. {10000, 10000, ..., 10000} (length = kU)
	prior_b2_1 = an array of the means for the normal prior for b2, e.g. {0, 0, ..., 0} (length = kU)
	prior_b2_2 = an array of the variance for the normal prior for b2, e.g. {10000, 10000, ..., 10000} (length = kU)
    note: the sig1 and sig2 in this program refer to sigma_1-squared and sigma_2-squared in the paper
	prior_c1 = an array of the first parameter of the priors for (mu1c, sig1c, mu2c, sig2c), e.g. {0, 0.001, 0, 0.001}
			mean of the normal prior for mu1c and mu2c; shape parameter of the inverse-gamma prior for sig1c and sig2c;
	prior_c2 = an array of the second parameter of the priors for (mu1c, sig1c, mu2c, sig2c), e.g. {10000, 0.001, 10000, 0.001}
			variance of the normal prior for mu1c and mu2c; scale parameter of the inverse-gamma prior for sig1c and sig2c; 
	infop_c1 = an array of the first parameter of the (slightly) informative priors for (mu1c, sig1c, mu2c, sig2c), e.g. {a0t, sig1t, b0t, sig2t}
			mean of the normal prior for mu1c and mu2c; mean of the inverse-gamma prior for sig1c and sig2c
	infop_c2 = an array of the second parameter of the (slightly) informative priors for (mu1c, sig1c, mu2c, sig2c), e.g. {10, 5, 10, 5}
			SD of the normal prior for mu1c and mu2c; SD of the inverse-gamma prior for sig1c and sig2c
	prior_v = a vector of the parameters of the priors for v (the base/strength parameter of Dirichlet prior)
			(v_low, v_up, v_w):  prior of v is proportional to ((v_up-v)/(v_up-v_low))^v_w  (Conley et al. 2008)
    width_b1 = width of random walk for b1
	width_a1 = an array of the width of random walk for a1
	width_a2 = an array of the width of random walk for a2
	width_b2 = an array of the width of random walk for b2
	width_c = an array of the width of random walk for (mu1c,sig1c,mu2c,sig2c,rhoc)
	width_v = width of random walk for v
	m = number of auxiliary parameters in Dirichlet process sampling, e.g. 10
	init_b1 = initial value for b1
	init_a1 = an array of initial values for a1 (length=kG)
	init_a2 = an array of initial values for a2 (length=kU)
	init_b2 = an array of initial values for b2 (length=kU)
	init_ncl = initial value of ncl (number of clusters)
	init_thetac = a double array of initial values for (mu1c,sig1c,mu2c,sig2c,rhoc) (dim: init_ncl*5)
	init_v = initial value for v
	init_freqcl = an array of initial frequencies of clusters
	init_ci = an array of initial cluster indicators (length=n)
	idum = seed number for random normal sample generator
	nsave = number of samples to save each time (for running on Hoffman2 cluster)
*/

void IV_DPM_MH_covar(double *X, double *T, int *d, int n, int kG, int kU, double **G, double **U, int ch, 
					 double prior_b1[], double *prior_a1_1,  double *prior_a1_2, double *prior_a2_1,  double *prior_a2_2,
					 double *prior_b2_1,  double *prior_b2_2, double prior_c1[], double prior_c2[], double infop_c1[], double infop_c2[], double prior_v[], 
					 double width_b1,  double width_a1[],  double width_a2[],  double width_b2[],  double width_c[], double width_v,
					 int m, double init_b1, double *init_a1, double *init_a2, double *init_b2, int init_ncl, 
					 double **init_thetac, double init_v, int *init_freqcl, int *init_ci, long *idum, int nsave)  
{

	int i, j, k, k2, k3, npar=7+kG+2*kU;  
	int maxcl=n+m; // maximum number of clusters (<=n+m) (could use proper number to save memory);
	double *mu1, *sig1, *mu2, *sig2, *rho, *mu1_c, *sig1_c, *mu2_c, *sig2_c, *rho_c, *pclust;
	int *ci,  *freqcl, *dc;  
	double *Xc, *Tc, **Gc, **Uc; // to store data for cluster c;
	int ncl, single; // ncl=number of clusters along the subject updates; single=singleton indicator; 
	double b1, *a1, *a2, *b2, v, b1n, *a1n, *a2n, *b2n, vn, mu1_cn, sig1_cn, mu2_cn, sig2_cn, rho_cn;   //current and proposed values of parameters;
	double v_low, v_up, v_w; //parameters for priors of v;
	double logL, logLn, logA, u; //log-likelihoods; log acceptance prob; uniform(0,1)
	double w, mu_p, sig_p, g1, g2; //tentative proposal width and prior parameters for M-H sampling
	double info_shape_sig1, info_scale_sig1, info_shape_sig2, info_scale_sig2; // parameters of (slightly) informative priors for (sig1, sig2)
	int acp, *accept, s_clust; // acp=M-H accept indicator; s_clust=cumulative number of clusters
	double *b1r, *vr; 
	int *nclr;
	int i_tmp; // index for temp results saving
	FILE  *temp;  // to save results of b1 posterior samples and cluster indicator ci (at the end of cluster status update)

	mu1  = (double *) malloc(n*sizeof(double)); // store (mu1,sig1,mu2,sig2,rho) for each subject
	sig1 = (double *) malloc(n*sizeof(double)); 
	mu2  = (double *) malloc(n*sizeof(double)); 
	sig2 = (double *) malloc(n*sizeof(double)); 
	rho  = (double *) malloc(n*sizeof(double)); 
	mu1_c  = (double *) malloc(maxcl*sizeof(double)); // store (mu1,sig1,mu2,sig2,rho) for each cluster
	sig1_c = (double *) malloc(maxcl*sizeof(double));
	mu2_c  = (double *) malloc(maxcl*sizeof(double));
	sig2_c = (double *) malloc(maxcl*sizeof(double));
	rho_c  = (double *) malloc(maxcl*sizeof(double));
	ci = (int *) malloc(n*sizeof(int)); // cluster indicator
	a1 = (double *) malloc(kG*sizeof(double));
	a2 = (double *) malloc(kU*sizeof(double));
	b2 = (double *) malloc(kU*sizeof(double));
	a1n = (double *) malloc(kG*sizeof(double));
	a2n = (double *) malloc(kU*sizeof(double));
	b2n = (double *) malloc(kU*sizeof(double));
	freqcl = (int *) malloc(maxcl*sizeof(int)); // manually keep track of the frequencies of the clusters
	pclust = (double *) malloc(maxcl*sizeof(double)); // store probabiliteis of re-assigning a subject to the clusters
	accept = (int *) malloc(npar*sizeof(int)); // overall acceptance recorder
	// store the data for cluster k
	Xc = (double *) malloc(n*sizeof(double));
	Tc = (double *) malloc(n*sizeof(double)); 
	dc = (int *) malloc(n*sizeof(int)); 
	Gc = (double **) malloc(n*sizeof(double *));
	for (i=0; i<n; i++) Gc[i]=(double *) malloc(kG*sizeof(double));
	Uc = (double **) malloc(n*sizeof(double *));
	for (i=0; i<n; i++) Uc[i]=(double *) malloc(kU*sizeof(double));
	// for saving temp results (for running on Hoffman2 cluster)
	b1r = (double *) malloc(nsave*sizeof(double));
	vr = (double *) malloc(nsave*sizeof(double));
	nclr = (int *) malloc(nsave*sizeof(int));
	
	for (i=0; i<npar; i++) accept[i]=0;

	// initial values;
	b1=init_b1;
	for (i=0; i<kG; i++) { a1n[i] = a1[i] = init_a1[i]; }
	for (i=0; i<kU; i++) { a2n[i] = a2[i] = init_a2[i]; 
						   b2n[i] = b2[i] = init_b2[i]; }
	v=init_v;
	// initial clusters
	ncl=init_ncl;   // inincl is an integer to record the number of clusters along the subject updates
	for (i=0; i<n; i++) ci[i]=init_ci[i];
	s_clust=ncl;
	for (i=0; i<ncl; i++) {
		mu1_c[i]=init_thetac[i][0];
		sig1_c[i]=init_thetac[i][1];
		mu2_c[i]=init_thetac[i][2];
		sig2_c[i]=init_thetac[i][3];
		rho_c[i]=init_thetac[i][4];
		freqcl[i]=init_freqcl[i];
	}
	for (i=ncl; i<maxcl; i++) freqcl[i]=0;


	// parameters of (slightly) informative priors for (sig1, sig2)
	info_shape_sig1 = infop_c1[1] * infop_c1[1] / infop_c2[1] / infop_c2[1] + 2;  
	info_scale_sig1 = (info_shape_sig1 - 1) * infop_c1[1];
	info_shape_sig2 = infop_c1[3] * infop_c1[3] / infop_c2[3] / infop_c2[3] + 2; 
	info_scale_sig2 = (info_shape_sig2 - 1) * infop_c1[3];
 

temp=fopen("temp.txt", "a");
	fprintf(temp, "b1	v	ncl\n");
	fprintf(temp, "%f	%f	%d\n", b1, v, ncl);
	fclose(temp);

	i_tmp=0;

	/************** start MCMC sampling **************/

	for (j=1; j<=ch; j++) 
	{
//printf("%d\n", j); 

		// save individual values of (mu1_i, sig1_i, mu2_i, sig2_i, rho_i)
		indvalues5(mu1, sig1, mu2, sig2, rho,  n, ci, mu1_c, sig1_c, mu2_c, sig2_c, rho_c);

		// Log-likelihood;
		logL = logL_IV_DPM_all_cov(T, d, X, n, G, kG, U, kU, b1, a1, a2, b2, mu1, sig1, mu2, sig2, rho);


		/******** Metropolis-Hastings sampling for b1 ********/
		w=width_b1;	
		mu_p = prior_b1[0]; sig_p = prior_b1[1]; 
		b1n = runif(b1-w, b1+w);
		logLn = logL_IV_DPM_all_cov(T, d, X, n, G, kG, U, kU, b1n, a1, a2, b2, mu1, sig1, mu2, sig2, rho);
		logA = logLn - logL + ((b1-mu_p)*(b1-mu_p)-(b1n-mu_p)*(b1n-mu_p))/sig_p/2.0;
		if(logA>=0) acp=1; else { u = log((double) rand() / ((double)RAND_MAX+1)) ; if (logA>=u) acp=1; else acp=0; }
		if (acp==1) { 
			b1 = b1n;	logL = logLn;
			accept[0]++;
		}
		b1r[i_tmp] = b1;


		/******** Metropolis-Hastings sampling for a1 ********/
		for (k=0; k<kG; k++){
			w=width_a1[k];
			mu_p = prior_a1_1[k]; sig_p = prior_a1_2[k]; 
			a1n[k] = runif(a1[k]-w, a1[k]+w);
			logLn = logL_IV_DPM_all_cov(T, d, X, n, G, kG, U, kU, b1, a1n, a2, b2, mu1, sig1, mu2, sig2, rho);
			logA = logLn - logL + ((a1[k]-mu_p)*(a1[k]-mu_p)-(a1n[k]-mu_p)*(a1n[k]-mu_p))/sig_p/2.0;
			if(logA>=0) acp=1; else { u = log((double) rand() / ((double)RAND_MAX+1)) ; if (logA>=u) acp=1; else acp=0; }	
			if (acp==1) { 
				a1[k] = a1n[k];	logL = logLn;
				accept[1+k]++;
			}else{
				a1n[k]=a1[k];
			}
		}

		/******** Metropolis-Hastings sampling for a2 ********/
		for (k=0; k<kU; k++){
			w=width_a2[k];
			mu_p = prior_a2_1[k]; sig_p = prior_a2_2[k]; 
			a2n[k] = runif(a2[k]-w, a2[k]+w);
			logLn = logL_IV_DPM_all_cov(T, d, X, n, G, kG, U, kU, b1, a1, a2n, b2, mu1, sig1, mu2, sig2, rho);
			logA = logLn - logL + ((a2[k]-mu_p)*(a2[k]-mu_p)-(a2n[k]-mu_p)*(a2n[k]-mu_p))/sig_p/2.0;
			if(logA>=0) acp=1; else { u = log((double) rand() / ((double)RAND_MAX+1)) ; if (logA>=u) acp=1; else acp=0; }	
			if (acp==1) { 
				a2[k] = a2n[k];	logL = logLn;
				accept[1+kG+k]++;
			}else{
				a2n[k]=a2[k];
			}
		}	

		/******** Metropolis-Hastings sampling for b2 ********/
		for (k=0; k<kU; k++){
			w=width_b2[k];
			mu_p = prior_b2_1[k]; sig_p = prior_b2_2[k]; 
			b2n[k] = runif(b2[k]-w, b2[k]+w);
			logLn = logL_IV_DPM_all_cov(T, d, X, n, G, kG, U, kU, b1, a1, a2, b2n, mu1, sig1, mu2, sig2, rho);
			logA = logLn - logL + ((b2[k]-mu_p)*(b2[k]-mu_p)-(b2n[k]-mu_p)*(b2n[k]-mu_p))/sig_p/2.0;
			if(logA>=0) acp=1; else { u = log((double) rand() / ((double)RAND_MAX+1)) ; if (logA>=u) acp=1; else acp=0; }
			if (acp==1) { 
				b2[k] = b2n[k];	logL = logLn;		
				accept[1+kG+kU+k]++;
			}else{
				b2n[k]=b2[k];
			}
		}

		/******** update ci using Neal (2000) algorithm 8 ********/
		// Update ci
		for (i=0; i<n; i++)	// update the cluster status of each subject
		{
			freqcl[ci[i]]--;  
			if (freqcl[ci[i]]>0) single=0; else single=1;
			if (single==0)	// subject i is not a singleton
			{
				for (k=ncl; k<ncl+m; k++)  // generate m auxiliary parameters
				{
					mu1_c[k] = infop_c1[0] + gasdev(idum) * infop_c2[0];
					sig1_c[k] = 1.0 / GetGamma(info_shape_sig1, 1.0/info_scale_sig1, idum);
					mu2_c[k] = infop_c1[2] + gasdev(idum) * infop_c2[2];
					sig2_c[k] = 1.0 / GetGamma(info_shape_sig2, 1.0/info_scale_sig2, idum);
					rho_c[k] = runif(-1.0, 1.0);
				}
			}else  // subject i is a singleton
			{
				// relabel the singleton cluster as the last one
				mu1_c[ncl] = mu1_c[ci[i]];  // temporarily save as the (ncl+1)-th cluster
				sig1_c[ncl] = sig1_c[ci[i]]; 
				mu2_c[ncl] = mu2_c[ci[i]]; 				
				sig2_c[ncl] = sig2_c[ci[i]]; 
				rho_c[ncl] = rho_c[ci[i]]; 
				for (k=ci[i]; k<ncl; k++)  // relabel
				{
					mu1_c[k] = mu1_c[k+1]; 
					sig1_c[k] = sig1_c[k+1]; 
					mu2_c[k] = mu2_c[k+1]; 
					sig2_c[k] = sig2_c[k+1];
					rho_c[k] = rho_c[k+1];
					freqcl[k] = freqcl[k+1];
				}
				for (k=ncl; k<ncl+m-1; k++) // generate m-1 auxiliary parameters
				{
					mu1_c[k] = infop_c1[0] + gasdev(idum) * infop_c2[0];
					sig1_c[k] = 1.0 / GetGamma(info_shape_sig1, 1.0/info_scale_sig1, idum);
					mu2_c[k] = infop_c1[2] + gasdev(idum) * infop_c2[2];
					sig2_c[k] = 1.0 / GetGamma(info_shape_sig2, 1.0/info_scale_sig2, idum);
					rho_c[k] = runif(-1.0, 1.0);
				}
				for (k=0; k<n; k++) {if (ci[k]>ci[i]) ci[k]--;} // relabel ci accordingly
				//ci[i]=ncl-1; // not useful, about to re-assign anyway
			}
			ncl = ncl + m - single;

			// cluster re-assign probabilities
			for (k=0; k < (ncl-m); k++)  // ncl-m equals k-minus in Neal (2000)
				pclust[k] = freqcl[k] * L_IV_ind_cov(T[i], d[i], X[i], G[i], kG, U[i], kU, b1, a1, a2, b2, mu1_c[k], sig1_c[k], mu2_c[k], sig2_c[k], rho_c[k]) ;
			for (k=ncl-m; k<ncl; k++)
				pclust[k] = v/m * L_IV_ind_cov(T[i], d[i], X[i], G[i], kG, U[i], kU, b1, a1, a2, b2, mu1_c[k], sig1_c[k], mu2_c[k], sig2_c[k], rho_c[k]) ;

			// randomly re-assign subject i to a cluster
			ci[i]=sampleint(pclust, ncl);
			if (ci[i] < ncl-m)  ncl=ncl-m;
			else {
				ncl = ncl-m+1;
				mu1_c[ncl-1] = mu1_c[ci[i]];
				sig1_c[ncl-1] = sig1_c[ci[i]];
				mu2_c[ncl-1] = mu2_c[ci[i]];
				sig2_c[ncl-1] = sig2_c[ci[i]];
				rho_c[ncl-1] = rho_c[ci[i]];
				ci[i] = ncl-1;
			}
			freqcl[ci[i]]++;
		}
		s_clust += ncl;  
		nclr[i_tmp] = ncl;

		/******** update (mu1_c, sig1_c, mu2_c, sig2_c, rho_c) by Metropolis-Hastings sampling ********/
		for (k=0; k<ncl; k++)
		{
			k2 = 0;
			for (i=0; i<n; i++)
			{
				if (ci[i]==k){
					Xc[k2]=X[i];
					Tc[k2]=T[i];
					dc[k2]=d[i];
					for (k3=0; k3<kG; k3++) Gc[k2][k3]=G[i][k3];
					for (k3=0; k3<kU; k3++) Uc[k2][k3]=U[i][k3];
					k2++; // final k2 should = freqcl[k]
				}
			}
			//if(k2!=freqcl[k]) printf("NOT match\n");

			// Log-likelihood;
			logL = logL_IV_normal_cov(Tc, dc, Xc, k2, Gc, kG, Uc, kU, mu1_c[k], a1, a2, sig1_c[k], mu2_c[k], b1, b2, sig2_c[k], rho_c[k]); 

			// Metropolis-Hastings sampling for mu1_c 
			w=width_c[0];
			mu_p = prior_c1[0]; sig_p = prior_c2[0]; 
			mu1_cn = runif(mu1_c[k]-w, mu1_c[k]+w);
			logLn = logL_IV_normal_cov(Tc, dc, Xc, k2, Gc, kG, Uc, kU, mu1_cn, a1, a2, sig1_c[k], mu2_c[k], b1, b2, sig2_c[k], rho_c[k]); 
			logA = logLn - logL + ((mu1_c[k]-mu_p)*(mu1_c[k]-mu_p)-(mu1_cn-mu_p)*(mu1_cn-mu_p))/sig_p/2.0;
			if(logA>=0) acp=1; else {	u = log((double) rand() / ((double)RAND_MAX+1)) ; if (logA>=u) acp=1; else acp=0; }
			if (acp==1) {
				mu1_c[k] = mu1_cn;
				logL = logLn;
				accept[npar-6]++;
			}

			// Metropolis-Hastings sampling for sig1_c
			w=width_c[1];
			g1 = prior_c1[1]; g2 = prior_c2[1]; 
			sig1_cn = runif(max2(0.0, sig1_c[k]-w), sig1_c[k]+w);
			logLn = logL_IV_normal_cov(Tc, dc, Xc, k2, Gc, kG, Uc, kU, mu1_c[k], a1, a2, sig1_cn, mu2_c[k], b1, b2, sig2_c[k], rho_c[k]); 
			logA = logLn - logL + log(min2(2*w, sig1_c[k]+w)) - log(min2(2*w, sig1_cn+w)) + (g1+1)*(log(sig1_c[k])-log(sig1_cn)) + g2*(1/sig1_c[k]-1/sig1_cn);
			if(logA>=0) acp=1; else {	u = log((double) rand() / ((double)RAND_MAX+1)) ; if (logA>=u) acp=1; else acp=0; }
			if (acp==1) {
				sig1_c[k] = sig1_cn;
				logL = logLn;
				accept[npar-5]++;
			}

			// Metropolis-Hastings sampling for mu2_c 
			w=width_c[2];
			mu_p = prior_c1[2]; sig_p = prior_c2[2]; 
			mu2_cn = runif(mu2_c[k]-w, mu2_c[k]+w);
			logLn = logL_IV_normal_cov(Tc, dc, Xc, k2, Gc, kG, Uc, kU, mu1_c[k], a1, a2, sig1_c[k], mu2_cn, b1, b2, sig2_c[k], rho_c[k]); 
			logA = logLn - logL + ((mu2_c[k]-mu_p)*(mu2_c[k]-mu_p)-(mu2_cn-mu_p)*(mu2_cn-mu_p))/sig_p/2.0;
			if(logA>=0) acp=1; else {	u = log((double) rand() / ((double)RAND_MAX+1)) ; if (logA>=u) acp=1; else acp=0; }
			if (acp==1) {
				mu2_c[k] = mu2_cn;
				logL = logLn;
				accept[npar-4]++;
			}

			// Metropolis-Hastings sampling for sig2_c 
			w=width_c[3];
			g1 = prior_c1[3]; g2 = prior_c2[3]; 
			sig2_cn = runif(max2(0.0, sig2_c[k]-w), sig2_c[k]+w);
			logLn = logL_IV_normal_cov(Tc, dc, Xc, k2, Gc, kG, Uc, kU, mu1_c[k], a1, a2, sig1_c[k], mu2_c[k], b1, b2, sig2_cn, rho_c[k]); 
			logA = logLn - logL + log(min2(2*w, sig2_c[k]+w)) - log(min2(2*w, sig2_cn+w)) + (g1+1)*(log(sig2_c[k])-log(sig2_cn)) + g2*(1/sig2_c[k]-1/sig2_cn);
			if(logA>=0) acp=1; else {	u = log((double) rand() / ((double)RAND_MAX+1)) ; if (logA>=u) acp=1; else acp=0; }
			if (acp==1) {
				sig2_c[k] = sig2_cn;
				logL = logLn;
				accept[npar-3]++;
			}

			// Metropolis-Hastings sampling for rho_c
			w=width_c[4];
			rho_cn = runif(max2(-1.0, rho_c[k]-w),  min2(1.0, rho_c[k]+w));
			logLn = logL_IV_normal_cov(Tc, dc, Xc, k2, Gc, kG, Uc, kU, mu1_c[k], a1, a2, sig1_c[k], mu2_c[k], b1, b2, sig2_c[k], rho_cn); 
			logA = logLn - logL + log(min2(1.0, rho_c[k]+w)-max2(-1.0, rho_c[k]-w)) - log(min2(1.0, rho_cn+w)-max2(-1.0, rho_cn-w)) ;
			if(logA>=0) acp=1; else {	u = log((double) rand() / ((double)RAND_MAX+1)) ; if (logA>=u) acp=1; else acp=0; }
			if (acp==1) {
				rho_c[k] = rho_cn;
				logL = logLn;
				accept[npar-2]++;
			}
		}

		/******** update v by Metropolis-Hastings sampling ********/
		w=width_v;	
		v_low = prior_v[0]; v_up = prior_v[1]; v_w = prior_v[2];
		vn = runif(max2(v_low, v-w),  min2(v_up, v+w));
		logA = log(min2(v_up, v+w) - max2(v_low, v-w)) - log(min2(v_up, vn+w) - max2(v_low, vn-w)) + v_w * (log(v_up - vn) - log(v_up - v)) + ncl* (log(vn) - log(v)) + Ln_Gamma_Function(vn) - Ln_Gamma_Function(vn+n) - Ln_Gamma_Function(v) + Ln_Gamma_Function(v+n) ;
		if(logA>=0) acp=1; 
		else {	u = log((double) rand() / ((double)RAND_MAX+1)) ;
				if (logA>=u) acp=1; else acp=0; }
		if (acp==1) { 
			v = vn;
			accept[npar-1]++;
		}
		vr[i_tmp] = v;

		i_tmp++;


		if (i_tmp==nsave) {
		// save the temp results
temp=fopen("temp.txt", "a");
			for (i=0; i<nsave; i++)	fprintf(temp, "%f	%f	%d\n", b1r[i], vr[i], nclr[i]);  
			fclose(temp);
		// record last values for next round 
temp=fopen("temp2_last.txt", "w");
			fprintf(temp, "%f\n", b1);
			for(i=0; i<kG; i++) fprintf(temp, "%f	", a1[i]); fprintf(temp,"\n");
			for(i=0; i<kU; i++) fprintf(temp, "%f	", a2[i]); fprintf(temp,"\n");
			for(i=0; i<kU; i++) fprintf(temp, "%f	", b2[i]); fprintf(temp,"\n");
			fprintf(temp, "%f\n", v);
			fprintf(temp, "%d\n", ncl); 
			for(i=0; i<ncl; i++) fprintf(temp, "%f	%f	%f	%f	%f	%d\n", mu1_c[i], sig1_c[i], mu2_c[i], sig2_c[i], rho_c[i], freqcl[i]);		
			for(i=0; i<n; i++) fprintf(temp, "%d	", ci[i]); fprintf(temp,"\n");
			fclose(temp);
		// reset index to 0 to save another round of temp results;
			i_tmp = 0;	
		}
	}


/*	// acceptance rate (remove when run on hoffman2 cluster)
	printf("acceptance rate:\n");
	printf("b1: %f\n", (double)accept[0]/ch);
	printf("a1:"); for(i=0; i<kG; i++) printf("%f	", (double)accept[1+i]/ch); printf("\n");
	printf("a2:"); for(i=0; i<kU; i++) printf("%f	", (double)accept[1+kG+i]/ch); printf("\n");
	printf("b2:"); for(i=0; i<kU; i++) printf("%f	", (double)accept[1+kG+kU+i]/ch); printf("\n");
	printf("theta_c:"); for(i=0; i<5; i++) printf("%f	", (double)accept[1+kG+2*kU+i]/s_clust); printf("\n");
	printf("v: %f\n", (double)accept[6+kG+2*kU]/ch); 
	printf("average num of clust: %f\n", (double)s_clust/ch);
*/

	// free memories;
	free(ci); free(mu1); free(sig1); free(mu2); free(sig2); free(rho); 
	free(mu1_c); free(sig1_c); free(mu2_c); free(sig2_c); free(rho_c); 	
	free(freqcl); free(pclust);
	free(Xc); free(Tc); free(dc); free(Gc); free(Uc); 
	free(a1); free(a2); free(b2); free(a1n); free(a2n); free(b2n);

}
