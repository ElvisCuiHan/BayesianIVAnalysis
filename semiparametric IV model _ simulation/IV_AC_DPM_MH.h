#include <math.h>
#include "gamma_function.c"
#include "ln_gamma_function.c"

/********* Function for IV linear regression for arbitrary-censored survival outcome with Dirichlet Process Mixtures *********/
/*  Revision date: 2014-02-16
	Use Algorithm 8 in Neal (2000)
	Use prior for strength parameter in DPM proposed by Conley et al. (2008)
	Note: This version only saves the posterior mean, SD and 95% CI for b1 (causal effect of true covariate W on outcome Y)
	and the posterior mean, SD and 95% CI for v (base/strength parameter of Dirichlet prior)
	X = a vector of univariate covariate (that might be subject to measurement error and/or unobserved confounding)
	L = a vector of observed left-censored time
	R = a vector of observed right-censored time
	d = a vector of censoring indicator (1=left-censored; 2=interval-censored; 3=right-censored; 4=event)
	G = a vector of instrumental variable (univariate)
	n = sample size (length of X, L, R, d and G)
	ch = number of iterations (length of the chain)
	burn = number of burn-ins
  	thin = thin the chain by taking every 'thin'-th sample
	prior_a1 = a vector of the mean and variance for the normal prior for a1, e.g. {0, 10000}
	prior_b1 = a vector of the mean and variance for the normal prior for b1, e.g. {0, 10000}
	note: the sig1 and sig2 in this program refer to sigma_1-squared and sigma_2-squared in the paper
	prior_c1 = a vector of the first parameter of the priors for (mu1c, sig1c, mu2c, sig2c), e.g. {0, 0.001, 0, 0.001}
			mean of the normal prior for mu1c and mu2c; shape parameter of the inverse-gamma prior for sig1c and sig2c;
	prior_c2 = a vector of the second parameter of the priors for (mu1c, sig1c, mu2c, sig2c), e.g. {10000, 0.001, 10000, 0.001}
			variance of the normal prior for mu1c and mu2c; scale parameter of the inverse-gamma prior for sig1c and sig2c; 
	infop_c1 = a vector of the first parameter of the (slightly) informative priors for (mu1c, sig1c, mu2c, sig2c), e.g. {a0t, sig1t, b0t, sig2t}
			mean of the normal prior for mu1c and mu2c; mean of the inverse-gamma prior for sig1c and sig2c
	infop_c2 = a vector of the second parameter of the (slightly) informative priors for (mu1c, sig1c, mu2c, sig2c), e.g. {5, 3, 10, 5}
			SD of the normal prior for mu1c and mu2c; SD of the inverse-gamma prior for sig1c and sig2c
	prior_v = a vector of the parameters of the priors for v (the base/strength parameter of Dirichlet prior)
			(v_low, v_up, v_w):  prior of v is proportional to ((v_up-v)/(v_up-v_low))^v_w  (Conley et al. 2008)
    width = a vector of the width of random walk for (a1,b1, mu1c,sig1c,mu2c,sig2c,rhoc, v)
	m = number of auxiliary parameters in Dirichlet process sampling, e.g. 10
	res_b1 = a vector to store results for b1 (length>=6)
	res_v = a vector to store results for v and nclust (length>=7)
	res_acp = a vector to store acceptance rate for (a1,b1,mu1c,sig1c,mu2c,sig2c,rhoc,v) and average number of clusters (length>=9)
	idum = seed number for random normal sample generator
	least-square estimate (a1,b1) or its 95% CI (mu1c,sig1c,mu2c,sig2c, with 2 initial clusters) will be used as initial values
*/


void IV_AC_DPM_MH(double *X, double *L, double *R, int *d, double *G, int n, int ch, int burn, int thin, double prior_a1[], double prior_b1[], 
		   double prior_c1[], double prior_c2[], double infop_c1[], double infop_c2[], double prior_v[], double width[], int m, 
		   double res_b1[], double res_v[], double res_acp[], long *idum)  
{
	int i, j, k, k2, acp, accept[]={0,0,0,0,0,0,0,0}; // acp=M-H accept indicator; accept=overall acceptance recorder; 
	int maxcl=n+m; // maximum number of clusters (<=n+m) (could use proper number to save memory);
	// for pointers, see 'malloc' part for annotations
	double *a1r, *b1r, *vr, *a1res, *b1res, *vres, *nclres;
	double *mu1, *sig1, *mu2, *sig2, *rho, *mu1_c, *sig1_c, *mu2_c, *sig2_c, *rho_c, *pclust;
	int *ci, *nclust, *freqcl, *dc;
	double *Xc, *Lc, *Rc, *Gc; // to store data for cluster c;
	int ncl, single, nclusttot; // ncl=number of clusters along the subject updates; single=singleton indicator; nclusttot=cumulative number of clusters (used to calculate overall acceptance); 
	double a1, b1, v, a1n, b1n, vn, mu1_cn, sig1_cn, mu2_cn, sig2_cn, rho_cn;   //current and proposed values of parameters;
	double v_low, v_up, v_w; //parameters for priors of v;
	double logL, logLn, logA, u; //log-likelihoods; log acceptance prob; uniform(0,1)
	double w, mu_p, sig_p, g1, g2; //tentative proposal width and prior parameters for M-H sampling
	int chfin = (ch-burn)/thin;   //final length of chain after burn-in and thinning
	double quar[3], Stmp; //quar=quartiles for b1; tentative variable 
	double a0, sig1h, a1se, a0se, b0, sig2h, b1se, b0se; // only used to calculate initial values
	double info_shape_sig1, info_scale_sig1, info_shape_sig2, info_scale_sig2; // parameters of (slightly) informative priors for (sig1, sig2)
	//FILE  *temp;  // to save results of b1 posterior samples and cluster indicator ci (at certain step). To use this, un-comment places about 'temp' below

	a1r = (double *) malloc(ch*sizeof(double));   // save results of a1 before burn-in and thinning
	b1r = (double *) malloc(ch*sizeof(double));   // save results of b1 before burn-in and thinning
	vr = (double *) malloc(ch*sizeof(double));   // save results of v before burn-in and thinning
	a1res = (double *) malloc(chfin*sizeof(double));	// final results of a1
	b1res = (double *) malloc(chfin*sizeof(double));	// final results of b1
	vres = (double *) malloc(chfin*sizeof(double));		// final results of v
	nclres = (double *) malloc(chfin*sizeof(double));  // final results of number of clusters
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
	ci = (int *) malloc(n*sizeof(int)); // cluster indicator for the n subjects
	nclust = (int *) malloc(ch*sizeof(int)); // store number of clusters (at the end of each update)
	freqcl = (int *) malloc(maxcl*sizeof(int)); // manually keep track of the frequencies of the clusters
	pclust = (double *) malloc(maxcl*sizeof(double)); // store probabiliteis of re-assigning a subject to the clusters
	// store the data for cluster k
	Xc = (double *) malloc(n*sizeof(double));
	Lc = (double *) malloc(n*sizeof(double)); 
	Rc = (double *) malloc(n*sizeof(double)); 
	dc = (int *) malloc(n*sizeof(int)); 
	Gc = (double *) malloc(n*sizeof(double)); 

	// use least-square estimates as initial values;
	a1 = a1r[0] = lm_beta1(X, G, n);  a0 = lm_beta0(X, G, n, a1); 	sig1h = lm_sig2(X, G, n, a0, a1);
	a1se = lm_b1SE(X, G, n, a1); 	a0se = lm_b0SE(G, n, a1se);
	b1 = b1r[0] = lm_beta1(R, X, n);  b0 = lm_beta0(R, X, n, b1);	sig2h = lm_sig2(R, X, n, b0, b1);
	b1se = lm_b1SE(R, X, n, b1); 	b0se = lm_b0SE(X, n, b1se);
	// two (arbitrary) initial clusters
	for (i=0; i<(int)n/2; i++) ci[i]=0; for (i=(int)n/2; i<n; i++) ci[i]=1;
	ncl=nclust[0]=2;   // ncl is an integer to record the number of clusters along the subject updates
	freqcl[0]=(int)n/2; freqcl[1]=n-(int)n/2;  	// frequencies of the two clusters
	for (i=2; i<n+m; i++) freqcl[i]=0;
	// bounds of 95% CI as initial values for mu_1, mu_2
	mu1_c[0] = a0 - 1.96*a0se; 	mu1_c[1] = a0 + 1.96*a0se;
	mu2_c[0] = b0 - 1.96*b0se; 	mu2_c[1] = b0 + 1.96*b0se;
	// .5 times and 2 times least-square estimate as initial values for sig1, sig2
	sig1_c[0] = sig1h*0.5;  sig1_c[1] = sig1h*2.0;
	sig2_c[0] = sig2h*0.5;  sig2_c[1] = sig2h*2.0;
	// arbitrary initial values for rho
	rho_c[0] = 0; rho_c[1] = 0;
	// arbitrary initial value for v
	v = vr[0] = 1.0; 

	// parameters of (slightly) informative priors for (sig1, sig2)
	info_shape_sig1 = infop_c1[1] * infop_c1[1] / infop_c2[1] / infop_c2[1] + 2;  
	info_scale_sig1 = (info_shape_sig1 - 1) * infop_c1[1];
	info_shape_sig2 = infop_c1[3] * infop_c1[3] / infop_c2[3] / infop_c2[3] + 2; 
	info_scale_sig2 = (info_shape_sig2 - 1) * infop_c1[3];


	/************** start MCMC sampling **************/

	for (j=1; j<ch; j++) 
	{
		// save individual values of (mu1_i, sig1_i, mu2_i, sig2_i, rho_i)
		indvalues5(mu1, sig1, mu2, sig2, rho,  n, ci, mu1_c, sig1_c, mu2_c, sig2_c, rho_c);

		// Log-likelihood;
		logL = logL_IV_AC_DPM_all(L, R, d, X, G, n, a1, b1, mu1, sig1, mu2, sig2, rho);

		/******** Metropolis-Hastings sampling for a1 ********/
		w=width[0];	mu_p = prior_a1[0]; sig_p = prior_a1[1]; 
		a1n = runif(a1-w, a1+w);
		logLn = logL_IV_AC_DPM_all(L, R, d, X, G, n, a1n, b1, mu1, sig1, mu2, sig2, rho);
		logA = logLn - logL + ((a1-mu_p)*(a1-mu_p)-(a1n-mu_p)*(a1n-mu_p))/sig_p/2.0;
		if(logA>=0) acp=1; 
		else {	u = log((double) rand() / ((double)RAND_MAX+1)) ;
				if (logA>=u) acp=1; else acp=0; }
		if (acp==1) { 
			a1 = a1r[j] = a1n;
			logL = logLn;
			accept[0]++;}
		else{ a1r[j] = a1; }

		/******** Metropolis-Hastings sampling for b1 ********/
		w=width[1];	mu_p = prior_b1[0]; sig_p = prior_b1[1]; 
		b1n = runif(b1-w, b1+w);
		logLn = logL_IV_AC_DPM_all(L, R, d, X, G, n, a1, b1n, mu1, sig1, mu2, sig2, rho);
		logA = logLn - logL + ((b1-mu_p)*(b1-mu_p)-(b1n-mu_p)*(b1n-mu_p))/sig_p/2.0;
		if(logA>=0) acp=1; 
		else {	u = log((double) rand() / ((double)RAND_MAX+1)) ;
				if (logA>=u) acp=1; else acp=0; }
		if (acp==1) { 
			b1 = b1r[j] = b1n;
			logL = logLn;
			accept[1]++;}
		else{ b1r[j] = b1; }


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
				pclust[k] = freqcl[k] * L_IV_AC_ind(L[i], R[i], d[i], X[i], G[i], a1, b1, mu1_c[k], sig1_c[k], mu2_c[k], sig2_c[k], rho_c[k]) ; 
			for (k=ncl-m; k<ncl; k++)
				pclust[k] = v/m * L_IV_AC_ind(L[i], R[i], d[i], X[i], G[i], a1, b1, mu1_c[k], sig1_c[k], mu2_c[k], sig2_c[k], rho_c[k]) ; 
		
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
		nclust[j]=ncl;  // record number of clusters after updating all subjects


		/******** update (mu1_c, sig1_c, mu2_c, sig2_c, rho_c) by Metropolis-Hastings sampling ********/
		for (k=0; k<ncl; k++)
		{
			k2 = 0;
			for (i=0; i<n; i++)
			{
				if (ci[i]==k){
					Xc[k2]=X[i];
					Lc[k2]=L[i];
					Rc[k2]=R[i];
					dc[k2]=d[i];
					Gc[k2]=G[i];
					k2++; // final k2 should = freqcl[k]
				}
			}
			//if(k2!=freqcl[k]) printf("NOT match\n");


			// Log-likelihood;
			logL = logL_IV_AC_normal(Lc, Rc, dc, Xc, Gc, k2, mu1_c[k], a1, sig1_c[k], mu2_c[k], b1, sig2_c[k], rho_c[k]);

			// Metropolis-Hastings sampling for mu1_c 
			w=width[2];
			mu_p = prior_c1[0]; sig_p = prior_c2[0]; 
			mu1_cn = runif(mu1_c[k]-w, mu1_c[k]+w);
			logLn = logL_IV_AC_normal(Lc, Rc, dc, Xc, Gc, k2, mu1_cn, a1, sig1_c[k], mu2_c[k], b1, sig2_c[k], rho_c[k]);
			logA = logLn - logL + ((mu1_c[k]-mu_p)*(mu1_c[k]-mu_p)-(mu1_cn-mu_p)*(mu1_cn-mu_p))/sig_p/2.0;
			if(logA>=0) acp=1; 
			else {	u = log((double) rand() / ((double)RAND_MAX+1)) ;
					if (logA>=u) acp=1; else acp=0; }
			if (acp==1) {
				mu1_c[k] = mu1_cn;
				logL = logLn;
				accept[2]++;
			}

			// Metropolis-Hastings sampling for sig1_c
			w=width[3];
			g1 = prior_c1[1]; g2 = prior_c2[1]; 
			sig1_cn = runif(max2(0.0, sig1_c[k]-w), sig1_c[k]+w);
			logLn = logL_IV_AC_normal(Lc, Rc, dc, Xc, Gc, k2, mu1_c[k], a1, sig1_cn, mu2_c[k], b1, sig2_c[k], rho_c[k]);
			logA = logLn - logL + log(min2(2*w, sig1_c[k]+w)) - log(min2(2*w, sig1_cn+w)) + (g1+1)*(log(sig1_c[k])-log(sig1_cn)) + g2*(1/sig1_c[k]-1/sig1_cn);
			if(logA>=0) acp=1; 
			else {	u = log((double) rand() / ((double)RAND_MAX+1)) ;
					if (logA>=u) acp=1; else acp=0; }
			if (acp==1) {
				sig1_c[k] = sig1_cn;
				logL = logLn;
				accept[3]++;
			}

			// Metropolis-Hastings sampling for mu2_c 
			w=width[4];
			mu_p = prior_c1[2]; sig_p = prior_c2[2]; 
			mu2_cn = runif(mu2_c[k]-w, mu2_c[k]+w);
			logLn = logL_IV_AC_normal(Lc, Rc, dc, Xc, Gc, k2, mu1_c[k], a1, sig1_c[k], mu2_cn, b1, sig2_c[k], rho_c[k]);
			logA = logLn - logL + ((mu2_c[k]-mu_p)*(mu2_c[k]-mu_p)-(mu2_cn-mu_p)*(mu2_cn-mu_p))/sig_p/2.0;
			if(logA>=0) acp=1; 
			else {	u = log((double) rand() / ((double)RAND_MAX+1)) ;
					if (logA>=u) acp=1; else acp=0; }
			if (acp==1) {
				mu2_c[k] = mu2_cn;
				logL = logLn;
				accept[4]++;
			}

			// Metropolis-Hastings sampling for sig2_c 
			w=width[5];
			g1 = prior_c1[3]; g2 = prior_c2[3]; 
			sig2_cn = runif(max2(0.0, sig2_c[k]-w), sig2_c[k]+w);
			logLn = logL_IV_AC_normal(Lc, Rc, dc, Xc, Gc, k2, mu1_c[k], a1, sig1_c[k], mu2_c[k], b1, sig2_cn, rho_c[k]);
			logA = logLn - logL + log(min2(2*w, sig2_c[k]+w)) - log(min2(2*w, sig2_cn+w)) + (g1+1)*(log(sig2_c[k])-log(sig2_cn)) + g2*(1/sig2_c[k]-1/sig2_cn);
			if(logA>=0) acp=1; 
			else {	u = log((double) rand() / ((double)RAND_MAX+1)) ;
					if (logA>=u) acp=1; else acp=0; }
			if (acp==1) {
				sig2_c[k] = sig2_cn;
				logL = logLn;
				accept[5]++;
			}

			// Metropolis-Hastings sampling for rho_c
			w=width[6];
			rho_cn = runif(max2(-1.0, rho_c[k]-w),  min2(1.0, rho_c[k]+w));
			logLn = logL_IV_AC_normal(Lc, Rc, dc, Xc, Gc, k2, mu1_c[k], a1, sig1_c[k], mu2_c[k], b1, sig2_c[k], rho_cn);
			logA = logLn - logL + log(min2(1.0, rho_c[k]+w)-max2(-1.0, rho_c[k]-w)) - log(min2(1.0, rho_cn+w)-max2(-1.0, rho_cn-w)) ;
			if(logA>=0) acp=1; 
			else {	u = log((double) rand() / ((double)RAND_MAX+1)) ;
					if (logA>=u) acp=1; else acp=0; }
			if (acp==1) {
				rho_c[k] = rho_cn;
				logL = logLn;
				accept[6]++;
			}
		}

		/******** update v by Metropolis-Hastings sampling ********/
		w=width[7];	
		v_low = prior_v[0]; v_up = prior_v[1]; v_w = prior_v[2];
		vn = runif(max2(v_low, v-w),  min2(v_up, v+w));
		logA = log(min2(v_up, v+w) - max2(v_low, v-w)) - log(min2(v_up, vn+w) - max2(v_low, vn-w)) + v_w * (log(v_up - vn) - log(v_up - v)) + ncl* (log(vn) - log(v)) + Ln_Gamma_Function(vn) - Ln_Gamma_Function(vn+n) - Ln_Gamma_Function(v) + Ln_Gamma_Function(v+n) ;
		if(logA>=0) acp=1; 
		else {	u = log((double) rand() / ((double)RAND_MAX+1)) ;
				if (logA>=u) acp=1; else acp=0; }
		if (acp==1) { 
			v = vr[j] = vn;
			accept[7]++;}
		else{ vr[j] = v; }

	//printf("%d %f %f %f\n", j, a1r[j], b1r[j], vr[j]);  //no printing when submitted to cluster
	
	}

	/*
	temp=fopen("temp.txt", "a");
	fprintf(temp, "a1	b1	v	nclust\n");
	for (i=0; i<ch; i++)
	fprintf(temp, "%f	%f	%f	%d\n", a1r[i], b1r[i], vr[i], nclust[i]);
	fclose(temp);

	temp=fopen("temp2.txt", "a");
	fprintf(temp, "ci\n");
	for (i=0; i<n;i++) 	fprintf(temp, "%d\n", ci[i]);
	fclose(temp); 
	*/


	for (i=0; i<chfin; i++) {
		b1res[i] = b1r[burn+(i+1)*thin-1];
		a1res[i] = a1r[burn+(i+1)*thin-1];
		vres[i] = vr[burn+(i+1)*thin-1];
		nclres[i] = (double)nclust[burn+(i+1)*thin-1];
	}
	
	res_b1[0]=mean(b1res, chfin); 
	res_b1[2]=sqrt(var(b1res, chfin));
 	qsort(b1res, chfin, sizeof(b1res[0]), compare_doubles);
	quartiles(b1res, chfin, quar);
	res_b1[1]=quar[1]; // median
	res_b1[3]=quar[0]; // lower bound of 95% CI
	res_b1[4]=quar[2]; // upper bound of 95% CI

	res_v[0]=mean(vres, chfin); 
	res_v[2]=sqrt(var(vres, chfin));
 	qsort(vres, chfin, sizeof(vres[0]), compare_doubles);
	quartiles(vres, chfin, quar);
	res_v[1]=quar[1]; // median
	res_v[3]=quar[0]; // lower bound of 95% CI
	res_v[4]=quar[2]; // upper bound of 95% CI
	res_v[5]=mean(nclres, chfin); 	// average number of clusters
	res_v[6]=sqrt(var(nclres, chfin));

	// acceptance rate
	for (i=0; i<2; i++) res_acp[i]=((double)accept[i])/ch;
	nclusttot=sumint(nclust,ch);
	for (i=2; i<7; i++) res_acp[i]=((double)accept[i])/nclusttot;
	res_acp[7]=((double)accept[7])/ch;
	// averege number of clusters
	Stmp=0; for (i=burn+thin-1; i<ch; i=i+thin) Stmp+=nclust[i];
	res_acp[8] =Stmp/chfin;  

	// free memories;
	free(a1r); free(b1r); free(a1res); free(b1res); free(vr); free(vres);
	free(ci); free(nclust); 
	free(mu1); free(sig1); free(mu2); free(sig2); free(rho); 
	free(mu1_c); free(sig1_c); free(mu2_c); free(sig2_c); free(rho_c); 	
	free(freqcl); free(pclust);
	free(Xc); free(Lc); free(Rc); free(dc); free(Gc); 

}
