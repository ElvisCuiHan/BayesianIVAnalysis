#include <math.h>
#define PI 3.14159265358979323846

double sum(double *a, int n)
{
	double s=0;
	int i;
	for(i=0;i<n; i++) s=s+a[i];
	return (s);
}

int sumint(int *a, int n)
{
	int s=0;
	int i;
	for(i=0;i<n; i++) s=s+a[i];
	return (s);
}


double mean(double *a, int n)
{
	double sum=0;
	int i;
	for(i=0;i<n; i++) sum=sum+a[i];
	sum=sum/n;
	return (sum);
}

/* sum of paired products of vector X and vector Y (dot product) */
double prodsum(double *X, double *Y, int n)
{
	double s=0;
	int i;
	for(i=0; i<n; i++) s=s+X[i]*Y[i];
	return (s);
}

/* mean of paired products of vector X and vector Y */
double prodmean(double *X, double *Y, int n)
{
	double s=0;
	int i;
	for(i=0; i<n; i++) s=s+X[i]*Y[i];
	s=s/n;
	return (s);
}


/* sum of Xi power m */
double arraypowsum(double *X, double m, int n)
{
	double s=0;
	int i;
	for(i=0; i<n; i++) s=s+pow(X[i], m);
	return (s);
}

/* mean of Xi power m */
double arraypowmean(double *X, double m, int n)
{
	double s=0;
	int i;
	for(i=0; i<n; i++) s=s+pow(X[i], m);
	s=s/n;
	return (s);
}

/* least square estimate of regression coefficient b1 in a simple linear regression Y~X */
double lm_beta1(double *Y, double *X, int n)
{
	double b1;
	b1 = (prodmean(X,Y,n) - mean(X,n)*mean(Y,n))/(arraypowmean(X,2,n)-pow(mean(X,n),2));
	return(b1);
}

/* intercept estimate b0, given b1 */
double lm_beta0(double *Y, double *X, int n, double b1)
{
	double b0;
	b0 = mean(Y,n)-mean(X,n)*b1;
	return(b0);
}

/* sample variance of vector X */
double var(double *X, int n)
{
	double v;
	v=arraypowsum(X,2,n)/(n-1)-pow(mean(X,n),2)*n/(n-1);
	return(v);
}


/* estimate of error variance */
double lm_sig2(double *Y, double *X, int n, double b0, double b1)
{
	double sig2, *eps;
	int i;
	eps = (double *) malloc(n*sizeof(double));
	for(i=0; i<n; i++) eps[i]=Y[i]-b0-b1*X[i];
	sig2 = arraypowsum(eps,2,n)/(n-2);
	return(sig2);
}


/* SE of b1 estimate */
double lm_b1SE(double *Y, double *X, int n, double b1)
{
	double Syy, Sy, Sxx, Sx, Se2, Sb1;
	Syy=arraypowsum(Y,2,n);
	Sxx=arraypowsum(X,2,n);
	Sy=sum(Y,n);
	Sx=sum(X,n);
	Se2=(n*Syy-pow(Sy,2)-pow(b1,2)*(n*Sxx-pow(Sx,2)))/n/(n-2);
	Sb1=sqrt(n*Se2/(n*Sxx-pow(Sx,2)));
	return(Sb1);
}

/* SE of b0 estimate */
double lm_b0SE(double *X, int n, double Sb1)
{
	double Sxx, Sb0;
	Sxx=arraypowsum(X,2,n);
	Sb0=Sb1*sqrt(Sxx/n);
	return(Sb0);
}

/* absolute value of double */
double absD(double x)
{
    if (x < 0.0)
        x = -x;
    return x;
}

/* compare function for sorting */
int compare_doubles (const void *a, const void *b)
{
	const double *da = (const double *) a;
    const double *db = (const double *) b;
    return (*da > *db) - (*da < *db);
}


/* generate a random number from 'sam' with probabilities proportional to 'prob' with length n */
double sample(double *sam, double *prob, int n)
{
	int x=0, i, found=0;
	double s=sum(prob, n), u, cump=0;
	
	for (i=0; i<n; i++) prob[i]=prob[i]/s; //standardize 'prob' to probabilities;
	
	u = ((double) rand() / ((double)RAND_MAX+1)) ; //a random uniform [0,1);

	do {
		cump=cump+prob[x];
		if (u<=cump) found=1;
		x=x+1;
	} while (found<1 && x<n);

	return(sam[x-1]);
}


/* generate a random integer from 0 to n-1 with probabilities proportional to 'prob' with length n */
int sampleint(double *prob, int n)
{
	int x=0, i, found=0;
	double s=sum(prob, n), u, cump=0;
	
	for (i=0; i<n; i++) prob[i]=prob[i]/s; //standardize 'prob' to probabilities;
	
	u = ((double) rand() / ((double)RAND_MAX+1)) ; //a random uniform [0,1);

	do {
		cump=cump+prob[x];
		if (u<=cump) found=1;
		x=x+1;
	} while (found<1 && x<n);

	return(x-1);
}


/* see if the jth element in a vector of cluster indicator ci (length n) is a singleton */
int singleton(int *ci, int j, int n)
{
	int i=0, single=0;
	do{
		if (ci[j]==ci[i] && j!=i) single=1;
		i++;	
	} while (single<1 && i<n);

	return(1-single);
}

/* pdf of normal distribution */
double dnorm(double x, double mu, double sig2)
{
	double d;
	d=exp(-(x-mu)*(x-mu)/(sig2*2))/sqrt(2.0*PI*sig2);
	return(d);
}


/* generate a random number by uniform distribution (lower, upper) */
double runif(double lower, double upper)
{
	double	u;
	u = (((double)rand())+1.0) / (((double)RAND_MAX)+2.0)  ; //a random uniform (0,1);
	u = lower + (upper-lower)*u; 
	return(u);
}



/* CDF function of standard normal distribution */
double pnorm(double d)
{
    double       A1 = 0.31938153;
    double       A2 = -0.356563782;
    double       A3 = 1.781477937;
    double       A4 = -1.821255978;
    double       A5 = 1.330274429;
    double RSQRT2PI = 0.39894228040143267793994605993438;

    double
    K = 1.0 / (1.0 + 0.2316419 * fabs(d));

    double
    cnd = RSQRT2PI * exp(- 0.5 * d * d) *
          (K * (A1 + K * (A2 + K * (A3 + K * (A4 + K * A5)))));

    if (d > 0)
        cnd = 1.0 - cnd;

    return cnd;
}


/* maximum and minimum of two double variables */
double max2(double x, double y)
{
	double z;
	if (x>=y) z=x; else z=y;
	return(z);
}
double min2(double x, double y)
{
	double z;
	if (x<=y) z=x; else z=y;
	return(z);
}