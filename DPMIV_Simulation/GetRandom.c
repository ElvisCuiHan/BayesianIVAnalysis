#include <math.h>


double GetGamma(double shape, double scale, long *idum)
{
	// use algorithm from SimpleRNG;
	double d, c, x, xsquared, v, u;
    if (shape >= 1.0)
    {
		d = shape - 1.0/3.0;
		c = 1.0/sqrt(9.0*d);
		for (;;)
		{
			do
			{
				x = gasdev(idum);
				v = 1.0 + c*x;
			}
			while (v <= 0.0);
			v = v*v*v;
			u = ((double) rand() / ((double)RAND_MAX+1)) ;
			xsquared = x*x;
			if (u < 1.0 -.0331*xsquared*xsquared || log(u) < 0.5*xsquared + d*(1.0 - v + log(v)))
				return scale*d*v;
		}
	}
	else
	{
		double g = GetGamma(shape+1.0, 1.0, idum);
		double w = ((double) rand() / ((double)RAND_MAX+1)) ;
		return scale*g*pow(w, 1.0/shape);
	}
}


double GetBeta(double a, double b, long *idum)
{
	double	u = GetGamma(a, 1.0,  idum);
	double  v = GetGamma(b, 1.0,  idum);
	return u/(u+v);
}