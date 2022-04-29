# Import IV analysis function; add route if necessary
source("IV_MH.R")

# simulation settings
n=500   # sample size
G1=rnorm(n) # instrument 1
G2=rnorm(n) # instrument 2
U_obs1=rnorm(n) # observed confounder 1
U_obs2=rnorm(n) # observed confounder 2
U_unobs=rnorm(n) # unobserved confounder
e1=rnorm(n, sd=1)# random error in intermediate covariate
a0=0; a1=c(0.5, 0.5); a2=c(0.5, 0.5, 0.5); # parameters of the first-stage model
W = a0 + cbind(G1,G2)%*%a1 + cbind(U_obs1, U_obs2, U_unobs)%*%a2 + e1;  # intermediate covariate
e2=rnorm(n, sd=0.5) # measurement error
X = W + e2;   # noisy suggorate X
b0=2; b1=1; b2=c(1, 1, 1); # parameters of the second-stage model; b1 is parameter of primary interest
e3=rnorm(n, sd=2) # random error in time-to-event
T = b0 + b1*W + cbind(U_obs1, U_obs2, U_unobs)%*%b2 + e3 # time-to-event
eC=rnorm(n, sd=2) # random error in censoring time
C = b0 + b1*W + cbind(U_obs1, U_obs2, U_unobs)%*%b2 + eC # censoring time
Y = pmin(T,C)  # observed time
d = as.numeric(T<=C) # censoring indicator
# Note: The observed data consists of: Y, d, X, U_obs1, U_obs2, G1, and G2

# IV model
fit = IV_MH(
Y=Y, 
d=d, 
X=X, 
G=cbind(G1,G2), 
U=cbind(U_obs1, U_obs2), 
m=11000, 
# adjust width to give reasonable acceptance rate
wid=c(0.25,0.25,0.25,0.25,0.25,0.4,0.5,0.4,0.6,0.6,1.7,0.2), 
# priors: N(0,10000) for a0,a1,a2,b0,b1; Inv-Gamma(0.001, 0.001) for sig1 and sig2
prior_1=c(rep(0, 5), .001, rep(0,4), .001),
prior_2=c(rep(10000,5), .001, rep(10000,4), .001)
) 

# acceptance rate for (a0,a1,a2,sig1,b0,b1,b2,sig2,rho)
fit$ap
# trace plot
plot(fit$b1r, type='l')
# posterior samples of b1 with 1000 burn-in
b1r_res=fit$b1r[1001:11000]
# histogram
hist(b1r_res)
# posterior means
mean(b1r_res)
# 95% credible interval
quantile(b1r_res,c(0.025,0.975))

