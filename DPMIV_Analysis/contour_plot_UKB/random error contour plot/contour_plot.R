library(MASS)

temp = read.csv("random error contour plot/cluster status_end of chain.csv",header=T)
attach(temp)

cov1 = sqrt(sig1*sig2)*rho

i=1
Sigma = matrix(c(sig1[i],cov1[i],cov1[i],sig2[i]),2,2)
e_simul = mvrnorm(n=n[i], c(mu1[i],mu2[i]),  Sigma)

for (i in 2:nrow(temp)) {
e_simul = rbind(e_simul, mvrnorm(n=n[i], c(mu1[i],mu2[i]),  Sigma))
}

plot(e_simul, xlab='e1', ylab='e2')


i=1
Sigma = matrix(c(sig1[i],cov1[i],cov1[i],sig2[i]),2,2)
e_simul = mvrnorm(n=n[i]*20, c(mu1[i],mu2[i]),  Sigma)

for (i in 2:nrow(temp)) {
e_simul = rbind(e_simul, mvrnorm(n=n[i]*20, c(mu1[i],mu2[i]),  Sigma))
}

fit = kde2d(e_simul[,1], e_simul[,2], n=100)
filled.contour(fit, color = rainbow)
filled.contour(fit, color = rainbow,xlim=c(-1.5,0.5),ylim=c(9,18))


