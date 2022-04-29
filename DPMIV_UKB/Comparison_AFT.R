rm(list=ls())
library(icenReg)
library(tidyverse)
library(bayesSurv)
library(coda)

#setwd("~/Downloads/Dr Li GSR/UKB Data/Interval_Censoring_Project_new/")
setwd("~/Desktop/coding/DPMIV_UKB-0323/")

filename <- "0322_DPMIV_UKB_less_3_6.dat"

#dm_filtered = read.delim(file = "DPMIV_test_data_eth1.dat", sep=" ", header=F)
#dm_filtered = read.delim(file = "0321_DPMIV_UKB_5U.dat", sep=" ", header=T)
dm_filtered = read.delim(file = filename, sep=" ", header=T)
ca = read_delim(file="0321_DPMIV_UKB_5U.dat")
colnames(dm_filtered) = colnames(ca)[-23]
dm_filtered$R[ dm_filtered$d == 3] = Inf

IND = matrix(NA, ncol=1, nrow=dim(dm_filtered)[1])
IND[dm_filtered$d == 2] = 3
IND[dm_filtered$d == 4] = 1
IND[dm_filtered$d == 3] = 0
IND[dm_filtered$d == 1] = 2
dm_filtered$d = IND

#dm_filtered = as.data.frame(dm_filtered)[1:300,]

################### Let's run some AFT models! ####################

### Frequentist AFT Model

mod1 = ic_par(cbind(exp(L), exp(R)) ~ First, model="aft", 
            data = dm_filtered, dist = "lnorm")
summary(mod1)

mod2 = ic_par(cbind(exp(L), exp(R)) ~ First  + SEX + AGE_RECRUITMENT +
                CHOLESTEROL + BMI + SMOKE, model="aft", 
              data = dm_filtered, dist = "lnorm")
summary(mod2)

mod3 = ic_par(cbind(exp(L), exp(R)) ~ ., data=dm_filtered[,-3])
summary(mod3)

### Two Stage AFT Model
# Noisy surrogate
X = as.matrix(dm_filtered[, c(4)])
# Instruments
G = as.matrix(dm_filtered[, c(5)])
U = as.matrix(dm_filtered[, c(6:10)])
fit1=lm(X~G+U, na.action=na.exclude);
pred=predict(fit1);
L = dm_filtered[, 1]
R = dm_filtered[, 2]
fit2=ic_par(cbind(exp(L), exp(R)) ~ pred + U, model="aft", 
            data = dm_filtered, dist = "lnorm")
summary(fit2)

### Semiparametric Bayesian AFT model

mod_null = survreg(Surv(time=exp(L), time2=exp(R), type="interval2")~First+ SEX + AGE_RECRUITMENT +
                     CHOLESTEROL + BMI + SMOKE, 
        data=dm_filtered, dist="lognormal")
summary(mod_null)

X <- bayessurvreg1(Surv(time=exp(L), time2=exp(R), type="interval2") ~ First + SEX + AGE_RECRUITMENT +
                      CHOLESTEROL + BMI + SMOKE, 
                  data = dm_filtered, onlyX = TRUE)
nregres <- dim(X)[2]

# Priors for the mixture
prior <- list()
prior$k.prior = "uniform"
prior$kmax <- 30
prior$k.prior <- "poisson"
prior$poisson.k <- 5
prior$dirichlet.w <- 1
prior$mean.mu <- 3.07
prior$var.mu <- 5.52^2
prior$shape.invsig2 <- 2
prior$shape.hyper.invsig2 <- 0.2
prior$rate.hyper.invsig2 <- 0.1
prior$pi.split <- c(1, rep(0.5, prior$kmax - 2), 0)
prior$pi.birth <- c(1, rep(0.5, prior$kmax - 2), 0)
prior$Eb0.depend.mix <- FALSE
print(prior)

# Priors for regression parameters β

##### Gibbs step only
prior.beta.gibbs <- list()
prior.beta.gibbs$mean.prior <- rep(0, nregres)
prior.beta.gibbs$var.prior <- rep(5000, nregres)
print(prior.beta.gibbs)

##### random walk MH step only
prior.beta.mh1 <- list()
prior.beta.mh1$mean.prior <- rep(0, nregres)
prior.beta.mh1$var.prior <- rep(1000, nregres)
prior.beta.mh1$blocks <- list()
prior.beta.mh1$blocks$ind.block <- list()
prior.beta.mh1$blocks$ind.block[[1]] <- 1:nregres
nblock <- length(prior.beta.mh1$blocks$ind.block)

# vars: proposal variances for each beta parameter
# cors: lower triangle of the proposal correlation matrix
# corsm: proposal correlation matrix itself
# covm: proposal covariance matrix
vars <- c(0.1, 0.38, 0.4, 0.305, 0.1186, 0.1735, 0.7024, 0.346, 0.3920)
cors <- c(1, 0.1, 0, 0.1, 0.15, 0, 0.4, 0.1, 0.2, 1, -0.15, 0.15,
          -0.2, -0.3, 0.2, -0.1, 0, 1, -0.2, 0.15, 0.2, 0.3, 0.2, 0.1,
          1, 0.2, -0.5, 0.2, -0.4, 0.4, 1, 0.15, 0.5, 0.3, 0.4, 1,
          0.15, 0.15, 0, 1, 0.35, 0.65, 1, 0.2, 1)

corsm <- diag(nregres)
corsm[lower.tri(corsm, diag = TRUE)] <- cors
corsm[upper.tri(corsm, diag = FALSE)] <- t(corsm)[upper.tri(t(corsm), diag = FALSE)]
covm <- diag(sqrt(vars)) %*% corsm %*% diag(sqrt(vars))

prior.beta.mh1$blocks$cov.prop <- list()
prior.beta.mh1$blocks$cov.prop[[1]] <- covm[lower.tri(covm, diag = TRUE)]
prior.beta.mh1$type.upd <- rep("random.walk.metropolis", nblock)
prior.beta.mh1$weight.unif <- rep(0.5, nblock)
prior.beta.mh1$half.range.unif <- c(1, 1, 1, .1, 1,
                                    1, 1, .1, .1)
print(prior.beta.mh1)

##### two blocks: Gibbs + MH steps
prior.beta.mh2 <- list()
prior.beta.mh2$mean.prior <- rep(0, nregres)
prior.beta.mh2$var.prior <- rep(1000, nregres)
prior.beta.mh2$blocks <- list()
prior.beta.mh2$blocks$ind.block <- list()
prior.beta.mh2$blocks$ind.block[[1]] <- 1:6
prior.beta.mh2$blocks$ind.block[[2]] <- 7:9
nblock <- length(prior.beta.mh2$blocks$ind.block)

# vars: proposal variances for each beta parameter
# cors: lower triangle of the proposal correlation matrix
# corsm: proposal correlation matrix itself
# covm: proposal covariance matrix
vars <- c(0.1, 0.35, 0.35)
cors <- c(1, 0.9, 0.9, 1, 0.9, 1)
corsm <- diag(3)
corsm[lower.tri(corsm, diag = TRUE)] <- cors
corsm[upper.tri(corsm, diag = FALSE)] <- t(corsm)[upper.tri(t(corsm), diag = FALSE)]
covm <- diag(sqrt(vars)) %*% corsm %*% diag(sqrt(vars))
print(corsm); print(covm)

prior.beta.mh2$blocks$cov.prop <- list()
prior.beta.mh2$blocks$cov.prop[[1]] <- NULL
prior.beta.mh2$blocks$cov.prop[[2]] <- covm[lower.tri(covm, diag = TRUE)]
prior.beta.mh2$type.upd <- c("gibbs", "random.walk.metropolis")

prior.beta.mh2$weight.unif <- c(0.05, 0.05)
prior.beta.mh2$half.range.unif <- c(0.25, 0.25, 0.01, 1, 0.15, 0.25, 0.3, 1, 1)
print(prior.beta.mh2)

# Initials values of MCMC
nobs <- dim(X)[1]
init1 <- list()
init1$iter <- 0
init1$mixture <- c(1, 1, rep(0, prior$kmax - 1), 3.07, 
                   rep(0, prior$kmax - 1), 1.8, rep(0, prior$kmax - 1))
init1$beta <- c(1.1, -0.66, 0.04, -1.76, 0.94, 1.03, 0.37, 1.22,
                 0.82)
init1$y <- NULL
init1$r <- rep(1, nobs)
init1$otherp <- rgamma(1, shape = prior$shape.hyper.invsig2, rate = prior$rate.hyper.invsig2)
init1$u <- c(runif(1), 0, 0, runif(3 * (prior$kmax - 1)))
print(init1)

# Running the MCMC!
store <- list(y = FALSE, r = FALSE, u = FALSE, regresres = FALSE)
#nsimul <- list(niter = 5000, nthin = 3, nburn = 500, nnoadapt = 0, nwrite = 500)
nsimul <- list(niter = 20000, nthin = 3, nburn = 500, nnoadapt = 0, nwrite = 500)
 
#dir.create("../../Codes/DPMIV_analysis/bayesianAFTtest1")
dirsim1test <- paste("../../Codes/DPMIV_analysis/bayesSurv/bayesianAFTtest1", sep = "")
simul1 <- bayessurvreg1(Surv(time=exp(L), time2=exp(R), type="interval2") ~ First+ SEX + AGE_RECRUITMENT +
                          CHOLESTEROL + HDL + SMOKE, 
                          data = dm_filtered, dir = dirsim1test, nsimul = nsimul, prior = prior,
                          prior.beta = prior.beta.gibbs,
                          init = init1, store = store)

ca = read.table("../../Codes/DPMIV_analysis/bayesSurv/bayesianAFTtest1/mixmoment.sim", header = T)
