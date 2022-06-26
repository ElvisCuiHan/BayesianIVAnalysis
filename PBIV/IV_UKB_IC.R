library("reshape2") 
library(tidyverse)
library(bayesplot)

rm(list=ls()) 

setwd("~/Downloads/Dr Li GSR/UKB Data/Interval_Censoring_Project_new/")

dm_filtered_ = read.delim(file = "0325_DPMIV_test_data_eth1_female_subsetsnp_2clust.dat", sep=" ", header=T)
dm_filtered_ = read.delim(file = "0330_DPMIV_test_data_eth1_male_subsetsnp_2clust.dat", sep=" ", header=T)

dm_filtered = dm_filtered_[c(1:6666), ]

setwd("~/Downloads/Dr Li GSR/Codes/PBIV")
#source("IV_MH_one_step.R")
source("IV_MH_Interval_Censored.R")

# Outcome
L = dm_filtered$L
R = dm_filtered$R
d = dm_filtered$d

# Noisy surrogate
X = dm_filtered$LOG_SBP

# Instruments
G = dm_filtered[, c(5:19)]
U = dm_filtered[, c(20:24)]

# wid = a vector of the random walk width for (a0,a1,a2,sig1,b0,b1,b2,sig2,rho)
Width = c(0.025, # a0
          rep(0.012, 15), # a1
          3 * c(0.008, 0.0048, 0.002, 0.004, 0.0026), # a2
          0.005, # sig1 (22)
          0.08, # b0
          .018, # b1 (24)
          0.08, 0.08, 0.12, 0.03, 0.052, # b2
          0.28, # sig2
          0.058) # rho
#INIT = c(rep(0, 21), .01, rep(0, 7), .1, 0.)
hi = read_delim("0625_female_9865_chain1.txt")
INIT = as.matrix(hi[dim(hi)[1], ])
ITER = 120000

# IV model
fit = IV_MH_IC(
  L=L,
  R=R,
  d=d, 
  X=X, 
  G=G, 
  U=U, 
  m=ITER, 
  # adjust width to give reasonable acceptance rate
  wid=5*Width, 
  init = INIT,
  # priors: N(0,10000) for a0,a1,a2,b0,b1,b2; Inv-Gamma(0.001, 0.001) for sig1 and sig2
  prior_1=c(rep(0, 21), .1, rep(0., 7), .1),
  prior_2=c(rep(0.8, 21), .1, rep(0.8, 7), .1)
) 

hi = as.data.frame(fit[1:9])
#write_delim(hi, "0625_female_9865_chain1.txt")

# acceptance rate for (a0,a1,a2,sig1,b0,b1,b2,sig2,rho)
round(fit$ap, 3)
# trace plot
acf(fit$b1r, lag.max = 5000)
plot(fit$b1r, type='l')
hi %>% mutate(index = 1:ITER) %>%
  pivot_longer(cols=b2r.1:b2r.2, values_to = c("value"), names_to = "para") %>%
  ggplot() + geom_line(aes(x=index, y=value, color=para), alpha=0.4) +
  #ylim(0, 0.02) +
  xlab("") + theme_bw()
# posterior samples of b1 with 1000 burn-in
b1r_res=fit$b1r[-c(1:0000)]
# histogram
hist(b1r_res)
# posterior means
mean(b1r_res)
sd(b1r_res)
# 95% credible interval
quantile(b1r_res,c(0.025,0.975))

