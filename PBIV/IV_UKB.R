library("reshape2") 
library(tidyverse)
library(bayesplot)

rm(list=ls()) 

setwd("~/Downloads/Dr Li GSR/UKB Data/Interval_Censoring_Project_new/")

dm_filtered_ = read.delim(file = "0325_DPMIV_test_data_eth1_female_colnames.dat", sep=" ", header=T) %>% 
  filter(., d!=2)
dm_filtered_ = read.delim(file = "0330_DPMIV_test_data_eth1_male_subsetsnp_2clust.dat", sep=" ", header=T) %>% 
  filter(., d!=2)
#ca = read_delim(file="DPMIV_test_data.dat")
#colnames(dm_filtered_) = colnames(ca)[-23]
dm_filtered_$d = dm_filtered_$d - 3

dm_filtered = dm_filtered_[c(1:1280), ]

setwd("~/Downloads/Dr Li GSR/Codes/parametric IV model")
source("IV_MH_one_step.R")

# Outcome
Y = (dm_filtered$R)
d = dm_filtered$d

# Noisy surrogate
X = dm_filtered$LOG_SBP

# Instruments
G = dm_filtered[, c(5:19)]
U = dm_filtered[, c(20:24)]

# wid = a vector of the random walk width for (a0,a1,a2,sig1,b0,b1,b2,sig2,rho)
Width = c(0.015, # a0
          rep(0.008, 15), # a1
          c(0.0048, 0.0048, 0.002, 0.004, 0.0026), # a2
          0.001, # sig1 (22)
          0.05, # b0
          .015, # b1 (24)
          0.04, 0.05, 0.08, 0.02, 0.042, # b2
          0.26, # sig2
          0.048) # rho
#INIT = c(rep(0, 21), .01, rep(0, 7), .1, 0.)
hi = read_delim("0330_male_burnin.txt")
INIT = as.matrix(hi[dim(hi)[1], ])
ITER = 120000

# IV model
fit = IV_MH(
  Y=Y, 
  d=d, 
  X=X, 
  G=G, 
  U=U, 
  m=ITER, 
  # adjust width to give reasonable acceptance rate
  wid=Width, 
  init = INIT,
  # priors: N(0,10000) for a0,a1,a2,b0,b1,b2; Inv-Gamma(0.001, 0.001) for sig1 and sig2
  prior_1=c(rep(0, 21), .1, rep(0, 7), .1),
  prior_2=c(rep(1, 21), .1, rep(1, 7), .1)
) 

hi = as.data.frame(fit[1:9])
#write_delim(hi, "0330_male_burnin.txt")

# acceptance rate for (a0,a1,a2,sig1,b0,b1,b2,sig2,rho)
round(fit$ap, 3)
# trace plot
acf(fit$b1r, lag.max = 5000)
plot(fit$b1r, type='l')
hi %>% mutate(index = 1:ITER) %>%
  pivot_longer(cols=a1r.1:a1r.15, values_to = c("value"), names_to = "para") %>%
  ggplot() + geom_line(aes(x=index, y=value, color=para), alpha=0.4) +
  #ylim(0, 0.02) +
  xlab("") + theme_bw()
# posterior samples of b1 with 1000 burn-in
b1r_res=fit$b1r[-c(1:000)]
# histogram
hist(b1r_res)
# posterior means
mean(b1r_res)
sd(b1r_res)
# 95% credible interval
quantile(b1r_res,c(0.025,0.975))

