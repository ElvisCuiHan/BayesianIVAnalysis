library("reshape2") 
library(tidyverse)
library(bayesplot)

rm(list=ls()) 

#setwd("~/Downloads/Dr Li GSR/UKB Data/Interval_Censoring_Project_new/")
setwd("/Users/david/Desktop/coding/DPMIV_UKB-0322")

dm_filtered_ = read.delim(file = "0321_DPMIV_UKB_5U.dat", sep=" ", header=T) %>% 
  filter(., d!=2)
#ca = read_delim(file="DPMIV_test_data.dat")
#colnames(dm_filtered_) = colnames(ca)[-23]
dm_filtered_$d = dm_filtered_$d - 3

dm_filtered = dm_filtered_[, ]

#setwd("~/Downloads/Dr Li GSR/Codes/parametric IV model")
source("IV_MH.R")

# Outcome
Y = (dm_filtered$R)
d = dm_filtered$d

# Noisy surrogate
X = dm_filtered$First

# Instruments
G = dm_filtered$Baseline
U = dm_filtered[, c(6:10)]

# wid = a vector of the random walk width for (a0,a1,a2,sig1,b0,b1,b2,sig2,rho)
Width = c(0.01, # a0
          c(0.003, 0.01, 0.0002, 0.002, 0.007, 0.02), # a1 a2
          0.001, # sig1
          0.5, # b0
          .1, # b1
          0.4, rep(0.02, 4), # b2
          1, # sig2
          0.1) # rho
Width[15] = 0.5; Width[12] = 0.005; Width[14] = .5
#INIT = c(rep(0, 10), .01, rep(-0, 10), .1, 0.1)
#INIT = as.matrix(hi[11000, ])
ITER = 110000

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
  #init = INIT,
  # priors: N(0,10000) for a0,a1,a2,b0,b1,b2; Inv-Gamma(0.001, 0.001) for sig1 and sig2
  prior_1=c(rep(0, 7), .001, rep(0, 7), .001),
  prior_2=c(rep(10000,7), .001, rep(10000,7), .001)
) 

hi = as.data.frame(fit[1:9])
write_delim(hi, "0322temp.txt")

# acceptance rate for (a0,a1,a2,sig1,b0,b1,b2,sig2,rho)
round(fit$ap, 3)
# trace plot
acf(fit$b1r, lag.max = 5000)
plot(fit$b1r, type='l')
hi %>% mutate(index = 1:ITER) %>%
  pivot_longer(cols=b2r.1:b2r.3, values_to = c("value"), names_to = "para") %>%
  ggplot() + geom_line(aes(x=index, y=value, color=para), alpha=0.4) +
  #ylim(0, 0.02) +
  xlab("")
# posterior samples of b1 with 1000 burn-in
b1r_res=fit$b1r[-c(1:1000)]
# histogram
hist(b1r_res)
# posterior means
mean(b1r_res)
sd(b1r_res)
# 95% credible interval
quantile(b1r_res,c(0.025,0.975))

