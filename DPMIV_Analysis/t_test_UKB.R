rm(list = ls())
library("reshape2") 
library(tidyverse)
library(jcolors)

output_beta = function(thinning=200, burn_in=2e5) {
  files = list.files("./")
  files = files[str_detect(files, "temp.txt")]
  #files = files[str_detect(files, "trace_a2")]
  fnum <- length(files)
  
  c1 = read_table(files[1])$b1
  ll <- length(c1)
  result <- tibble(c1 = read_table(files[1])$b1[1:ll])
  for (i in 2:fnum){
    ca = read_table(files[i])$b1[1:ll]
    name = paste("c", as.character(i), sep="")
    print(name)
    result[name] = ca
  }
  
  out = result[burn_in:ll, 1:51] %>% slice(which(row_number() %% thinning == 1)) %>%
    mutate(., "ind"=c(1:trunc((ll-burn_in)/thinning+1.0))) %>% pivot_longer(c1:c51) %>%
    select(value)
  rm(rsult)
  out
}

setwd("~/Downloads/Dr Li GSR/Results/DPMIV_MCMC_results/0406-run-51chain-result-female/")

female = output_beta(thinning = 1000)$vale

setwd("~/Downloads/Dr Li GSR/Results/DPMIV_MCMC_results/0414-run-51chain-result-male/male-temp/")

male = output_beta(thinning = 1000)$value

####### Histogram #######

hist(female, breaks=75, freq = F, xlab="", main="", col="white")

hist(male, breaks=75, freq = F, xlab="", main="", col="white")

####### Test #######

qqnorm(female)

qqnorm(male)

t.test(female, male)
