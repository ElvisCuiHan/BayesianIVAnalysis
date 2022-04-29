library("reshape2") 
library(tidyverse)
library(bayesplot)

rm(list=ls()) 
setwd("~/Downloads/Dr Li GSR/Codes/DPMIV_analysis")

#ca = read.table("DPMIV_MCMC_results/run-40chain-result/0315-162020_ch=220K_temp.txt", header = T)
#ca = read.table("../DPMIV_UKB/DPMIV_UKB-0318/0322-161631_ch=1100K_temp.txt", header=T)[-c(1:10000),]
#ca = read.table("../../../Results/DPMIV_MCMC_results/0330-run-51chain-result/0326-153659_ch=3200K_b1=-3.000_temp.txt", header=T)
ca = read.table("../DPMIV_UKB-0325-SNP/0331-005226_ch=1100K_b1=1.744_temp.txt", header=T)[,]

hist((ca$b1), breaks=100, main = expression(paste("Histograms of the posterior samples of ", beta[1])),
     xlab = expression(paste(beta[1])))

acf(ca$b1, lag.max = 5000, main="Autocorrelation plot")

thinning = 200
mean(ca$b1[c(rep(FALSE, thinning), TRUE)])
sd(ca$b1[c(rep(FALSE, thinning), TRUE)])

CI = quantile(ca$b1[c(rep(FALSE, thinning), TRUE)], probs=c(0.025,0.5,0.975))
ggplot() + geom_histogram(aes(ca$b1[c(rep(FALSE, thinning), TRUE)]), bins = 45,
                          color="lightblue", fill="orange") +
  xlab("") + ylab(expression(paste(beta[1]))) +
  ggtitle(expression(paste("Trace plots of the posterior samples of ", beta[1], 
                           ", thinning=", 20))) +
  geom_errorbarh(aes(xmin=CI[1], xmax=CI[3], y=-0.5), col="#0094EA", size=2.2) +
  theme_bw() +
  geom_vline(xintercept = CI[2], cex=1.2, color="purple")

plot(ca$b1[c(rep(FALSE, thinning), TRUE)], type="l", ylab = expression(paste(beta[1])),
     xlab = "", main=expression(paste("Trace plots of the posterior samples of ", beta[1], 
                                      ", thinning=", 20)))

hist(ca$ncl, breaks=100, main = expression(paste("Histograms of the number of clusters")),
     xlab = "ncluster")
plot(ca$ncl[c(rep(FALSE, thinning), TRUE)], type="l")

############ Mixing Diagonasis ##################

rm(list=ls())

library(viridis)

setwd("~/Downloads/Dr Li GSR/Codes/DPMIV_analysis/DPMIV_MCMC_results/0316-run-40chain-result/")
setwd("~/Downloads/Dr Li GSR/Codes/DPMIV_UKB/DPMIV_UKB-0325-SNP/")
files = list.files("./")
files = files[str_detect(files, "temp.txt")]

result = tibble(c1 = read_table(files[1], col_names = T)$b1)
ll = length(result$c1)
for (i in 2:length(files)){
  ca = read_table(files[i])$b1[1:ll]
  name = paste("c", as.character(i), sep="")
  print(name)
  result[name] = ca
}

thinning = 2
burn_in = 000

to_plot = result[burn_in:ll, 1:5] %>% slice(which(row_number() %% thinning == 1)) %>%
  mutate(., "ind"=c(1:trunc((ll-burn_in)/thinning+1.0))) %>% pivot_longer(c1:c5)

ggplot(aes(x=ind, y=value), data=to_plot) + geom_line(aes(color=name, alpha=0.01)) +
  #scale_color_brewer(palette = "Spectral") +
  theme_bw() + theme(legend.position="none") +
  ggtitle(expression(paste("Trace plots of the posterior samples of ", beta[1], 
                           ", thinning=", 20))) + #ylim(-2.5, 1.5) +
  xlab("Diagonasis of Chain Mixing")

############ INDIVIDUAL LEVEL RESULT ############

ind <- read.delim("ind_result.txt")
View(ind[329000:329999,])

get.ESS <- function(x, autocor = 1, ESS = 1){
  M <- length(x)
  # Re-centre at zero
  x <- x - mean(x)
  # Autocorrelations
  rho <- rep(NA, M-2)
  for(lag in 0:(M-2)){
    # Split time-series in two, separated by 'lag'
    theta.1 <- x[1:(M-lag)]
    theta.2 <- x[(1+lag):M]
    # Autocorrelation
    if(autocor %in% c(1, 4)){
      rho[lag + 1] <- cor(theta.2, theta.1)
    }else if(autocor == 2){
      rho[lag + 1] <- mean(theta.1 * theta.2) / var(x)
    }else if(autocor == 5){
      rho[lag + 1] <- cov(theta.2, theta.1)
    }
  }
  if(autocor == 3){
    rho <- stats::acf(x = x, lag.max = M-2, plot = FALSE)$acf
  }else if(autocor == 4){
    # 'Relative' autocorrelation
    rho <- rho * (M - 0:(M-2)) / M
  }
  
  # Effective sample size
  if(ESS == 1){
    E <- M / (1 + 2 * sum(rho))
  }else if(ESS == 2){
    E <- M / (1 + 2 * sum(abs(rho))) 
  }else if(ESS == 3){     # To be used with autocor = 5
    lambda.sq <- M * var(x)
    sigma.sq <- lambda.sq + 2 * sum(rho)
    E <- M * lambda.sq / sigma.sq
  }
  return(E)
}
get.ESS(ca$b1[3000:4000], 3, 3)
mcmcse::ess(ca$b1)
coda::effectiveSize(ca$b1)
