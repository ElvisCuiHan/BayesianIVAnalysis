library("reshape2") 
library(tidyverse)
library(bayesplot)

rm(list=ls()) 
#setwd("~/Downloads/Dr Li GSR/Codes/DPMIV_analysis")
setwd("/Users/david/Desktop/coding/DPMIV_UKB-0323-40times/")

filename <- "0324-092500_ch=1100K_temp.txt"
#ca = read_delim(filename, col_names = F)
#ca = read.delim("acceptrate_a1.txt", sep=" ", header=F)
ca = read.table(filename, header = T)

#ca = ca[-c(1:1000),]

hist((ca$b1), breaks=100, main = expression(paste("Histograms of the posterior samples of ", beta[1])),
     xlab = expression(paste(beta[1])))

acf(ca$b1, lag.max = 10000, main="Autocorrelation plot")

thinning = 20
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

plot(ca$b1, type="l")


setwd("/Users/david/Desktop/coding/DPMIV_UKB-0323-40times/")

filename <- "0324-233029_ch=1100K_b1=-0.500_temp.txt"
ca = read_table(filename, col_names=TRUE)
plot(ca$b1, type="l", ylab ="b1 trace")
plot(ca$b2, type="l", ylab ="b2 trace")
plot(ca$a1, type="l", ylab ="a1 trace")
plot(ca$a2, type="l", ylab ="a2 trace")

filename <- "0324-092500_ch=1100K_trace_a1.txt"
ca = read_table(filename, col_names=FALSE)
plot(ca$X1, type="l", ylab ="a1 trace")

filename <- "0324-092500_ch=1100K_trace_a2.txt"
ca = read_table(filename, col_names=FALSE)
plot(ca$X1, type="l", ylab ="a2 trace")

filename <- "0324-092500_ch=1100K_trace_b2.txt"
ca = read_table(filename, col_names=FALSE)
plot(ca$X1, type="l", ylab ="b2 trace")


############ Mixing Diagonasis ##################

thinning = 20
burnin = 10000

ca6 = read_delim("temp_n_30000_trial_6.txt")[-c(1:burnin), ]
ca5 = read_delim("temp_n_30000_trial_5.txt")[-c(1:burnin), ]
ca4 = read_delim("temp_n_30000_trial_4.txt")[-c(1:burnin), ]
ca3 = read_delim("temp_n_30000_trial_3.txt")[-c(1:burnin), ]
ca2 = read_delim("temp_n_30000_trial_2.txt")[-c(1:burnin), ]
ca1 = read_delim("temp_n_30000_trial_1.txt")[-c(1:burnin), ]

t1 = c(1:length(ca3$b1))[c(rep(FALSE, thinning), TRUE)]
hi = tibble(t = 1:length(t1),
            c1 = ca1$b1[c(rep(FALSE, thinning), TRUE)],
            c2 = ca2$b1[c(rep(FALSE, thinning), TRUE)],
            c3 = ca3$b1[c(rep(FALSE, thinning), TRUE)],
            c4 = ca4$b1[c(rep(FALSE, thinning), TRUE)],
            c5 = ca5$b1[c(rep(FALSE, thinning), TRUE)],
            c6 = ca6$b1[c(rep(FALSE, thinning), TRUE)]) %>%
  pivot_longer(cols = c1:c6)

ggplot(aes(x=t,y=value), data = hi) +
  geom_line(aes(color=name), cex=0.3) +
  xlab("Iteration") + ylab(expression(paste(beta[1]))) +
  ylim(-30, 6) +
  theme_bw() +
  theme(axis.title.y=element_text(size=14)) +
  ggtitle("Diagonasis of Chain Mixing")

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
