rm(list = ls())
library("reshape2") 
library(tidyverse)
library(jcolors)

setwd("~/Downloads/Dr Li GSR/Results/DPMIV_MCMC_results/0406-run-51chain-result-female/")
setwd("~/Downloads/Dr Li GSR/Results/DPMIV_MCMC_results/0414-run-51chain-result-male/male-temp/")
#setwd("~/Downloads/Dr Li GSR/Results/DPMIV_MCMC_results/0328-run-51chain-result/SNP-results(n=760K)-trace-a2b2/")

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

c1 = read_table(files[1], col_names = F)$X5
ll <- length(c1)
result <- tibble(c1 = read_table(files[1], col_names = F)$X5[1:ll])
for (i in 2:fnum){
  ca = read_table(files[i], col_names = F)$X5[1:ll]
  name = paste("c", as.character(i), sep="")
  print(name)
  result[name] = ca
}

thinning = 200
burn_in = 5e4

to_plot = result[burn_in:ll, 1:51] %>% slice(which(row_number() %% thinning == 1)) %>%
  mutate(., "ind"=c(1:trunc((ll-burn_in)/thinning+1.0))) %>% pivot_longer(c1:c51)

ggplot(aes(x=ind, y=value), data=to_plot) + 
  geom_line(aes(color=(name)), alpha=0.5, cex=0.7) +
  #scale_colour_stata("s2color") +
  #scale_color_brewer(palette = "stata_colors") +
  theme_bw() + 
  theme(legend.position="none") +
  ggtitle(expression(paste("Trace plots of the posterior samples of ", beta[1], 
                           ", thinning=", 20))) + #ylim(-2, 2.) +
  xlab("Diagonasis of Chain Mixing")

plot_mixing <- function(thinning = 100, low=-2, upp=2, lwd=1) {
  plot(result_burn_in$c1[c(rep(FALSE, thinning-1), TRUE)], type="l", col="red", lwd=lwd, ylim = c(low, upp),
       main=expression(paste("Trace plots of the posterior samples of ", beta[1])),
       xlab = (paste("Diagonasis of Chain Mixing, Thinning = ", thinning)),
       ylab = "")
  lines(result_burn_in$c2[c(rep(FALSE, thinning-1), TRUE)], type="l", col="orange", lwd=lwd)
  lines(result_burn_in$c3[c(rep(FALSE, thinning-1), TRUE)], type="l", col="yellow", lwd=lwd)
  lines(result_burn_in$c4[c(rep(FALSE, thinning-1), TRUE)], type="l", col="green", lwd=lwd)
  lines(result_burn_in$c5[c(rep(FALSE, thinning-1), TRUE)], type="l", col="dodgerblue", lwd=lwd)
  lines(result_burn_in$c6[c(rep(FALSE, thinning-1), TRUE)], type="l", col="blue", lwd=lwd)
  lines(result_burn_in$c7[c(rep(FALSE, thinning-1), TRUE)], type="l", col="purple", lwd=lwd)
  lines(result_burn_in$c8[c(rep(FALSE, thinning-1), TRUE)], type="l", col="grey", lwd=lwd)
  lines(result_burn_in$c9[c(rep(FALSE, thinning-1), TRUE)], type="l", col="black", lwd=lwd)
  lines(result_burn_in$c10[c(rep(FALSE, thinning-1), TRUE)], type="l", col="red", lwd=lwd)
  lines(result_burn_in$c11[c(rep(FALSE, thinning-1), TRUE)], type="l", col="orange", lwd=lwd)
  lines(result_burn_in$c12[c(rep(FALSE, thinning-1), TRUE)], type="l", col="gold", lwd=lwd)
  lines(result_burn_in$c13[c(rep(FALSE, thinning-1), TRUE)], type="l", col="lightgreen", lwd=lwd)
  lines(result_burn_in$c14[c(rep(FALSE, thinning-1), TRUE)], type="l", col="dodgerblue", lwd=lwd)
  lines(result_burn_in$c15[c(rep(FALSE, thinning-1), TRUE)], type="l", col="black", lwd=lwd)
  lines(result_burn_in$c16[c(rep(FALSE, thinning-1), TRUE)], type="l", col="orange", lwd=lwd)
  lines(result_burn_in$c17[c(rep(FALSE, thinning-1), TRUE)], type="l", col="yellow", lwd=lwd)
  lines(result_burn_in$c18[c(rep(FALSE, thinning-1), TRUE)], type="l", col="green", lwd=lwd)
  lines(result_burn_in$c19[c(rep(FALSE, thinning-1), TRUE)], type="l", col="dodgerblue", lwd=lwd)
  lines(result_burn_in$c20[c(rep(FALSE, thinning-1), TRUE)], type="l", col="blue", lwd=lwd)
  lines(result_burn_in$c21[c(rep(FALSE, thinning-1), TRUE)], type="l", col="purple", lwd=lwd)
  lines(result_burn_in$c22[c(rep(FALSE, thinning-1), TRUE)], type="l", col="grey", lwd=lwd)
  lines(result_burn_in$c23[c(rep(FALSE, thinning-1), TRUE)], type="l", col="black", lwd=lwd)
  lines(result_burn_in$c24[c(rep(FALSE, thinning-1), TRUE)], type="l", col="red", lwd=lwd)
  lines(result_burn_in$c25[c(rep(FALSE, thinning-1), TRUE)], type="l", col="orange", lwd=lwd)
  lines(result_burn_in$c26[c(rep(FALSE, thinning-1), TRUE)], type="l", col="gold", lwd=lwd)
  lines(result_burn_in$c27[c(rep(FALSE, thinning-1), TRUE)], type="l", col="lightgreen", lwd=lwd)
  lines(result_burn_in$c28[c(rep(FALSE, thinning-1), TRUE)], type="l", col="dodgerblue", lwd=lwd)
  lines(result_burn_in$c29[c(rep(FALSE, thinning-1), TRUE)], type="l", col="purple", lwd=lwd)
  lines(result_burn_in$c30[c(rep(FALSE, thinning-1), TRUE)], type="l", col="orange", lwd=lwd)
  lines(result_burn_in$c31[c(rep(FALSE, thinning-1), TRUE)], type="l", col="yellow", lwd=lwd)
  lines(result_burn_in$c32[c(rep(FALSE, thinning-1), TRUE)], type="l", col="green", lwd=lwd)
  lines(result_burn_in$c33[c(rep(FALSE, thinning-1), TRUE)], type="l", col="dodgerblue", lwd=lwd)
  lines(result_burn_in$c34[c(rep(FALSE, thinning-1), TRUE)], type="l", col="blue", lwd=lwd)
  lines(result_burn_in$c35[c(rep(FALSE, thinning-1), TRUE)], type="l", col="purple", lwd=lwd)
  lines(result_burn_in$c36[c(rep(FALSE, thinning-1), TRUE)], type="l", col="grey", lwd=lwd)
  lines(result_burn_in$c37[c(rep(FALSE, thinning-1), TRUE)], type="l", col="black", lwd=lwd)
  lines(result_burn_in$c38[c(rep(FALSE, thinning-1), TRUE)], type="l", col="red", lwd=lwd)
  lines(result_burn_in$c39[c(rep(FALSE, thinning-1), TRUE)], type="l", col="orange", lwd=lwd)
  lines(result_burn_in$c40[c(rep(FALSE, thinning-1), TRUE)], type="l", col="gold", lwd=lwd)
  lines(result_burn_in$c41[c(rep(FALSE, thinning-1), TRUE)], type="l", col="lightgreen", lwd=lwd)
  lines(result_burn_in$c42[c(rep(FALSE, thinning-1), TRUE)], type="l", col="brown", lwd=lwd)
  lines(result_burn_in$c43[c(rep(FALSE, thinning-1), TRUE)], type="l", col="red", lwd=lwd)
  lines(result_burn_in$c44[c(rep(FALSE, thinning-1), TRUE)], type="l", col="dodgerblue", lwd=lwd)
  lines(result_burn_in$c45[c(rep(FALSE, thinning-1), TRUE)], type="l", col="purple", lwd=lwd)
  lines(result_burn_in$c46[c(rep(FALSE, thinning-1), TRUE)], type="l", col="gold", lwd=lwd)
  lines(result_burn_in$c47[c(rep(FALSE, thinning-1), TRUE)], type="l", col="lightgreen", lwd=lwd)
  lines(result_burn_in$c48[c(rep(FALSE, thinning-1), TRUE)], type="l", col="blue", lwd=lwd)
  lines(result_burn_in$c49[c(rep(FALSE, thinning-1), TRUE)], type="l", col="purple", lwd=lwd)
  lines(result_burn_in$c50[c(rep(FALSE, thinning-1), TRUE)], type="l", col="black", lwd=lwd)
  lines(result_burn_in$c51[c(rep(FALSE, thinning-1), TRUE)], type="l", col="yellow", lwd=lwd)
  #lines(result_burn_in$c52[c(rep(FALSE, thinning-1), TRUE)], type="l", col="darkgreen", lwd=lwd)
  #lines(result_burn_in$c53[c(rep(FALSE, thinning-1), TRUE)], type="l", col="gold", lwd=lwd)
}

burn_in = 0e4
result_burn_in = result[burn_in:ll, 1:51]
plot_mixing(200, -1.2, .65, 2)
#plot_mixing(1000, min(to_plot$value, na.rm=T),
#                  max(to_plot$value, na.rm=T), 2.)

png(file="../../Figures_for_paper/trace_female.png",
    width=1500, height=600)
plot_mixing(200, -1.2, .65, 2)
dev.off()

png(file="../../../Figures_for_paper/trace_male.png",
    width=1500, height=600)
plot_mixing(200, -1.5, 1., 2)
dev.off()

########## ANALYSIS OF CHAIN CUTTING ################

thinning = 200
kaishi = 2e6
jieshu = ll

to_cal = result[kaishi:jieshu, 1:51] %>% slice(which(row_number() %% thinning == 1)) %>%
  mutate(., "ind"=c(1:trunc((jieshu - kaishi)/thinning+1.0))) %>% pivot_longer(c1:c51)

mean(to_cal$value, na.rm = T); sd(to_cal$value, na.rm = T)
quantile(to_cal$value, probs = c(0.025, 0.5, 0.975), na.rm=T)
