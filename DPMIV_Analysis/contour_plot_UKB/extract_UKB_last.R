rm(list=ls())

setwd("~/Downloads/Dr Li GSR/Results/DPMIV_MCMC_results/0406-run-51chain-result-female/")
setwd("~/Downloads/Dr Li GSR/Results/DPMIV_MCMC_results/0414-run-51chain-result-male/male-temp/")
#setwd("~/Downloads/Dr Li GSR/Results/DPMIV_MCMC_results/0328-run-51chain-result/SNP-results(n=760K)-trace-a2b2/")

files = list.files("./")
files = files[str_detect(files, "last")]
#files = files[str_detect(files, "trace_a2")]
fnum <- length(files)

ca = read.table(files[1], sep=",", header = T)
ca = matrix(as.numeric(unlist(strsplit(ca[6:(length(ca[,1])-1), ], " "))),
            ncol=6, byrow=TRUE)
#colnames(ca) =  c("mu1",	"sig1",	"mu2",	"sig2",	"rho",	"n")

result <- tibble(mu1 = ca[,1], sig1 = ca[,2], mu2 = ca[,3], sig2 = ca[,4], rho=ca[,5], n=ca[,6])
for (i in 2:fnum){
  ca = read.table(files[i], sep=",", header = T)
  ca = matrix(as.numeric(unlist(strsplit(ca[6:(length(ca[,1])-1), ], " "))),
              ncol=6, byrow=TRUE)
  colnames(ca) =  c("mu1",	"sig1",	"mu2",	"sig2",	"rho",	"n")
  ca = as.data.frame(ca)
  print(i)
  result = result %>% add_row(ca)
}

write_delim(result, "../../../../Codes/DPMIV_Analysis/contour_plot_UKB/male_last_result.txt")
