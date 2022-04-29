library("reshape2") 
library(tidyverse)
library(bayesplot)
#library(ggplot)

setwd("/Users/david/Desktop/coding/DPMIV_UKB-0325-SNP/")
filename1 <- "0326-131424_ch=110K_b1=-0.300_temp.txt"
filename2 <- "0326-131428_ch=110K_b1=-0.200_temp.txt"
filename3 <- "0326-131432_ch=110K_b1=-0.100_temp.txt"
hi1 = read_table(filename1)$b1
hi2 = read_table(filename2)$b1
hi3 = read_table(filename3)$b1
ggplot() +
  geom_line(aes(x=1:length(hi1),y=hi1),color="red") +
  geom_line(aes(x=1:length(hi2),y=hi2),color="green") +
  geom_line(aes(x=1:length(hi3),y=hi3),color="blue") +
  theme_bw()


filename1 <- "0326-144802_ch=3200K_b1=-1.356_temp.txt"
filename2 <- "0326-144735_ch=3200K_b1=-0.090_temp.txt"
filename3 <- "0326-144733_ch=3200K_b1=2.226_temp.txt"
filename4 <- "0326-145953_ch=3200K_b1=-2.694_temp.txt"
hi1 = read_table(filename1)$b1
hi2 = read_table(filename2)$b1
hi3 = read_table(filename3)$b1
hi4 = read_table(filename4)$b1
ggplot() +
  geom_line(aes(x=1:length(hi1),y=hi1),color="red") +
  geom_line(aes(x=1:length(hi2),y=hi2),color="green") +
  geom_line(aes(x=1:length(hi3),y=hi3),color="blue") +
  geom_line(aes(x=1:length(hi4),y=hi4),color="black") +
  theme_bw()



filename1 <- ""
filename2 <- ""
filename3 <- ""

ca = read_table(filename, col_names=TRUE)
plot(ca$b1, type="l", ylab ="b1 trace")
plot(ca$b2, type="l", ylab ="b2 trace")
plot(ca$a1, type="l", ylab ="a1 trace")
plot(ca$a2, type="l", ylab ="a2 trace")