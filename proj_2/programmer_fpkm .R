library(ggplot2)
library(ggpubr)
library(gplots)
library(tidyverse)

set.seed(5) 
# read data
data=read.table("/projectnb/bf528/users/lava_lamp/project_5/divya_pr5/P0_1_cufflinks/genes.fpkm_tracking", header=TRUE)


# removing zero fpkm values
genefpkm_trimmed <- data[data$FPKM!=0, ]



# checking FPKM field data, remove duplicates
fpkm  <-as.data.frame(genefpkm_trimmed)
fpkm  <-fpkm %>% select(FPKM,gene_short_name,gene_id)
fpkm2 <-fpkm

fpkm2 <-fpkm2 %>% distinct(fpkm2$gene_short_name, .keep_all = TRUE)

# Total genes after removing duplicates
nrow(fpkm2)


# Histogram 
# png("genefpkm_density.png")
fpkmplotdata<-as.numeric(genefpkm_trimmed$FPKM, center=TRUE)
max(fpkm)
fpkm_sd = sd(fpkm)
fpkm_mean = mean(fpkm)
fpkm_median = median(fpkm)
hist(fpkm,probability=T,xlab="FPKM values",main=paste(c("Mean=",fpkm_mean,";","SD=",fpkm_sd),collapse=""),border="blue")
lines(density(fpkm,bw=10),col='red')
# dev.off()
