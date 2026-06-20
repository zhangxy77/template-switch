setwd("./20_mix/")
library(tidyverse)
library(dplyr)
library(data.table)
library(phylotools)
load("./F123_16S2_20_mix.RData")
ts_data <- umi_ts[umi_ts$ts==1,]
save(ts_data,file="./ts_Data.RData")
for(i in seq(1:nrow(ts_data)))
{
  
  output_file <- paste("./muscle_data/",i,".fasta",sep = "")
  seq_data <- matrix(0,nrow = 6,ncol = 1)
  seq_data[1,1] <- paste(">",ts_data$umi3[i],sep = "")
  seq_data[2,1] <- ts_data$mg_minor[i]
  seq_data[3,1] <- paste(">",ts_data$ms_major[[i]],"_A",sep = "")
  seq_data[4,1] <- ts_data$ref_major[[i]]
  seq_data[5,1] <- paste(">",ts_data$ms_minor[[i]],"_B",sep = "")
  seq_data[6,1] <- ts_data$ref_minor[[i]]
  write.table(seq_data,output_file,quote = F,col.names = F,row.names = F)

}
