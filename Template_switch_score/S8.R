setwd("./16S2_2reads_minor_identity/")
library(tidyverse)
library(dplyr)
library(data.table)
library(ggplot2)
library(stringr)
library(Biostrings)
library(phylotools)
#20
load("./umi_ts_20.RData")
ts_data <- umi_ts_20
ts_data <- ts_data[ts_data$ts==1,]
ts_data$index <- seq(1,length(ts_data[,1]))
ts_data_fil <- ts_data[ts_data$ms_major!=ts_data$ms_minor,]
load("./20_mix/ts_location_10_2.RData")
ts_location$index <- ts_data_fil$index
ts_location <- ts_location[order(ts_location$index,decreasing = F),]
sum(ts_location$index == ts_data_fil$index)
ts_location_part <- ts_location[ts_location$ts_length>1,]
ts_data_fil <- ts_data_fil[ts_data_fil$index %in% ts_location_part$index == T,]
S8_20_1 <- data.frame(cycle=20,umi3=ts_data_fil$umi3,umi5=ts_data_fil$umi5_minor,
                    ref_3=ts_data_fil$ms_major,ref_5=ts_data_fil$ms_minor,Chimeras=ts_data_fil$mg_minor,ts_region=ts_location_part$change_seq)
ts_data_3 <- ts_data[ts_data$ms_major==ts_data$ms_minor,]
S8_20_2 <- data.frame(cycle=20,umi3=ts_data_3$umi3,umi5=ts_data_3$umi5_minor,
                      ref_3=ts_data_3$ms_major,ref_5=ts_data_3$ms_minor,Chimeras=ts_data_3$mg_minor,ts_region=ts_data_3$mg_minor)
S8_20 <- rbind(S8_20_1,S8_20_2)

#25
load("./umi_ts_25.RData")
ts_data <- umi_ts_25
ts_data <- ts_data[ts_data$ts==1,]
ts_data$index <- seq(1,length(ts_data[,1]))
ts_data_fil <- ts_data[ts_data$ms_major!=ts_data$ms_minor,]
load("./25_mix/ts_location_10_2.RData")
ts_location$index <- ts_data_fil$index
ts_location <- ts_location[order(ts_location$index,decreasing = F),]
sum(ts_location$index == ts_data_fil$index)
ts_location_part <- ts_location[ts_location$ts_length>1,]
ts_data_fil <- ts_data_fil[ts_data_fil$index %in% ts_location_part$index == T,]
S8_25_1 <- data.frame(cycle=25,umi3=ts_data_fil$umi3,umi5=ts_data_fil$umi5_minor,
                    ref_3=ts_data_fil$ms_major,ref_5=ts_data_fil$ms_minor,Chimeras=ts_data_fil$mg_minor,ts_region=ts_location_part$change_seq)
ts_data_3 <- ts_data[ts_data$ms_major==ts_data$ms_minor,]
S8_25_2 <- data.frame(cycle=25,umi3=ts_data_3$umi3,umi5=ts_data_3$umi5_minor,
                      ref_3=ts_data_3$ms_major,ref_5=ts_data_3$ms_minor,Chimeras=ts_data_3$mg_minor,ts_region=ts_data_3$mg_minor)
S8_25 <- rbind(S8_25_1,S8_25_2)

#30
load("./umi_ts_30.RData")
ts_data <- umi_ts_30
ts_data <- ts_data[ts_data$ts==1,]
ts_data$index <- seq(1,length(ts_data[,1]))
ts_data_fil <- ts_data[ts_data$ms_major!=ts_data$ms_minor,]
load("./30_mix/ts_location_10_2.RData")
ts_location$index <- ts_data_fil$index
ts_location <- ts_location[order(ts_location$index,decreasing = F),]
sum(ts_location$index == ts_data_fil$index)
ts_location_part <- ts_location[ts_location$ts_length>1,]
ts_data_fil <- ts_data_fil[ts_data_fil$index %in% ts_location_part$index == T,]
S8_30_1 <- data.frame(cycle=30,umi3=ts_data_fil$umi3,umi5=ts_data_fil$umi5_minor,
                    ref_3=ts_data_fil$ms_major,ref_5=ts_data_fil$ms_minor,Chimeras=ts_data_fil$mg_minor,ts_region=ts_location_part$change_seq)
ts_data_3 <- ts_data[ts_data$ms_major==ts_data$ms_minor,]
S8_30_2 <- data.frame(cycle=30,umi3=ts_data_3$umi3,umi5=ts_data_3$umi5_minor,
                      ref_3=ts_data_3$ms_major,ref_5=ts_data_3$ms_minor,Chimeras=ts_data_3$mg_minor,ts_region=ts_data_3$mg_minor)
S8_30 <- rbind(S8_30_1,S8_30_2)


S8 <- rbind(S8_20,S8_25,S8_30)

Ref <- read.fasta("./all_16S_ref_PCR.fasta")
Ref_species <- data.frame(Ref=Ref$seq.name,species=0)
Ref_species$species <- c("L. monocytogenes_16S_1","L. monocytogenes_16S_2","L. monocytogenes_16S_3",
                         "P. aeruginosa_16S_1","B. subtilis_16S_1","B. subtilis_16S_2",
                         "B. subtilis_16S_4","B. subtilis_16S_5","B. subtilis_16S_6",
                         "B. subtilis_16S_7","E. coli_16S_1","E. coli_16S_2",
                         "E. coli_16S_4","S. enterica_16S_1","S. enterica_16S_2",
                         "L. fermentum_16S_1","L. fermentum_16S_3","L. fermentum_16S_5",
                         "E. faecalis_16S_1","S. aureus_16S_1","S. aureus_16S_2","S. aureus_16S_3")
for(i in seq(1,nrow(Ref_species)))
{
  S8$ref_3 <- gsub(Ref_species$Ref[[i]], Ref_species$species[[i]], S8$ref_3)
  S8$ref_5 <- gsub(Ref_species$Ref[[i]], Ref_species$species[[i]], S8$ref_5)
}


write.csv(S8,"S8.csv",quote = F,row.names = F)
save(S8,file = "16S_2reads_minor_identity_S8.RData")
