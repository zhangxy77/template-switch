setwd("./20_mix/")
library(tidyverse)
library(dplyr)
library(data.table)
library(ggplot2)
library(Biostrings)
library(muscle)
library(phylotools)
library(ape)
library(parallel)
library(dplyr)
Ref <- read.fasta("all_16S_ref_PCR.fasta")
umi_read_1 <- read.table("n16S1_20.txt")
umi_read_2 <- read.table("n16S2_20.txt")
umi_read <- rbind(umi_read_1,umi_read_2)
colnames(umi_read) <- c("hole","index","read","umi_5","umi_3")
load("umi_3_species.RData")
species_frequecy <- read.csv("16S_freq.csv")


umi3_count_2reads <- umi_read %>%
  dplyr::group_by(umi_3)%>%
  dplyr::filter(n()>1)%>%
  dplyr::ungroup()

umi3_valid <- umi3_count_2reads %>%
  group_by(umi_3) %>%
  filter(n_distinct(umi_5) >= 2) %>%
  ungroup()

colnames(umi3_valid)[[3]] <- "Read"

umi3_identity_species <- umi3_valid %>%
  left_join(umi_3_species, by = c("umi_3", "Read"))


umi_3_major <- umi3_identity_species %>%
  group_by(umi_3) %>%
  filter(Identity_Score == max(Identity_Score)) %>%  
  filter(n() == 1) %>%  
  mutate(major_umi5 = umi_5) %>%  
  ungroup()
umi_3_major$umi3_5 <- paste(umi_3_major$umi_3,umi_3_major$umi_5,sep = "-")

umi_3_minor <- umi3_identity_species[umi3_identity_species$umi_3 %in% umi_3_major$umi_3 ==T,]
umi_3_minor$umi3_5 <- paste(umi_3_minor$umi_3,umi_3_minor$umi_5,sep = "-")
umi_3_minor <- umi_3_minor[umi_3_minor$umi3_5 %in% umi_3_major$umi3_5==F,]


umi3_valid_2 <- umi3_count_2reads %>%
  group_by(umi_3) %>%
  filter(n_distinct(umi_5) >= 1) %>%
  ungroup()

colnames(umi3_valid_2)[[3]] <- "Read"

umi3_identity_species_2 <- umi3_valid_2 %>%
  left_join(umi_3_species, by = c("umi_3", "Read"))


umi_3_major_2 <- umi3_identity_species_2 %>%
  group_by(umi_3) %>%
  filter(Identity_Score == max(Identity_Score)) %>%  
  filter(n() == 1) %>%  
  mutate(major_umi5 = umi_5) %>%  
  ungroup()

umi_3_major_2 <- umi_3_major_2[umi_3_major_2$umi_3 %in% umi_3_major$umi_3==F,]
umi_ts_3 <- data.frame(umi3=umi_3_major_2$umi_3,umi5=umi_3_major_2$umi_5,type="major",
                       umi5_major=umi_3_major_2$umi_5,umi5_minor=0,mg_major=umi_3_major_2$Read,mg_minor=0,
                       ref_major=umi_3_major_2$Best_Reference,ref_minor=0,
                       ms_major=umi_3_major_2$Best_species,ms_minor=0,score_major=umi_3_major_2$Identity_Score,score_minor=0,ts=0)
umi_3_minor_unique <- umi_3_minor %>%
  group_by(umi_5, umi_3) %>%
  sample_n(1) %>%  
  ungroup()
umi_3_minor <- umi_3_minor_unique
umi_ts <- data.frame(umi3=rep(0,length(umi_3_minor$umi_3)),umi5=umi_3_minor$umi_5,type="minor",
                     umi5_major=0,umi5_minor=0,mg_major=0,mg_minor=0,ref_major=0,ref_minor=0,
                     ms_major=0,ms_minor=0,score_major=0,score_minor=0,ts=1)
for(i in seq(1,length(umi_3_minor$umi_3)))
{
  umi3 <- umi_3_minor$umi_3[[i]]
  umi_ts$umi3[[i]] <- umi3
  umi_ts$umi5_minor[[i]] <- umi_3_minor$umi_5[[i]]
  umi_ts$mg_minor[[i]] <- umi_3_minor$Read[[i]]
  umi_ts$ref_minor[[i]] <- umi_3_minor$Best_Reference[[i]]
  umi_ts$ms_minor[[i]] <- umi_3_minor$Best_species[[i]]
  umi_ts$score_minor[[i]] <- umi_3_minor$Identity_Score[[i]]
  loc_3 <- which(umi_3_major$umi_3==umi3)
  umi_ts$umi5_major[[i]] <- umi_3_major$umi_5[[loc_3]]
  umi_ts$mg_major[[i]] <- umi_3_major$Read[[loc_3]]
  umi_ts$ref_major[[i]] <- umi_3_major$Best_Reference[[loc_3]]
  umi_ts$ms_major[[i]] <- umi_3_major$Best_species[[loc_3]]
  umi_ts$score_major[[i]] <- umi_3_major$Identity_Score[[loc_3]]
  if(i %% 1000 ==0)
  {save(umi_ts,file = "F123_16S2_20_mix.RData");cat(paste("umi_ts:",i,"/",length(umi_ts[,1]),"\n",sep = ""))}
}
umi_ts_1 <- umi_ts

umi_ts_2 <- data.frame(umi3=umi_3_major$umi_3,umi5=umi_3_major$umi_5,type="major",
                       umi5_major=umi_3_major$umi_5,umi5_minor=0,mg_major=umi_3_major$Read,mg_minor=0,
                       ref_major=umi_3_major$Best_Reference,ref_minor=0,
                       ms_major=umi_3_major$Best_species,ms_minor=0,score_major=umi_3_major$Identity_Score,score_minor=0,ts=0)
umi_ts <- rbind(umi_ts_1,umi_ts_2,umi_ts_3)

umi_ts_20 <- umi_ts
save(umi_ts_20,file = "umi_ts_20.RData")
save(umi3_identity_species,file = "umi3_identity_species.RData")