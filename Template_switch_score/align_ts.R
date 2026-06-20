setwd("./20_mix/muscle_align/")
library(tidyverse)
library(dplyr)
library(data.table)
library(ggplot2)
library(Biostrings)
library(phylotools)
library(splines)  
load("./20_mix/ts_Data.RData")
load("./umi_ts_20.RData")
umi_ts_1 <- umi_ts_20[umi_ts_20$ts==1,]
ts_data$index <- seq(1,length(ts_data[,1]))
ts_data_fil <- ts_data[ts_data$ms_major!=ts_data$ms_minor,]
ts_data <- ts_data_fil

ts_data$key14 <- apply(ts_data[, 1:14], 1, paste, collapse = "_")
umi_ts_1$key14 <- apply(umi_ts_1[, 1:14], 1, paste, collapse = "_")


umi_keys <- unique(umi_ts_1$key14)


ts_data_tagged <- ts_data %>%
  filter(key14 %in% umi_keys) %>%
  distinct(key14, .keep_all = TRUE) %>%
  mutate(tag = 1)
ts_data_fil_2 <- ts_data_tagged[ts_data_tagged$ms_major!=ts_data_tagged$ms_minor,]


load("./Ref_identity.RData")
identity_scores <- as.data.frame(identity_scores) %>% 
  mutate(across(everything(), ~replace_na(.x, 100)))
Ref <- read.fasta("./all_16S_ref_PCR.fasta")
Ref$seq_length <- nchar(Ref$seq.text)
Ref_len <- unlist(strsplit(Ref[1,2],split=""))
ts_data <- ts_data_fil_2
ts_location <- data.frame(hole=ts_data$umi3,read=ts_data$mg_minor,ms_5=ts_data$ms_major,ms_3=ts_data$ms_minor,
                          Q_sequence=0,A_sequence=0,B_sequence=0,change_index=0,Ref_ABchange_index=0,change_seq=0,ts_index=0)
for(i in seq(1,nrow(ts_data)))
{
  if(i %%100 == 0)
    {cat(paste("align:",i,"/",nrow(ts_data),"\n",sep = ""))}
  filename <- paste(ts_data$index[[i]],"_aligned.fasta",sep="")
  aligndata <- read.fasta(filename)
  A_loc <- which(grepl("_A$", aligndata$seq.name))
  B_loc <- which(grepl("_B$", aligndata$seq.name))
  A_sequence <- aligndata$seq.text[[A_loc]]
  B_sequence <- aligndata$seq.text[[B_loc]]
  Q_loc <- which(c(1,2,3) %in% c(A_loc,B_loc)==F)
  Q_sequence <- aligndata$seq.text[[Q_loc]]
  ts_location$Q_sequence[[i]] <- Q_sequence
  ts_location$A_sequence[[i]] <- A_sequence
  ts_location$B_sequence[[i]] <- B_sequence
  Ref_A <- unlist(strsplit(ts_data$ref_major[[i]],split=""))
  Ref_B <- unlist(strsplit(ts_data$ref_minor[[i]],split=""))
  ts_location$Ref_ABchange_index[i] <- list(which(Ref_A != Ref_B))
  
  window_size <- 10
  
  
  start_positions <- seq(1, nchar(Q_sequence) - window_size + 1, by = 1)
  
  
  windows <- tibble(
    window_number = seq_along(start_positions),
    Q_window = map_chr(start_positions, ~ substr(Q_sequence, .x, .x + window_size - 1)),
    A_window = map_chr(start_positions, ~ substr(A_sequence, .x, .x + window_size - 1)),
    B_window = map_chr(start_positions, ~ substr(B_sequence, .x, .x + window_size - 1)),
  )
  
  windows <- as.data.frame(windows)
  
  
  calculate_score <- function(seq1, seq2) {
    
    chars1 <- str_split(seq1, "", simplify = TRUE)
    chars2 <- str_split(seq2, "", simplify = TRUE)
    
    
    score <- sum(ifelse(chars1 == chars2, 1, ifelse(chars1 == "-" | chars2 == "-", -2, -1)))
    return(score)
  }
  
  
  windows <- windows %>%
    mutate(
      QA_score = map2_dbl(Q_window, A_window, calculate_score),
      QB_score = map2_dbl(Q_window, B_window, calculate_score)
    )
  windows$delta_S <- windows$QA_score-windows$QB_score
  change_index <- which(diff(sign(windows$delta_S)) != 0)
  change_index_down <- which(diff(sign(windows$delta_S)) < 0)
  ts_location$change_index[i] <- list(change_index)
  Ms_5 <- sapply(strsplit(ts_data$ms_major[[i]], split = "_"), function(x) paste(x[1:2], collapse = "_"))
  Ms_3 <- sapply(strsplit(ts_data$ms_minor[[i]], split = "_"), function(x) paste(x[1:2], collapse = "_"))
  if(length(change_index_down)!=0)
  {
    if(all(windows$delta_S[1:10] >= 0))
    {
      if(sum(windows$delta_S>0)>0&sum(windows$delta_S<0)>0)
      {
        for(j in seq(1,length(change_index_down)))
        {
          k <- change_index_down[[j]]
          if(j < length(change_index_down))
          {
            d <- change_index_down[[j+1]]
            A_rate <- sum(windows$delta_S[1:k]<0)/length(windows$delta_S[1:k])
            B_rate <- sum(windows$delta_S[k:length(windows$delta_S)]>0)/length(windows$delta_S[k:length(windows$delta_S)])
            if(A_rate < 0.05 &B_rate < 0.05) 
            {
              ts_location$change_seq[[i]] <- substr(Q_sequence, k, d);ts_location$ts_index[[i]]<- paste(k,d,sep = ":");break}
          }
        }
        
      }
      if(sum(windows$delta_S<0)==0)
      {
        len <- length(change_index_down)
        k <- change_index_down[[len]]
        ts_location$change_seq[[i]] <- substr(Q_sequence, k, length(windows$delta_S))
        ts_location$ts_index[[i]]<- paste(k,length(windows$delta_S),sep = ":")
      }
      if(sum(windows$delta_S>0)==0)
      {
        k <- change_index_down[[1]]
        ts_location$change_seq[[i]] <- substr(Q_sequence, 1, k)
        ts_location$ts_index[[i]]<- paste(1,k,sep = ":")
      }
      
    }
  }
      
}
ts_location$ts_length <- nchar(ts_location$change_seq)
save(ts_location,file="./20_mix/ts_location_10_2.RData")
save.image(file="./20_mix/align_2.RData")
