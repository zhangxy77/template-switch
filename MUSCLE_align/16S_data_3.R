setwd("./20_mix/")
library(tidyverse)
library(muscle)
library(phylotools)
library(ape)
library(Biostrings)
library(parallel)
muscle_path <- "./miniconda3/envs/rna/bin/muscle"
Ref <- read.fasta("all_16S_ref_PCR.fasta")
umi_read_1 <- read.table("n16S1_20.txt")
umi_read_2 <- read.table("n16S2_20.txt")
umi_read <- rbind(umi_read_1,umi_read_2)
colnames(umi_read) <- c("hole","index","read","umi_5","umi_3")

#use 
umi3_count_2reads <- umi_read %>%
  dplyr::group_by(umi_3)%>%
  dplyr::filter(n()>1)%>%
  dplyr::ungroup()
umi3_major <- data.frame(umi_3=umi3_count_2reads$umi_3,majorgenotype=umi3_count_2reads$read)
umi3_major <- unique(umi3_major)

save(umi3_major,file="umi3_major.RData")



read_sequences <- as.matrix(umi3_major[,2])
reference_sequences <- Ref[,2]  



identity_score <- function(alignment, reference_seq_length) {
  aln <- as.data.frame(as.character(alignment))
  aln <- as.data.frame(do.call(rbind, lapply(aln[,1], function(x) unlist(strsplit(x, "")))))
  matches <- sum(aln[1,] == aln[2,])
  return(matches / reference_seq_length * 100)  
}


output_dir <- "./umi3_ref_20"
dir.create(output_dir, showWarnings = FALSE)  


process_read <- function(i) {
  tryCatch({
    best_score <- 0
    best_reference <- ""
    best_species <- ""
    
    read_seq <- as.character(read_sequences[i, 1])  
    
    for (j in 1:length(reference_sequences)) {
      ref_seq <- as.character(reference_sequences[j])  
      
      
      temp_in_file <- file.path(output_dir, paste0("read_", i, "_ref_", j, ".fasta"))
      temp_out_file <- file.path(output_dir, paste0("alignment_", i, "_ref_", j, ".fasta"))
      
      
      writeXStringSet(DNAStringSet(c(read_seq, ref_seq)), temp_in_file)
      
      
      system(paste(muscle_path, 
                   "-in", temp_in_file, 
                   "-out", temp_out_file, 
                   "-maxiters 1 -diags -quiet"))
      
      
      alignment <- readDNAStringSet(temp_out_file, format = "fasta")
      
      
      score_1 <- identity_score(alignment, nchar(ref_seq))
      score_2 <- identity_score(alignment, nchar(read_seq))
      score <- min(score_1,score_2)
      
      if (score > best_score) {
        best_score <- score
        best_reference <- reference_sequences[j]
        best_species <- Ref[j, 1]  
      }
      
      
      unlink(temp_in_file) 
      unlink(temp_out_file) 
    }
    
    
    return(list(read = read_sequences[i], best_ref = best_reference, best_spec = best_species, best_score = best_score))
    
  }, error = function(e) {
    
    message("Error in processing read ", i, ": ", e$message)
    message("Read sequence: ", as.character(read_sequences[i, 1]))
    return(list(read = read_sequences[i], best_ref = NA, best_spec = NA, best_score = NA))  
  })
}




umi_3_species <- data.frame(
  umi_3 = umi3_major$umi_3,
  Read = umi3_major$majorgenotype,
  Best_Reference = NA,
  Best_species = NA,
  Identity_Score = NA,
  stringsAsFactors = FALSE
)


total_sequences <- length(read_sequences)
load("umi_3_species.RData")

chunk_size <- 100 
for (start in seq(1, total_sequences, by = chunk_size)) {
  end <- min(start + chunk_size - 1, total_sequences)  
  
  
  result_list <- mclapply(start:end, process_read, mc.cores = 30, mc.preschedule = FALSE)
  
  
  for (i in seq_along(result_list)) {
    idx <- start + i - 1 
    umi_3_species[idx, "Best_Reference"] <- result_list[[i]]$best_ref
    umi_3_species[idx, "Best_species"] <- result_list[[i]]$best_spec
    umi_3_species[idx, "Identity_Score"] <- result_list[[i]]$best_score
  }
  
  
  save(umi_3_species, file = paste0("umi_3_species", ".RData"))
  cat(paste("Saved umi_3_species up to sequence:",end,"/",total_sequences,"\n",sep = ""))
}

