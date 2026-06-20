
library(parallel)

needle_fun <- function(file_path){

  
  setwd(file_path)
  files <- list.files(pattern = "^my")
  
  
  dir.create(paste(file_path,"log",sep = ""),recursive = TRUE)
  dir.create(paste(file_path,"result_sam",sep = ""),recursive = TRUE)
  
  
  mclapply(files,function(i){
    name = sub(".fasta","",i)
    system(paste('cd ~/anaconda3/bin/;source ./activate;nohup time needle ',file_path,'needle_ref_804ura.fasta ',file_path,i,
    ' -gapopen 10 -gapextend 0.5 -aformat sam -outfile ',file_path,'result_sam/',name,'.sam > ',file_path,'log/',name,'.log 2>&1 &', sep = ""))
    
  },mc.cores = 50)

}

run <- needle_fun("/mnt/data6/disk/hxx/Moleculars_Labeled/TUY_D0/2.needle/")
