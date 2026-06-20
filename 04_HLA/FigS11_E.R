library(plyr)
library(dplyr)
library(stringr)
library(ggbreak)
library(readxl)
library(ggplot2)
library(patchwork)



my_theme <- theme(legend.position = "bottom",
                  plot.title=element_text(hjust=0.5,size=8),
                  axis.text = element_text(size = 6),
                  strip.text = element_text(size = 8),
                  axis.title = element_text(size = 8),
                  legend.title = element_text(size = 8),
                  legend.text = element_text(size = 6)) 




#S11_E:##################



read_and_combine_txt_files <- function(folder_path) {

  txt_files <- list.files(path = folder_path, 
                          pattern = "\\.txt$", 
                          full.names = TRUE)
  
  data_list <- list()

  for (file in txt_files) {

    file_name <- basename(file)
    
    cycles <- str_extract(file_name, "(?<=_)\\d+(?=_)")
    
    bio <- str_extract(file_name, "\\d+(?=\\.txt)")
    
    if (is.na(cycles) || is.na(bio)) {
      warning(paste("Could not extract cycles or bio from filename:", file_name))
      next
    }
    

    df <- tryCatch({
      read.table(file) %>% 
        mutate(cycles = cycles, bio = bio)
    }, error = function(e) {
      warning(paste("Error reading file:", file, "\n", e$message))
      NULL
    })
    
    if (!is.null(df)) {
      data_list[[file_name]] <- df
    }
  }
  

  if (length(data_list) > 0) {
    combined_data <- bind_rows(data_list)
    return(combined_data)
  } else {
    stop("No valid data could be read from the files.")
  }
}



folder_path <- "/2.split_index_pass3/splited_txt_file_pass3_only_C_gene"
combined_data <- read_and_combine_txt_files(folder_path)



combined_data <- combined_data %>%
  dplyr::rename(
    "zmw_hole" = "V1",
    "index" = "V2",
    "geno" = "V3",
    "umi5" = "V4",
    "umi3" = "V5"
  )%>%
  select(c(umi3,umi5,cycles,bio))


all_reads_result <- combined_data %>%

  group_by(umi3, umi5,cycles,bio) %>%
  summarise(
    umi5_count = n(),  
    .groups = "drop"
  ) %>%

  group_by(umi3, cycles,bio) %>%
  mutate(
    umi3_count = sum(umi5_count) 
  ) %>%
  ungroup() 



#筛选reads>=2的umi3
umi3_reads_filter <- all_reads_result%>%filter(.,umi3_count>=2)






load_and_combine_F1_bio <- function(base_path,bior) {
  
  file_var_mapping <- list(
    "umi_ts_20_F1.RData" = "umi_ts_20",
    "umi_ts_25_F1.RData" = "umi_ts_25",
    "umi_ts_30_F1.RData" = "umi_ts_30"
  )
  

  data_list <- list()
  

  for (file in names(file_var_mapping)) {
    file_path <- file.path(base_path, file)
    var_name <- file_var_mapping[[file]]
    
    if (file.exists(file_path)) {

      env <- new.env()
      load(file_path, envir = env)
      

      cycle <- gsub("umi_ts_(\\d+)_F\\d+\\.RData", "\\1", file)  
      df <- get(var_name, envir = env) %>%
        mutate(cycles = cycle,bio=bior)
      

      data_list[[file]] <- df
    } else {
      warning(paste("File not found:", file_path))
    }
  }
  

  if (length(data_list) == 3) {
    F1_comb<- bind_rows(data_list)
    
    return(F1_comb)
  } else {
    stop("No valid data files were loaded")
  }
}


F1_bio1_data <- load_and_combine_F1_bio("/2reads_minor_identity/1","1")

F1_bio2_data <- load_and_combine_F1_bio("/2reads_minor_identity/2","2")

F1_bio_all <- rbind(F1_bio1_data,F1_bio2_data)%>%select(c(umi3,umi5,type,cycles,bio))


merged_df <- F1_bio_all %>%
  left_join(
    umi3_reads_filter,
    by = c("umi3", "umi5", "cycles","bio")  
  )



all_5umi_reads_sum <- merged_df %>%
  group_by(cycles, bio) %>%
  summarise(umi5_reads_sum = sum(umi5_count)) %>%
  ungroup()

minor_5umi_reads_sum <- merged_df %>%
  filter(.,type=="minor")%>%
  group_by(cycles, bio) %>%
  summarise(minor_umi5_reads_sum = sum(umi5_count)) %>%
  ungroup()



ratio_data <- all_5umi_reads_sum%>%
  left_join(
    minor_5umi_reads_sum,
    by = c("cycles","bio")  
  )%>%
  mutate(.,Ratio=minor_umi5_reads_sum/umi5_reads_sum,
         major_reads=umi5_reads_sum-minor_umi5_reads_sum)

write.csv(ratio_data,file = paste0(outpath,"S11_e_data.csv"),row.names=FALSE)


S11e_data <- ratio_data%>%
  dplyr::group_by(cycles)%>%
  summarise(ratio_mean=mean(Ratio),
            ratio_sd=sd(Ratio),
            sample_size = n(),
            ratio_se=ratio_sd/sqrt(sample_size))

S11e_data$cycles <- factor(S11e_data$cycles)


S11e_f <- ggplot(S11e_data,aes(x=cycles,y=ratio_mean,fill=cycles))+
  geom_bar(stat = "identity",width=0.5)+
  geom_errorbar(aes(ymin = ratio_mean-ratio_se, ymax = ratio_mean+ratio_se,color=cycles), width = 0.3) +
  labs(x=NULL,y="Minor Perc")+ 
  theme_classic()+
  my_theme

ggsave(S11e_f,filename = paste0(outpath,"S11e_f.pdf"),
       width = 8.5,
       height = 8.5,units = "cm")

