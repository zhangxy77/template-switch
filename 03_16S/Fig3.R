library(stringr)
library(ggplot2)
library(plyr)
library(readxl)
library(dplyr)
library(ggrepel)
library(patchwork)
library(reshape2)
library(ggsignif)
library(ggbreak)
library(forcats)
library(RColorBrewer)
library(splines)



my_theme <- theme(
  plot.title=element_text(hjust=0.5,size=8),
  axis.text = element_text(size = 6),
  strip.text = element_text(size = 8),
  axis.title = element_text(size = 8),
  legend.title = element_text(size = 8),
  legend.text = element_text(size = 6)) 


#Fig3#######################
load_and_combine_F1 <- function(base_path) {
  

  target_files <- c(
    "umi_ts_20.RData",
    "umi_ts_25.RData",
    "umi_ts_30.RData"
  )
  

  data_list <- list()
  

  for (file in target_files) {
    file_path <- file.path(base_path, file)
    
    if (file.exists(file_path)) {

      env <- new.env()
      load(file_path, envir = env)
      

      cycle <- gsub(".*_|\\..*", "", file) 
      df <- env$umi_ts %>%
        mutate(cycle = cycle)
      

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


F1_combined <- load_and_combine_F1("/3reads_minor/combine")



F1_major <- F1_combined%>%
  filter(.,type=="major")


F1_data<- F1_major%>%
  group_by(cycle,species) %>% 
  summarise(
    count = n(),
    .groups = "drop"                           
  )




standard_freq_data <- read_excel("16S_Combind_data_2_reads.xlsx", sheet = 5)%>%
  mutate(.,frequency2=frequency,
         short_name=str_replace(species,"^([A-Za-z])[a-z]+_([A-Za-z])[a-z]+.*", "\\1\\2") %>%toupper())%>%
  subset(.,select=c(-frequency))%>%
  arrange(frequency2)

F1_factor <- factor(standard_freq_data$short_name,levels = unique(standard_freq_data$short_name))


F1_data2 <- F1_data%>%
  group_by(cycle) %>%  
  dplyr::group_modify(~ {  
    merge(.x, standard_freq_data, 
          by = "species",
          all = TRUE) 
  })%>%
  ungroup()%>% 
  mutate(
    count = ifelse(is.na(count), 0, count)
  )


F1_data3 <- F1_data2%>%
  group_by(cycle) %>%  
  mutate(cycle_count = sum(count)) %>% 
  ungroup()%>% 
  mutate(.,F1=count/cycle_count,
         F1_perc=F1*100)

F1_data3$short_name <- factor(F1_data3$short_name,levels=unique(F1_factor))
F1_data3$cycle <- factor(F1_data3$cycle)


#Fig3_C##################
F1_20_cor <- F1_data3%>%
  filter(cycle==20)%>%
  select(short_name,F1_perc,frequency2)%>%
  filter(.,F1_perc!=0)%>%
  arrange(frequency2,F1_perc)%>%
  mutate(.,type=rep(c("down", "up"),length.out = n()),
         combined_legend = paste(short_name, type, sep = " - "))


F1_20_cor_result <- cor.test(F1_20_cor$frequency2,F1_20_cor$F1_perc, method = "pearson")
F1_20_cor_coef <- F1_20_cor_result$estimate  
F1_20_p_value <- F1_20_cor_result$p.value  


F1_20_cor$short_name <-factor(F1_20_cor$short_name ,levels=unique(F1_factor))

F1_20_cor$type <-factor(F1_20_cor$type)


F1_20_cor_pic2 <- ggplot(F1_20_cor,aes(x=frequency2,y=F1_perc))+
  geom_abline(intercept = 0,slope = 1,col="grey",lty = 2)+
  geom_point(aes(fill = short_name, shape = type), size = 2, color = "black")+
  scale_x_log10(limits=c(0.0001,100),breaks = c(0.0001,0.001,0.01,0.1,1,10,100),labels = c("0.0001","0.001","0.01","0.1","1","10","100"))+#通过breaks和labels参数手动设置了刻度的位置和标签，以禁用科学计数法并以指定格式显示坐标轴上的数值。
  scale_y_log10(limits=c(0.0001,100),breaks = c(0.0001,0.001,0.01,0.1,1,10,100),labels = c("0.0001","0.001","0.01","0.1","1","10","100"))+
  scale_fill_brewer(palette = "Set3")+
  scale_shape_manual(
    values = c("up" = 25, "down" = 24),
    name = "Species"
  ) +
  guides(fill = guide_legend(override.aes = list(shape = c(24, 25)[as.numeric(factor(F1_20_cor$type))])),
         shape = "none")+ 
  annotate("text", x = 0.0001, y = 100,
           label =  paste("R == ",ifelse(F1_20_cor_coef>0.99,"0.99",F1_20_cor_coef), 
                          " * ', ' ~ italic(P) == '", format(signif(F1_20_p_value, digits = 3), scientific = TRUE),"'"), 
           hjust = 0, vjust = 1, size = 2, parse = TRUE)+
  labs(
    x="Expected frequency (%)",y="Observed frequency (%)",
    fill = "Species")+
  theme_classic()+
  my_theme+
  theme(legend.text = element_text(face = "italic"))



ggsave(F1_20_cor_pic2,filename = paste0(output,"triangle_16S_F1_20cycles_3reads_cor_Frequency_minor_v2.pdf"), 
       width = 8.5,
       height = 6,units = "cm")




data_Fig3_C <- F1_20_cor%>%
  select(short_name,F1_perc,frequency2)

colnames(data_Fig3_C) <- c("Species","Observed Frequency","Expected Frequency")

write.csv(data_Fig3_C,file=paste0(output,"data_Fig3_C.csv"),row.names = F)





#FigS9##############################


load_and_combine_F1_bio <- function(base_path) {
  

  file_var_mapping <- list(
    "umi_ts_20.RData" = "umi_ts_20",
    "umi_ts_25.RData" = "umi_ts_25",
    "umi_ts_30.RData" = "umi_ts_30"
  )
  

  data_list <- list()
  

  for (file in names(file_var_mapping)) {
    file_path <- file.path(base_path, file)
    var_name <- file_var_mapping[[file]]
    
    if (file.exists(file_path)) {

      env <- new.env()
      load(file_path, envir = env)
      

      cycle <- gsub(".*_|\\..*", "", file) 
      df <- get(var_name, envir = env) %>%
        mutate(cycle = cycle)
      

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


F1_bio1_data <- load_and_combine_F1_bio("/16S2/3reads_minor/1")

F1_bio1_data2 <- F1_bio1_data%>%filter(.,type=="major")%>%mutate(.,biorepeat="repeat1")


F1_bio2_data <- load_and_combine_F1_bio("/16S2/3reads_minor/2")

F1_bio2_data2 <- F1_bio2_data%>%filter(.,type=="major")%>%mutate(.,biorepeat="repeat2")

F1_bio_all <- rbind(F1_bio1_data2,F1_bio2_data2)


F1_bio_all2 <- F1_bio_all%>%
  group_by(biorepeat,cycle,species) %>%  
  summarise(
    count = n(),
    .groups = "drop"                             
  )



F1_bio_all3 <- F1_bio_all2%>%
  group_by(biorepeat,cycle) %>% 
  mutate(bio_cycle_count = sum(count)) %>%  
  ungroup()%>%  
  mutate(.,F1=count/bio_cycle_count,
         F1_perc=F1*100)



F1_bio_all4 <- reshape2::dcast(F1_bio_all3,cycle+species~biorepeat,value.var = "F1_perc")%>%
  mutate(.,short_name=str_replace(species,"^([A-Za-z])[a-z]+_([A-Za-z])[a-z]+.*", "\\1\\2") %>%toupper())


F1_bio_all5 <- F1_bio_all4%>%
  left_join(F1_20_cor %>% select(short_name, type), by = "short_name")%>%
  filter(.,repeat1>0 &repeat2>0)




cor_results <- F1_bio_all5 %>%
  group_by(cycle) %>%
  summarise(
    R = cor.test(repeat1, repeat2)$estimate,
    p_value = cor.test(repeat1, repeat2)$p.value,
    .groups = "drop"
  ) %>%
  mutate(label = paste("italic(R) == ", ifelse(R>0.99,"0.99",R), 
                       " * ', ' ~ italic(P) == '", format(signif(p_value, digits = 3), scientific = TRUE),"'"))





F1_bio_all5$short_name <- factor(F1_bio_all5$short_name, levels = unique(F1_factor))
F1_bio_all5$type <- factor(F1_bio_all5$type, levels = c("up", "down"))


set3_colors <- brewer.pal(12, "Set3")[2:7]  




F1_bio_pic <- ggplot(data=F1_bio_all5,aes(x=repeat1,y=repeat2))+
  geom_abline(intercept = 0,slope = 1,col="grey",lty = 2)+
  geom_point(aes(fill = short_name, shape = type), size = 2, color = "black")+
  facet_wrap(~cycle,scales = "free",
             labeller = labeller(cycle = function(x) paste0(x, " cycles")))+
  scale_x_log10(limits=c(0.0001,100),breaks = c(0.0001,0.001,0.01,0.1,1,10,100),labels = c("0.0001","0.001","0.01","0.1","1","10","100"))+#通过breaks和labels参数手动设置了刻度的位置和标签，以禁用科学计数法并以指定格式显示坐标轴上的数值。
  scale_y_log10(limits=c(0.0001,100),breaks = c(0.0001,0.001,0.01,0.1,1,10,100),labels = c("0.0001","0.001","0.01","0.1","1","10","100"))+
  scale_fill_manual(values = set3_colors) +
  scale_shape_manual(
    values = c("up" = 25, "down" = 24),
    name = "Species"
  ) +
  guides(
    fill = guide_legend(override.aes = list(shape = ifelse(levels(F1_bio_all5$type) == "up", 25, 24))),  
    shape = "none"  
  )+
  geom_text(data = cor_results, aes(x = 0.0001, y = 100, label = label), 
            size=2,hjust = 0, vjust = 1,parse = TRUE,inherit.aes = FALSE) +
  labs(
    x="Observed frequency in technical replicate1 (%)",
    y="Observed frequency in technical replicate2 (%)",fill="Species")+
  theme_classic()+
  my_theme+
  theme(legend.text = element_text(face = "italic"),
        legend.position = "bottom",
        axis.text.x = element_text(angle = 45,hjust=1))




ggsave(F1_bio_pic,filename = paste0(output,"triangle_16S_F1_3reads_biorepeat_Frequency_minor_v4.pdf"), 
       width = 19,
       height = 8,units="cm") 


#Fig3_H#################################

F1_top3_data <- F1_bio_all3%>%
  group_by(biorepeat,cycle) %>%
  arrange(desc(F1_perc)) %>%  
  slice_head(n = 3) %>%  
  ungroup()%>%  
  select(-count,-bio_cycle_count,-F1)

F1_top3_data_20 <- F1_top3_data%>%
  filter(cycle=="20")%>%
  dplyr::rename(F1_perc_20 = F1_perc)%>%  
  select(-cycle)


F1_top3_data_ratio <-F1_top3_data %>%
  left_join(F1_top3_data_20, by = c("biorepeat","species"))%>%
  mutate(ratio = F1_perc / F1_perc_20) %>%
  select(-F1_perc,-F1_perc_20) 


F1_top3_data_ratio2 <- F1_top3_data_ratio%>%
  group_by(cycle,species)%>%
  summarise(SD=sd(ratio),
            Mean=mean(ratio),
            SE = sd(ratio) / sqrt(n()),   
            N = n(),                      
            .groups = "drop")%>%
  mutate(.,
         short_name=str_replace(species,"^([A-Za-z])[a-z]+_([A-Za-z])[a-z]+.*", "\\1\\2") %>%toupper(),
         Mean_perc=Mean*100,
         SD_perc=SD*100,
         SE_perc = SE * 100)


F1_top3_data_ratio2 <- F1_top3_data_ratio2%>%
  left_join(standard_freq_data,by="short_name")%>%
  arrange(desc(frequency2))

F1_top3_data_ratio2$short_name <- factor(F1_top3_data_ratio2$short_name,level=unique(F1_factor))



top3_p <- ggplot(F1_top3_data_ratio2, aes(x = fct_rev(short_name), y = Mean_perc, fill = factor(cycle))) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.5) +  
  geom_hline(yintercept = 100, linetype = "dashed", color = "grey", linewidth = 0.4) + 
  labs(
    x = NULL,  
    y = "Relative observed frequency (%)",   
    fill = "Cycles"
  ) +
  scale_y_continuous(limits = c(0, 130), breaks = c(0, 50, 100))+
  geom_errorbar(aes(ymin = Mean_perc - SE_perc, ymax = Mean_perc + SE_perc), 
                position = position_dodge(width = 0.8),    
                width = 0.1, color = "black")+
  scale_x_discrete(labels = c(
    "L. monocytogenes(95.9%)",
    "P. aeruginosa(2.8%)",
    "B. subtilis(1.2%)"
  ))+ 
  theme_classic()+
  my_theme+
  theme(legend.position = "none")


t_test_results <- F1_top3_data_ratio %>%
  group_by(species) %>%
  summarise(
    p_value_25 = t.test(ratio[cycle == "25"], ratio[cycle == "20"])$p.value,  # 25 vs 20
    p_value_30 = t.test(ratio[cycle == "20"], ratio[cycle == "30"])$p.value,   # 20 vs 30
    .groups = "drop"
  ) %>%
  mutate(
    short_name=str_replace(species,"^([A-Za-z])[a-z]+_([A-Za-z])[a-z]+.*", "\\1\\2") %>%toupper(),
    label_25 = paste0("italic(P)==",format(signif(p_value_25, digits = 3), scientific = TRUE)), 
    label_30 = paste0("italic(P)==",format(signif(p_value_30, digits = 3), scientific = TRUE))
  )

t_test_results2 <- t_test_results%>%
  mutate(.,
         p_value_25_onetail=p_value_25/2,
         p_value_30_onetail=p_value_30/2,
         label_25_onetail = paste0("italic(P)=='",format(signif(p_value_25_onetail, digits = 3), scientific = TRUE),"'"),
         label_30_onetail = paste0("italic(P)=='",format(signif(p_value_30_onetail, digits = 3), scientific = TRUE),"'"))




species_order <- unique(F1_top3_data_ratio2$short_name)



for (i in seq_along(t_test_results2$short_name)) {
  species <- t_test_results2$short_name[i]
  species_index <- which(species_order == species)  
  
  top3_p <- top3_p + geom_signif(
    annotations = t_test_results2$label_25_onetail[i],
    y_position = max(F1_top3_data_ratio2$Mean_perc) * 1.1,  
    xmin = species_index - 0.26, 
    xmax = species_index,
    tip_length = 0.01,
    vjust = 0,
    textsize = 2,
    parse = TRUE   
  )
  
  top3_p <- top3_p + geom_signif(
    annotations =t_test_results2$label_30_onetail[i],
    y_position = max(F1_top3_data_ratio2$Mean_perc) * 1.2, 
    xmin = species_index - 0.26 ,  
    xmax = species_index + 0.26,
    tip_length = 0.01,
    vjust = 0,
    textsize = 2,
    parse = TRUE   
  )
}


ggsave(top3_p,filename = paste0(output,"se_16S_F1_top3_ratio_3reads_minor_onetail_v2.pdf"), 
       width = 9.5,
       height = 4.5,units = "cm") 



data_Fig3_H <- F1_top3_data_ratio%>%
  mutate(.,species=str_replace(species, "^([A-Z])[a-z]+_([a-z]+)", "\\1. \\2"),
         ratio=ratio*100)

colnames(data_Fig3_H) <- c("Repeat","Cycles","Species","Relative observed Frequency")

write.csv(data_Fig3_H,file=paste0(output,"data_Fig3_H.csv"),row.names = F)




#Fig3_D##########################
load("16S_3reads_minor_S8.RData")

S8_2 <- S8%>%
  subset(.,select=c(-Chimeras,-ts_region))%>%
  mutate(species_3 = str_split(ref_3, "_", simplify = TRUE)[, 1],
         species_5 = str_split(ref_5, "_", simplify = TRUE)[, 1])


comp_data_by_group <- S8_2 %>%
  group_by(cycle) %>%  
  summarise(
    same16S = sum(ref_3 == ref_5),               
    samesp = sum((species_3 == species_5)& (ref_3 != ref_5)),        
    diffsp = sum(species_3 != species_5),         
    all_events = same16S+samesp+diffsp,           
    .groups = "drop"                             
  ) %>%
  mutate(
    fre_same16S = same16S / all_events,           
    fre_samesp = samesp / all_events,            
    fre_diffsp = diffsp / all_events              
  )


comp_data_long <- comp_data_by_group %>%
  pivot_longer(
    cols = c(same16S, samesp, diffsp),
    names_to = "category",
    values_to = "value"
  ) %>%
  mutate(
    fre = value / all_events  
  )

comp_data_long$category <- factor(comp_data_long$category,levels = c("same16S","samesp","diffsp"))


pie_plot <- ggplot(comp_data_long, aes(x = "", y = fre, fill = category)) +
  geom_bar(stat = "identity", width = 1, alpha = 0.9) +
  coord_polar("y", start = 0) +
  geom_text(
    aes(label = paste0(value, "\n", sprintf("%.2f%%", fre * 100))),
    position = position_stack(vjust = 0.5),
    size = 4
  ) +
  scale_fill_manual(
    values = c("#619CFF", "#F8766D", "#00BA38"),  
    labels = c("Within the same 16s rRNA", "Between 16s rRNAs within the same species", "Between 16s rRNAs within different species")  
  ) +
  facet_wrap(~ cycle,
             labeller = labeller(cycle = function(x) paste0(x, " cycles"))) +  
  labs(fill = "Category") +
  theme_void() +
  theme(
    strip.text = element_text(size = 15), 
    legend.position = "bottom",             
    legend.title = element_blank(),
    legend.text = element_text(size = 10)
  )


#保存图片
ggsave(pie_plot,filename = paste0(output,"16S_3reads_pie.pdf"), 
       width = 10,
       height = 5) 



data_Fig3_D <- comp_data_by_group%>%
  select(cycle,same16S,samesp,diffsp,all_events)

colnames(data_Fig3_D) <- c("Cycles","Number of same 16s rRNA","Number of same species","Number of different species","All events")

write.csv(data_Fig3_D,file=paste0(output,"data_Fig3_D.csv"),row.names = F)


###############################

load("/3reads_minor/1/16S_3reads_minor_S8.RData")

S8_bio1 <- S8%>%mutate(.,biorepeat="1")


load("/3reads_minor/2/16S_3reads_minor_S8.RData")

S8_bio2 <- S8%>%mutate(.,biorepeat="2")

S8_bio <- rbind(S8_bio1,S8_bio2)%>%
  subset(.,select=c(-Chimeras,-ts_region))%>%
  mutate(species_3 = str_split(ref_3, "_", simplify = TRUE)[, 1],
         species_5 = str_split(ref_5, "_", simplify = TRUE)[, 1])

#Fig3_E###########################
F2_data <- S8_bio%>%
  filter(.,ref_3=="L. monocytogenes_16S_1")


F2_by_group <- F2_data %>%
  group_by(cycle,biorepeat,ref_3,ref_5) %>%  
  summarise(
    count = n(),
    .groups = "drop"                             
  ) %>%
  mutate(ref_5_short_name = str_replace(ref_5, "_16S", ""))




load("/16S2/Ref_identity.RData")
identity_scores_2 <- identity_scores%>%as.data.frame()


identity_scores_2[is.na(identity_scores_2)] <- 100

lis_1_identity <- identity_scores_2%>%
  subset(.,select=c("Listeria_monocytogenes_16S_1"))


lis_1_identity <- lis_1_identity %>% 
  mutate(row_id = row.names(.),
         short_name = str_replace(row_id,
                                  "^([A-Z])[a-z]+_([a-z]+)_16S_(\\d+)",  
                                  "\\1. \\2_\\3"                   
         )
  )%>%
  dplyr::rename(score = Listeria_monocytogenes_16S_1)%>%
  subset(.,select=c(-row_id))

rownames(lis_1_identity) <- NULL 


F2_merged <- F2_by_group %>%
  group_by(cycle,biorepeat) %>%  
  dplyr::group_modify(~ {  
    merge(.x, lis_1_identity, 
          by.x = "ref_5_short_name", 
          by.y = "short_name", 
          all = TRUE)  
  }) %>%
  ungroup()%>%  
  subset(.,select=c(-ref_3,-ref_5))%>%
  mutate(
    count = ifelse(is.na(count), 0, count)
  )


F2_merged2 <- F2_merged%>%
  group_by(cycle,biorepeat) %>%  
  mutate(cycle_biorepeat_count = sum(count)) %>% 
  ungroup()%>%  
  mutate(.,F2=count/cycle_biorepeat_count,
         F2_perc=F2*100)%>%
  group_by(cycle,ref_5_short_name,score) %>%
  summarise(
    SD = sd(F2_perc),                
    Mean = mean(F2_perc),           
    SE = sd(F2_perc) / sqrt(n()),    
    N = n(),                       
    .groups = "drop"              
  )


F2_merged3 <- F2_merged2%>%
  mutate(Label = paste(ref_5_short_name, cycle, sep = ";")) %>%
  arrange(cycle,desc(score),desc(Mean))

F2_merged3$Label <- factor(F2_merged3$Label, levels = unique(F2_merged3$Label))


custom_labeller <- function(x) {
  sapply(x, function(label) {
    unlist(strsplit(label, ";"))[1]
  })
}


custom_labeller2 <- function(x) {
  sapply(x, function(label) {
    part1 <- unlist(strsplit(label, ";"))[1]
    num <- gsub(".*_(\\d+)$", "\\1", part1)
    return(num)
  })
}

F2_plot_pic <-  ggplot(data=F2_merged3,aes(x=Label,y=Mean))+
  geom_blank() +
  geom_bar(stat = "identity",width = 0.5,fill="#FB8072",data = ~ subset(.x, Mean != 0))+ 
  facet_wrap(~cycle,scales = "free",
             labeller = labeller(cycle = function(x) paste0(x, " cycles")))+
  geom_errorbar(aes(ymin = Mean-SE, ymax = Mean + SE), width = 0.1, color = "black") +
  geom_smooth(aes(group = 1), method = "loess",  span = 0.2,method.args = list(degree = 1),se = FALSE, color = "grey", linewidth = 0.5) +
  scale_x_discrete(labels = custom_labeller2)+ 
  scale_y_continuous(breaks = scales::pretty_breaks(n = 5))+
  labs(
    x="Sequence similarity",
    y="Percentage of template switching events (%)")+
  theme_classic()+
  my_theme  




ggsave(F2_plot_pic,filename = paste0(output_se_F23,"se_16S_F2_3reads_reads_minor_v2_noname.pdf"), 
       width = 10,
       height = 5,units = "cm")



data_Fig3_E <- F2_merged%>%
  group_by(cycle,biorepeat) %>%  
  mutate(cycle_biorepeat_count = sum(count)) %>%  
  ungroup()%>%  
  mutate(.,F2=count/cycle_biorepeat_count,
         F2_perc=F2*100)%>%
  mutate(ref_5_short_name=str_replace(ref_5_short_name, 
                                      "^([A-Z])[a-z]*\\.?\\s*([a-z])[a-z]*_([0-9]+)", 
                                      "\\1\\2_\\3") %>%toupper())%>%
  arrange(cycle,biorepeat,desc(score),desc(F2_perc))%>%
  select(cycle,biorepeat,ref_5_short_name,count,cycle_biorepeat_count,F2_perc,score)

colnames(data_Fig3_E) <- c("Cycles","Repeat","16S rRNA","Number of events","Number of all events","Percentage","Sequence similarity")

write.csv(data_Fig3_E,file=paste0(output_se_F23,"data_Fig3_E.csv"),row.names = F)

#Fig3_F#####################
F3_data <- S8_bio%>%
  group_by(cycle,biorepeat,species_3) %>%  
  summarise(
    F3_count = n(),
    .groups = "drop"                              
  ) %>%
  mutate(.,species_3_2=str_replace_all(species_3, "(\\b[A-Za-z])[A-Za-z]*", "\\1") %>%
           str_remove_all("[^A-Za-z]") %>%
           toupper())



F3_merged <- F3_data %>%
  group_by(cycle,biorepeat) %>%  
  dplyr::group_modify(~ {   
    merge(.x%>%select(-species_3), standard_freq_data%>%select(-species), 
          by.x = "species_3_2", 
          by.y = "short_name", 
          all = TRUE) 
  }) %>%
  ungroup()%>%
  mutate(
    F3_count = ifelse(is.na(F3_count), 0, F3_count)
  )


F3_merged2 <- F3_merged%>%
  group_by(cycle,biorepeat) %>%  
  mutate(cycle_biorepeat_count = sum(F3_count)) %>%  
  ungroup()%>%  
  mutate(.,F3=F3_count/cycle_biorepeat_count,
         F3_perc=F3*100)%>%
  group_by(cycle,species_3_2,frequency2) %>%
  summarise(
    SD = sd(F3_perc),              
    Mean = mean(F3_perc),           
    SE = sd(F3_perc) / sqrt(n()),   
    N = n(),                       
    .groups = "drop"             
  )%>%
  arrange(cycle,desc(frequency2),desc(Mean))


F3_merged2$species_3_2 <- factor(F3_merged2$species_3_2,levels = rev(unique(F1_factor)))


#改为means ± s.e.!!

F3_plot_pic2 <- ggplot(data=F3_merged2,aes(x=species_3_2,y=Mean))+
  geom_blank() +
  geom_bar(stat = "identity",width = 0.5,fill="#FB8072",data = ~ subset(.x, Mean != 0))+
  facet_wrap(~cycle,scales = "free",
             labeller = labeller(cycle = function(x) paste0(x, " cycles")))+
  geom_errorbar(aes(ymin = Mean-SE, ymax = Mean + SE), width = 0.1, color = "black") +
  geom_smooth(aes(group = 1), method = "loess", span = 0.7,method.args = list(degree = 2),se = FALSE, color = "grey", linewidth = 0.5) +
  scale_y_continuous(limits = c(0, 100), breaks = c(0, 20, 40, 60, 80,100))+
  labs(
    x="Expected frequency",
    y="Percentage of template switching events (%)")+
  theme_classic()+
  theme(axis.text.x = element_text(face = "italic",angle = 45,hjust=1))+  
  my_theme


ggsave(F3_plot_pic2,filename = paste0(output_se_F23,"se_16S_F3_3reads_minor_v2_standard.pdf"), 
       width = 10,
       height = 5.5,units = "cm") 



data_Fig3_F <- F3_merged%>%
  group_by(cycle,biorepeat) %>% 
  mutate(cycle_biorepeat_count = sum(F3_count)) %>%  
  ungroup()%>%  
  mutate(.,F3=F3_count/cycle_biorepeat_count,
         F3_perc=F3*100)%>%
  arrange(cycle,biorepeat,desc(frequency2),desc(F3_perc))%>%
  select(cycle,biorepeat,species_3_2,F3_count,cycle_biorepeat_count,F3_perc,frequency2)

colnames(data_Fig3_F) <- c("Cycles","Repeat","Species","Number of events","Number of all events","Percentage","Expected frequency")

write.csv(data_Fig3_F,file=paste0(output_se_F23,"data_Fig3_F.csv"),row.names = F)





#Fig3_G#####################

load("/16S2_3reads_minor/plot2.RData")

combined_plot <- ggplot(combined_data, aes(x = sequence_index, y = Count)) +  

  geom_vline(
    data = filter(vline_data, type == "region"),
    aes(xintercept = xintercept),
    linetype = "dashed",
    color = "gray"
  ) +

  geom_vline(
    data = filter(vline_data, type == "Lis_region"),
    aes(xintercept = xintercept),
    linetype = "dashed",
    color = "#FB8072"
  )  +  
  geom_point(size=0.5) + 
  facet_wrap(~cycle, nrow = 1, scales = "free_y") + 
  labs(x = "Position in 16S rRNA sequence (bp)", y = "Template switching score") +  
  theme_classic() +  
  scale_y_continuous(limits = c(0, 220), breaks = c(0, 50, 100, 150, 200)) +  
  my_theme


ggsave(combined_plot,filename = paste0(output,"16S_3reads_ts_score.pdf"), 
       width = 10,
       height = 5.5,units = "cm") 

data_fig3_G1 <- combined_data%>%mutate(.,cycle=str_extract(cycle, "\\d+"))%>%
  select(cycle,sequence_index,Count)
colnames(data_fig3_G1) <- c("Cycles","Position", "Template Switch Score" )

write.csv(data_fig3_G1,file=paste0(output,"data_fig3_G1.csv"),row.names = F)






