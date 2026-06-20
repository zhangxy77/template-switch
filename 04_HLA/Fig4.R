library(RColorBrewer)
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
library(stringr)
library(splines)


my_theme <- theme(
  plot.title=element_text(hjust=0.5,size=8),
  axis.text = element_text(size = 6),
  strip.text = element_text(size = 8),
  axis.title = element_text(size = 8),
  legend.title = element_text(size = 8),
  legend.text = element_text(size = 6)) 


#Fig4####################

load_and_combine_F1_bio <- function(base_path) {
  
  # 定义目标文件名和对应的变量名（即导入数据框的名字）映射
  file_var_mapping <- list(
    "umi_ts_20_F1.RData" = "umi_ts_20",
    "umi_ts_25_F1.RData" = "umi_ts_25",
    "umi_ts_30_F1.RData" = "umi_ts_30"
  )
  
  # 初始化空列表存储数据
  data_list <- list()
  
  # 循环加载每个文件
  for (file in names(file_var_mapping)) {
    file_path <- file.path(base_path, file)
    var_name <- file_var_mapping[[file]]
    
    if (file.exists(file_path)) {
      # 加载数据到临时环境
      env <- new.env()
      load(file_path, envir = env)
      
      # 提取umi_ts数据并添加cycle列
      cycle <- gsub("umi_ts_(\\d+)_F\\d+\\.RData", "\\1", file)  # 从文件名提取数字(20/25/30)
      df <- get(var_name, envir = env) %>%
        mutate(cycle = cycle)
      
      # 添加到列表
      data_list[[file]] <- df
    } else {
      warning(paste("File not found:", file_path))
    }
  }
  
  # 合并所有数据并筛选major类型，三个循环数的数据都存在list之后再合并
  if (length(data_list) == 3) {
    F1_comb<- bind_rows(data_list)
    
    return(F1_comb)
  } else {
    stop("No valid data files were loaded")
  }
}

F1_combined <- load_and_combine_F1_bio("/2reads_minor_identity/combine")



F1_combined2 <- F1_combined%>%
  mutate(species=ms_major)



F1_major <- F1_combined2%>%
  filter(.,type=="major")


F1_data<- F1_major%>%
  group_by(cycle,species) %>% 
  summarise(
    count = n(),
    .groups = "drop" 
  )


standard_freq_data <- read_excel("/HLA/HLA_Combind_data_2_reads.xlsx", sheet = 5)%>%
  mutate(.,frequency2=frequency,
         short_name=gsub(".*(C\\*\\d+:\\d+:\\d+).*", "\\1", species))%>%
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
         F1_perc=F1*100)%>%
  arrange(cycle,frequency2,F1_perc)


F1_data3$short_name <- factor(F1_data3$short_name,levels=unique(F1_factor))
F1_data3$cycle <- factor(F1_data3$cycle)


#Fig4_C##################

F1_20_cor <- F1_data3%>%
  filter(cycle==20)%>%
  select(short_name,F1_perc,frequency2)%>%
  mutate(frequency3=frequency2*100)%>%
  filter(.,F1_perc!=0)




F1_20_cor_result <- cor.test(F1_20_cor$frequency3,F1_20_cor$F1_perc, method = "pearson")
F1_20_cor_coef <- F1_20_cor_result$estimate 
F1_20_p_value <- F1_20_cor_result$p.value 



F1_20_cor$short_name <-factor(F1_20_cor$short_name ,levels=unique(F1_factor))




F1_20_cor <- F1_20_cor[order(F1_20_cor$short_name), ]%>%
  mutate(.,type=rep(c("down", "up"),length.out = n()),
         combined_legend = paste(short_name, type, sep = " - "))


F1_20_cor$type <-factor(F1_20_cor$type)
shape_map <- c("up" = 25, "down" = 24)



F1_20_cor_pic2 <- ggplot(F1_20_cor, aes(x = frequency3, y = F1_perc)) +
  geom_abline(intercept = 0, slope = 1, col = "grey", lty = 2) +
  geom_point(aes(fill = short_name, shape = type), 
             size = 2, color = "black") +
  scale_y_log10(limits=c(0.005,500),breaks = c(0.005,0.05,0.5,5,50,500),labels = c("0.005","0.05","0.5","5","50","500"))+
  scale_x_log10(limits=c(0.005,500),breaks = c(0.005,0.05,0.5,5,50,500),labels = c("0.005","0.05","0.5","5","50","500"))+
  scale_fill_brewer(palette = "Set3", name = "Alleles") +
  scale_shape_manual(
    values = c("up" = 25, "down" = 24),
    name = "Alleles"
  ) +
  guides(
    fill = guide_legend(override.aes = list(shape = shape_map[levels(F1_20_cor$type)])),
    shape = "none")+ 
  annotate("text", x =0.005, y = 500,
           label =  paste("R == ", round(F1_20_cor_coef, 2), 
                          " * ', ' ~ italic(P) == '", format(signif(F1_20_p_value, digits = 3), scientific = TRUE),"'"), 
           hjust = 0, vjust = 1, size = 2, parse = TRUE)+
  labs(
    x="Expected frequency (%)",y="Observed frequency (%)",
    fill = "Allele")+
  theme_classic()+
  my_theme+
  theme(legend.text = element_text(face = "italic"))


ggsave(F1_20_cor_pic2,filename = paste0(output,"triangle_HLA_F1_20cycles_2reads_cor_Frequency_minor_v3.pdf"),
       height = 6,units = "cm")



data_Fig4_C <- F1_20_cor%>%
  select(short_name,F1_perc,frequency3)

colnames(data_Fig4_C) <- c("Alleles","Observed Frequency","Expected Frequency")

write.csv(data_Fig4_C,file=paste0(output,"data_Fig4_C.csv"),row.names = F)






#FigS11_F:##############################




F1_bio1_data <- load_and_combine_F1_bio("/HLA2_C/2reads_minor_identity/1")

F1_bio1_data2 <- F1_bio1_data%>%filter(.,type=="major")%>%mutate(.,biorepeat="repeat1")



F1_bio2_data <- load_and_combine_F1_bio("/HLA2_C/2reads_minor_identity/2")

F1_bio2_data2 <- F1_bio2_data%>%filter(.,type=="major")%>%mutate(.,biorepeat="repeat2")

F1_bio_all <- rbind(F1_bio1_data2,F1_bio2_data2)



F1_bio_all2 <- F1_bio_all%>%mutate(species=ms_major)



F1_bio_all3 <- F1_bio_all2%>%
  group_by(biorepeat,cycle,species) %>% 
  summarise(
    count = n(),
    .groups = "drop"
  )


F1_bio_all4 <- F1_bio_all3%>%
  group_by(biorepeat,cycle) %>%
  mutate(bio_cycle_count = sum(count)) %>%
  ungroup()%>% 
  mutate(.,F1=count/bio_cycle_count,
         F1_perc=F1*100)


F1_bio_all5 <- reshape2::dcast(F1_bio_all4,cycle+species~biorepeat,value.var = "F1_perc")%>%
  mutate(.,short_name=gsub(".*(C\\*\\d+:\\d+:\\d+).*", "\\1", species))


F1_bio_all6 <- F1_bio_all5%>%
  left_join(F1_20_cor %>% select(short_name, type), by = "short_name")%>%
  filter(.,repeat1>0 &repeat2>0)



cor_results <- F1_bio_all6 %>%
  group_by(cycle) %>%
  summarise(
    R = cor.test(repeat1, repeat2)$estimate,
    p_value = cor.test(repeat1, repeat2)$p.value,
    .groups = "drop"
  ) %>%
  mutate(label = paste("italic(R) == ", ifelse(R>0.99,"0.99",R), 
                       " * ', ' ~ italic(P) == '", format(signif(p_value, digits = 3), scientific = TRUE),"'"))



F1_bio_all6$short_name <- factor(F1_bio_all6$short_name, levels = unique(F1_factor))
F1_bio_all6$type <- factor(F1_bio_all6$type)



F1_bio_pic <- ggplot(data=F1_bio_all6,aes(x=repeat1,y=repeat2))+
  geom_abline(intercept = 0,slope = 1,col="grey",lty = 2)+
  geom_point(aes(fill = short_name, shape = type), size = 2, color = "black")+
  facet_wrap(~cycle,scales = "free",
             labeller = labeller(cycle = function(x) paste0(x, " cycles")))+
  scale_y_log10(limits=c(0.0005,500),breaks = c(0.0005,0.005,0.05,0.5,5,50,500),labels = c("0.0005","0.005","0.05","0.5","5","50","500"))+
  scale_x_log10(limits=c(0.0005,500),breaks = c(0.0005,0.005,0.05,0.5,5,50,500),labels = c("0.0005","0.005","0.05","0.5","5","50","500"))+

  scale_fill_brewer(palette = "Set3")+
  scale_shape_manual(
    values = c("up" = 25, "down" = 24),
    name = "Species"
  ) +
  guides(
    fill = guide_legend(override.aes = list(shape = shape_map[levels(F1_bio_all6$type)])),
    shape = "none"  # 隐藏独立的 shape 图例
  )+
  geom_text(data = cor_results, aes(x = 0.0005, y = 500, label = label), 
            size=2,hjust = 0, vjust = 1,parse = TRUE,inherit.aes = FALSE) +
  labs(
    x="Observed frequency in technical replicate1 (%)",
    y="Observed frequency in technical replicate2 (%)",fill="Alleles")+
  theme_classic()+#用classic模板
  my_theme+
  theme(legend.text = element_text(face = "italic"),
        legend.position = "bottom")


ggsave(F1_bio_pic,filename = paste0(output,"triangle_HLA_F1_2reads_biorepeat_Frequency_minor_v3.pdf"),
       width = 19,
       height = 8,units = "cm")




#Fig4_H:#################################

F1_top3_data <- F1_bio_all4%>%
  mutate(.,short_name=gsub(".*(C\\*\\d+:\\d+:\\d+).*", "\\1", species))%>%
  filter(short_name %in% F1_factor[5:10])%>% 
  select(-count,-bio_cycle_count,-F1)




F1_top3_data_20 <- F1_top3_data%>%
  filter(cycle=="20")%>%
  dplyr::rename(F1_perc_20 = F1_perc)%>%
  select(-cycle)

F1_top3_data_ratio <-F1_top3_data %>%
  left_join(F1_top3_data_20, by = c("biorepeat","short_name","species"))%>%
  mutate(ratio = F1_perc / F1_perc_20) %>%
  select(-F1_perc,-F1_perc_20) 


F1_top3_data_ratio2 <- F1_top3_data_ratio%>%
  group_by(cycle,short_name)%>%
  summarise(SD=sd(ratio),
            Mean=mean(ratio),
            SE = sd(ratio) / sqrt(n()),   
            N = n(),                       
            .groups = "drop")%>%
  mutate(.,Mean_perc=Mean*100,
         SD_perc=SD*100,
         SE_perc = SE * 100)


F1_top3_data_ratio2$short_name <- factor(F1_top3_data_ratio2$short_name,level=unique(F1_factor))



top3_p <- ggplot(F1_top3_data_ratio2, aes(x = fct_rev(short_name), y = Mean_perc, fill = factor(cycle))) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.5) + 
  geom_hline(yintercept = 100, linetype = "dashed", color = "grey", linewidth = 0.4) +
  labs(
    x = NULL, 
    y = "Relative observed frequency (%)", 
    fill = "Cycles" # 图例标题
  ) +
  scale_y_break(c(150, 360),scales = 1)+
  scale_y_continuous(limits = c(0, 800),breaks = c(0,50,100,150,400,600,800))+
  geom_errorbar(aes(ymin = Mean_perc - SE_perc, ymax = Mean_perc + SE_perc), 
                position = position_dodge(width = 0.8),  
                width = 0.1, color = "black")+
  scale_x_discrete(labels = c(
    "C*04:03:01"="C*04:03:01(44.445%)",
    "C*07:04:01"="C*07:04:01(44.445%)",
    "C*15:05:02"="C*15:05:02(5%)",
    "C*15:02:01"="C*15:02:01(5%)",
    "C*06:02:01"="C*06:02:01(0.5%)",
    "C*08:01:01"="C*08:01:01(0.5%)"
  )) +
  theme_classic()+
  my_theme+
  theme(
    axis.line.x.top = element_blank(),  
    axis.ticks.x.top = element_blank(),   
    axis.text.x.top = element_blank(),     
    axis.line.y.right = element_blank(),  
    axis.text.y.right = element_blank(),  
    axis.ticks.y.right = element_blank(),  
    legend.position = "none"
  )


t_test_results <- F1_top3_data_ratio %>%
  group_by(species) %>%
  summarise(
    p_value_25 = t.test(ratio[cycle == "25"], ratio[cycle == "20"])$p.value,  # 25 vs 20
    p_value_30 = t.test(ratio[cycle == "20"], ratio[cycle == "30"])$p.value,   # 20 vs 30
    .groups = "drop"
  ) %>%
  mutate(
    short_name=gsub(".*(C\\*\\d+:\\d+:\\d+).*", "\\1", species),
    label_25 = paste0("italic(P)==",format(signif(p_value_25, digits = 3), scientific = TRUE)), 
    label_30 = paste0("italic(P)==",format(signif(p_value_30, digits = 3), scientific = TRUE))
  )

t_test_results2 <- t_test_results%>%
  mutate(.,
         p_value_25_onetail=p_value_25/2,
         p_value_30_onetail=p_value_30/2,
         label_25_onetail = ifelse(p_value_25_onetail < 0.01,
                                   paste0("italic(P)=='", format(signif(p_value_25_onetail, digits = 3), scientific = TRUE), "'"),
                                   paste0("italic(P)=='", format(round(p_value_25_onetail, 3), nsmall = 3), "'")),
         label_30_onetail = ifelse(p_value_30_onetail < 0.01,
                                   paste0("italic(P)=='", format(signif(p_value_30_onetail, digits = 3), scientific = TRUE), "'"),
                                   paste0("italic(P)=='", format(round(p_value_30_onetail, 3), nsmall = 3), "'")))



species_order <- rev(F1_factor)


for (i in seq_along(t_test_results2$short_name)) {
  species <- t_test_results2$short_name[i]
  species_index <- which(species_order == species) 
  

  top3_p <- top3_p + geom_signif(
    annotations = t_test_results2$label_25_onetail[i],
    y_position = max(F1_top3_data_ratio2$Mean_perc) * 1.3,
    xmin = species_index - 0.26,
    xmax = species_index,
    tip_length = 0.01,
    vjust = 0,
    textsize = 2,
    parse = TRUE
  )
  

  top3_p <- top3_p + geom_signif(
    annotations =t_test_results2$label_30_onetail[i],
    y_position = max(F1_top3_data_ratio2$Mean_perc) * 1.5, 
    xmin = species_index - 0.26 ,
    xmax = species_index + 0.26,
    tip_length = 0.01,
    vjust = 0,
    textsize = 2,
    parse = TRUE  
  )
}



ggsave(top3_p,filename = paste0(output,"se_HLA_F1_top3_ratio_2reads_minor_onetail_v3.pdf"), 
       width = 10,
       height = 5,units = "cm") 



data_Fig4_H <- F1_top3_data_ratio%>%
  mutate(.,ratio=ratio*100)%>%
  select(-species)

colnames(data_Fig4_H) <- c("Repeat","Cycles","Alleles","Relative observed Frequency")

write.csv(data_Fig4_H,file=paste0(output,"data_Fig4_H.csv"),row.names = F)


#Fig4_D##########################
load("/combine/HLA_C_2reads_minor_identity_S8.RData")


S8_2 <- S8%>%
  subset(.,select=c(-Chimeras,-ts_region))%>%
  mutate(species_3 = gsub(".*(C\\*\\d+:\\d+:\\d+).*", "\\1", ref_3),
         species_5 = gsub(".*(C\\*\\d+:\\d+:\\d+).*", "\\1", ref_5))


comp_data_by_group <- S8_2 %>%
  group_by(cycle) %>%  
  summarise(
    samesp = sum(species_3 == species_5),    
    diffsp = sum(species_3 != species_5),    
    all_events = samesp+diffsp,           
    .groups = "drop"                          
  ) %>%
  mutate(
    fre_samesp = samesp / all_events,           
    fre_diffsp = diffsp / all_events            
  )


comp_data_long <- comp_data_by_group %>%
  pivot_longer(
    cols = c(samesp, diffsp),
    names_to = "category",
    values_to = "value"
  ) %>%
  mutate(
    fre = value / all_events 
  )


comp_data_long$category <- factor(comp_data_long$category,levels = c("samesp","diffsp"))


pie_plot <- ggplot(comp_data_long, aes(x = "", y = fre, fill = category)) +
  geom_bar(stat = "identity", width = 1, alpha = 0.9) +
  coord_polar("y", start = 0) +
  geom_text(
    aes(label = paste0(value, "\n", sprintf("%.2f%%", fre * 100))),
    position = position_stack(vjust = 0.5),
    size = 4
  ) +
  scale_fill_manual(
    values = c("#619CFF", "#F8766D"),  
    labels = c("Within the same HLA-C allele", "Between different HLA-C alleles") 
  ) +
  facet_wrap(~ cycle,
             labeller = labeller(cycle = function(x) paste0(x, " cycles"))) + 
  labs(fill = "Category") +
  theme_void() +
  theme(
    strip.text = element_text(size = 15),  
    legend.position = "bottom",            
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  )

ggsave(pie_plot,filename = paste0(output,"HLA_2reads_pie.pdf"), 
       width = 10,
       height = 5) 


data_Fig4_D <- comp_data_by_group%>%
  select(cycle,samesp,diffsp,all_events)

colnames(data_Fig4_D) <- c("Cycles","Number of same HLA-C allele","Number of different HLA-C alleles","All events")

write.csv(data_Fig4_D,file=paste0(output,"data_Fig4_D.csv"),row.names = F)






  
###############################

load("/2reads_minor_identity/1/HLA_C_2reads_minor_identity_S8.RData")

S8_bio1 <- S8%>%mutate(.,biorepeat="1")

load("/2reads_minor_identity/2/HLA_C_2reads_minor_identity_S8.RData")

S8_bio2 <- S8%>%mutate(.,biorepeat="2")

S8_bio <- rbind(S8_bio1,S8_bio2)%>%
  subset(.,select=c(-Chimeras,-ts_region))%>%
  mutate(species_3 = gsub(".*(C\\*\\d+:\\d+:\\d+).*", "\\1", ref_3),
         species_5 = gsub(".*(C\\*\\d+:\\d+:\\d+).*", "\\1", ref_5))  
  
#Fig4_E###########################

F2_data <- S8_bio%>%
  filter(.,ref_3=="K32C1_C*04:03:01:01")
  
  
F2_by_group <- F2_data %>%
  group_by(cycle,biorepeat,species_3,species_5) %>% 
  summarise(
    count = n(),
    .groups = "drop"                              
  ) 


load("/HLA/identity_heatmap/ref_identity_score_gap.RData")
identity_scores_2 <- identity_scores_3


lis_1_identity <- identity_scores_2%>%
  subset(.,select=c("K32C1_C*04:03:01:01"))


lis_1_identity <- lis_1_identity %>% 
  mutate(row_id = row.names(.),
         short_name = gsub(".*(C\\*\\d+:\\d+:\\d+).*", "\\1", row_id))%>%
  dplyr::rename(score = "K32C1_C*04:03:01:01")%>%
  subset(.,select=c(-row_id))

rownames(lis_1_identity) <- NULL  #



F2_merged <- F2_by_group %>%
  group_by(cycle,biorepeat) %>% 
  dplyr::group_modify(~ {   
    merge(.x, lis_1_identity, 
          by.x = "species_5", 
          by.y = "short_name", 
          all = TRUE) 
  }) %>%
  ungroup()%>%  
  subset(.,select=c(-species_3))%>%
  mutate(
    count = ifelse(is.na(count), 0, count)
  )

F2_merged2 <- F2_merged%>%
  group_by(cycle,biorepeat) %>% 
  mutate(cycle_biorepeat_count = sum(count)) %>% 
  ungroup()%>% 
  mutate(.,F2=count/cycle_biorepeat_count,
         F2_perc=F2*100)%>%
  group_by(cycle,species_5,score) %>%
  summarise(
    SD = sd(F2_perc),               
    Mean = mean(F2_perc),            
    SE = sd(F2_perc) / sqrt(n()),   
    N = n(),                      
    .groups = "drop"             
  )


F2_merged3 <- F2_merged2%>%
  mutate(Label = paste(species_5, cycle, sep = ";")) %>%
  arrange(cycle,desc(score),desc(Mean))

F2_merged3$Label <- factor(F2_merged3$Label, levels = unique(F2_merged3$Label))


custom_labeller <- function(x) {
  sapply(x, function(label) {
    unlist(strsplit(label, ";"))[1]
  })
}


F2_plot_pic <-  ggplot(data=F2_merged3,aes(x=Label,y=Mean))+
  geom_blank() +
  geom_bar(stat = "identity",width = 0.5,fill="#FB8072",data = ~ subset(.x, Mean != 0))+  
  facet_wrap(~cycle,scales = "free",
             labeller = labeller(cycle = function(x) paste0(x, " cycles")))+
  geom_errorbar(aes(ymin = Mean-SE, ymax = Mean + SE), width = 0.1, color = "black") +
  geom_smooth(aes(group = 1), method = "loess",  span = 0.5,method.args = list(degree = 2),se = FALSE, color = "grey", linewidth = 0.5) + 
  scale_x_discrete(labels = custom_labeller)+ 
  scale_y_continuous(limits = c(0, 90), breaks = c(0, 20, 40, 60, 80))+
  labs(
    x="Sequence similarity",
    y="Percentage of template switching events (%)")+
  theme_classic()+
  my_theme+
  theme(axis.text.x = element_text(face = "italic",angle = 45,hjust=1,size = 5))

ggsave(F2_plot_pic,filename = paste0(output_se_F23,"se_HLA_F2_2reads_reads_minor_v3.pdf"), 
       width = 10,
       height = 6,units = "cm")


data_Fig4_E <- F2_merged%>%
  group_by(cycle,biorepeat) %>% 
  mutate(cycle_biorepeat_count = sum(count)) %>%  
  ungroup()%>% 
  mutate(.,F2=count/cycle_biorepeat_count,
         F2_perc=F2*100)%>%
  arrange(cycle,biorepeat,desc(score),desc(F2_perc))%>%
  select(cycle,biorepeat,species_5,count,cycle_biorepeat_count,F2_perc,score)

colnames(data_Fig4_E) <- c("Cycles","Repeat","Alleles","Number of events","Number of all events","Percentage","Sequence similarity")

write.csv(data_Fig4_E,file=paste0(output_se_F23,"data_Fig4_E.csv"),row.names = F)




#Fig4_F#####################


F3_data <- S8_bio%>%
  group_by(cycle,biorepeat,species_3) %>% 
  summarise(
    F3_count = n(),
    .groups = "drop"                           
  ) 



F3_merged <- F3_data %>%
  group_by(cycle,biorepeat) %>% 
  dplyr::group_modify(~ {   
    merge(.x, standard_freq_data%>%select(-species), 
          by.x = "species_3", 
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
  group_by(cycle,species_3,frequency2) %>%
  summarise(
    SD = sd(F3_perc),                
    Mean = mean(F3_perc),           
    SE = sd(F3_perc) / sqrt(n()),  
    N = n(),                     
    .groups = "drop"              
  )%>%
  arrange(cycle,desc(frequency2),desc(Mean))



F3_merged2$species_3 <- factor(F3_merged2$species_3,levels = rev(unique(F1_factor)))


#改为means ± s.e.!!


F3_plot_pic2 <- ggplot(data=F3_merged2,aes(x=species_3,y=Mean))+
  geom_blank() +
  geom_bar(stat = "identity",width = 0.5,fill="#FB8072",data = ~ subset(.x, Mean != 0))+
  facet_wrap(~cycle,scales = "free",
             labeller = labeller(cycle = function(x) paste0(x, " cycles")))+
  geom_errorbar(aes(ymin = Mean-SE, ymax = Mean + SE), width = 0.1, color = "black") +
  geom_smooth(aes(group = 1), method = "loess", span = 0.55,method.args = list(degree = 1),se = FALSE, color = "grey", linewidth = 0.5) +
  scale_x_discrete(labels = custom_labeller)+
  scale_y_continuous(limits = c(0, 68), breaks = c(0, 20, 40, 60))+
  labs(
    x="Expected frequency",
    y="Percentage of template switching events (%)")+
  theme_classic()+#用classic模板
  my_theme+
  theme(axis.text.x = element_text(face = "italic",angle = 45,hjust=1,size = 5)) 


ggsave(F3_plot_pic2,filename = paste0(output_se_F23,"se_HLA_F3_2reads_minor_v3_standard.pdf"),
       width = 10,
       height = 6,units = "cm") 



data_Fig4_F <- F3_merged%>%
  group_by(cycle,biorepeat) %>%  
  mutate(cycle_biorepeat_count = sum(F3_count)) %>%  
  ungroup()%>%  
  mutate(.,F3=F3_count/cycle_biorepeat_count,
         F3_perc=F3*100)%>%
  arrange(cycle,biorepeat,desc(frequency2),desc(F3_perc))%>%
  select(cycle,biorepeat,species_3,F3_count,cycle_biorepeat_count,F3_perc,frequency2)

colnames(data_Fig4_F) <- c("Cycles","Repeat","Alleles","Number of events","Number of all events","Percentage","Expected frequency")

write.csv(data_Fig4_F,file=paste0(output_se_F23,"data_Fig4_F.csv"),row.names = F)



#Fig4_G#####################


load("/ts_plot/HLA_C_2reads_minor_identity/plot2.RData")


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
    color = "tomato"
  )  +  
  geom_point(size=0.5) +  
  facet_wrap(~cycle, nrow = 1, scales = "free_y") +  
  labs(x = "Position in HLA-C sequence (bp)", y = "Template switching score") +  
  theme_classic() +  
  scale_y_continuous(limits = c(0, 22), breaks = c(0, 5,10,15,20)) +  
  my_theme

ggsave(combined_plot,filename = paste0(output,"HLAC_2reads_ts_score.pdf"), 
       width = 10,
       height = 5.5,units = "cm")



data_fig4_G1 <- combined_data%>%mutate(.,cycle=str_extract(cycle, "\\d+"))%>%
  select(cycle,sequence_index,Count)
colnames(data_fig4_G1) <- c("Cycles","Position", "Template Switch Score" )

write.csv(data_fig4_G1,file=paste0(output,"data_fig4_G1.csv"),row.names = F)













