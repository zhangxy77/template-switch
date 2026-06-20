library(plyr)
library(dplyr)
library(ggbreak)
library(readxl)
library(ggplot2)

my_theme <- theme(legend.position = "bottom",
                  plot.title=element_text(hjust=0.5,size=8),
                  axis.text = element_text(size = 6),
                  strip.text = element_text(size = 8),
                  axis.title = element_text(size = 8),
                  legend.title = element_text(size = 8),
                  legend.text = element_text(size = 6)) 


#S8c#############


S8c_geno <- read.table("/Umi3_geno_fig1_reads_2/Ura_diff_length_All_OtoM_ratio_3umi_2.txt",
                       header = T)%>%
  filter(.,Template=="500" & Length=="1000" & Enzyme=="Q" & Cycles == "20")%>%
  mutate(.,Group="geno")%>%
  select(c(Sample,Template,Cycles,Biorepeat,Ratio,Group))


S8c_umi5 <- read.table("/Umi3_umi5_fig1_reads_2/16S_50ng_All_OtoM_ratio_3umi_5umi_2.txt",
                       header = T)%>%
  filter(.,Cycles == "20")%>%
  mutate(.,Group="umi5")%>%
  select(c(Sample,Template,Cycles,Biorepeat,Ratio,Group))

S8c_data <- rbind(S8c_geno,S8c_umi5)


S8c_data2 <- S8c_data%>%
  dplyr::group_by(Cycles,Group)%>%
  summarise(ratio_mean=mean(Ratio),
            ratio_sd=sd(Ratio),
            sample_size = n(),
            ratio_se=ratio_sd/sqrt(sample_size),
            .groups = "drop")

S8c_data2$Cycles <- factor(S8c_data2$Cycles)


S8c_t_test_results <- S8c_data %>%
  group_by(Cycles) %>%
  summarise(
    p_value = t.test(Ratio ~ Group)$p.value,
    .groups = "drop"
  ) %>%
  mutate(
    label_two_sided=paste("P = ",format(p_value, digits = 3,scientific = TRUE)),
    label_one_sided=paste("P = ",format(p_value/2, digits = 3,scientific = TRUE))
  )


write.csv(S8c_t_test_results,file=paste0(outpath,"figS8c_p_value_t_test_v3.csv"))


S8c_pic <- ggplot(S8c_data2,aes(x=Group,y=ratio_mean,fill=Group))+
  geom_col(
    position = position_dodge(width = 0.8), 
    width = 0.5)+ 
  geom_errorbar(aes(ymin = ratio_mean-ratio_se, ymax = ratio_mean+ratio_se,color=Group), width = 0.2,position = position_dodge(width = 0.8))+
  labs(x = NULL) +
  scale_x_discrete(labels=c("DMS","MMC"))+
  scale_y_continuous(
    breaks =c(0,1,2,3),
    name = "MtoO ratio1", 
    limits = c(0, max(S8c_data2$ratio_mean + S8c_data2$ratio_se) * 1.2),
    expand = expansion(mult = c(0.05, 0.1)), 
    sec.axis = sec_axis(
      ~ .,  
      breaks =c(0,1,2,3),
      name = "MtoO ratio2"
    )
  ) +
  theme_classic()+
  my_theme


ggsave(S8c_pic,filename = paste0(outpath,"S8c_pic_v3.pdf"), 
       width = 7.15,
       height = 7.5,units = "cm")





#S8d#############

S8d_geno <- read.table("/Umi3_geno_fig1_reads_2/Ura_diff_length_minor_pro_3umi_2.txt",
                       header = T)%>%
  filter(.,Template=="500" & Length=="1000" & Enzyme=="Q" & Cycles == "20")%>%
  mutate(.,Group="geno")%>%
  select(c(Sample,Template,Cycles,Biorepeat,Minor_proportion,Group))


S8d_umi5 <- read.table("/Umi3_umi5_fig1_reads_2/16S_50ng_minor_pro_3umi_5umi_2.txt",
                       header = T)%>%
  filter(.,Cycles == "20")%>%
  mutate(.,Group="umi5")%>%
  select(c(Sample,Template,Cycles,Biorepeat,Minor_proportion,Group))


S8d_data <- rbind(S8d_geno,S8d_umi5)



S8d_data2 <- S8d_data%>%
  dplyr::group_by(Cycles,Group)%>%
  summarise(ratio_mean=mean(Minor_proportion),
            ratio_sd=sd(Minor_proportion),
            sample_size = n(),
            ratio_se=ratio_sd/sqrt(sample_size),
            .groups = "drop")

S8d_data2$Cycles <- factor(S8d_data2$Cycles)



S8d_t_test_results <- S8d_data %>%
  group_by(Cycles) %>%
  summarise(
    p_value = t.test(Minor_proportion ~ Group)$p.value,
    .groups = "drop"
  ) %>%
  mutate(
    label_two_sided=paste("P = ",format(p_value, digits = 3,scientific = TRUE)),
    label_one_sided=paste("P = ",format(p_value/2, digits = 3,scientific = TRUE))
  )


write.csv(S8d_t_test_results,file=paste0(outpath,"figS8d_p_value_t_test_v3.csv"))



S8d_pic <- ggplot(S8d_data2,aes(x=Group,y=ratio_mean,fill=Group))+
  geom_col(
    position = position_dodge(width = 0.8),  
    width = 0.5)+ 
  geom_errorbar(aes(ymin = ratio_mean-ratio_se, ymax = ratio_mean+ratio_se,color=Group), width = 0.2,position = position_dodge(width = 0.8))+
  labs(x = NULL) +
  scale_x_discrete(labels=c("DMS","MMC"))+
  scale_y_continuous(
    name = "Minor ratio1",  
    limits = c(0, max(S8d_data2$ratio_mean + S8d_data2$ratio_se) * 1.2),
    expand = expansion(mult = c(0.05, 0.1)),  
    sec.axis = sec_axis(
      ~ .,  
      name = "Minor ratio2" 
    )
  ) +
  theme_classic()+
  my_theme

ggsave(S8d_pic,filename = paste0(outpath,"S8d_pic_v3.pdf"), 
       width = 7.5,
       height = 7.5,units = "cm")









