library(plyr)
library(dplyr)
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


#S7a##################


S7a_data <- read.table("/Umi3_umi5_fig1_reads_2/Com_repeat_16S_50ng_3umi_umi5_types_3umi_2.txt",header = T)%>%
  filter(.,Geno_type_num <=5)

S7a_data$Cycles <- as.factor(S7a_data$Cycles)


S7a_f <- ggplot(S7a_data,aes(x=Cycles,y=Geno_type_num,fill=Cycles))+
  geom_violin(scale = "width",trim=FALSE)+
  labs(
    x=NULL,y="The diversity of 5' UMIs in each 3' UMI")+ 
  theme_classic()+
  my_theme

ggsave(S7a_f,filename = paste0(outpath,"S7a_f.pdf"),
       width = 8.5,
       height = 8.5,units = "cm")


#S7c:#################

S7c_data <- read.table("/Umi3_umi5_fig1_reads_2/Com_repeat_16S_50ng_3umi_umi5_major_pro_3umi_2.txt",header = T)

S7c_data$Cycles <- as.factor(S7c_data$Cycles)


S7c_f <- ggplot(S7c_data,aes(x=Cycles,y=Major_proportion,fill=Cycles))+
  geom_violin(scale = "width",trim=FALSE)+
  scale_y_continuous(breaks = c(0,0.25,0.5,0.75,1.0,1.25))+
  labs(
    x=NULL,
    y="The frequency of HiFi reads of major 5' UMI in each 3' UMI")+ 
  theme_classic()+
  my_theme

ggsave(S7c_f,filename = paste0(outpath,"S7c_f.pdf"), 
       width = 8.5,
       height = 8.5,units = "cm")


#S7b:#################

S7b_data <- read.table("/Umi3_umi5_fig1_reads_2/16S_50ng_All_OtoM_ratio_3umi_5umi_2.txt",header = T)


S7b_data2 <- S7b_data%>%
  dplyr::group_by(Cycles)%>%
  summarise(ratio_mean=mean(Ratio),
            ratio_sd=sd(Ratio),
            sample_size = n(),
            ratio_se=ratio_sd/sqrt(sample_size))

S7b_data2$Cycles <- factor(S7b_data2$Cycles)


S7b_f <- ggplot(S7b_data2,aes(x=Cycles,y=ratio_mean,fill=Cycles))+
  geom_bar(stat = "identity",width=0.5)+
  geom_errorbar(aes(ymin = ratio_mean-ratio_se, ymax = ratio_mean+ratio_se,color=Cycles), width = 0.3) +
  labs(x=NULL,y="MtoO Ratio")+ 
  theme_classic()+
  my_theme


ggsave(S7b_f,filename = paste0(outpath,"S7b_f.pdf"), 
       width = 8.5,
       height = 8.5,units = "cm")



#S7d:#################

S7d_data <- read.table("/Umi3_umi5_fig1_reads_2/16S_50ng_minor_pro_3umi_5umi_2.txt",header = T)


S7d_data2 <- S7d_data%>%
  dplyr::group_by(Cycles)%>%
  summarise(ratio_mean=mean(Minor_proportion),
            ratio_sd=sd(Minor_proportion),
            sample_size = n(),
            ratio_se=ratio_sd/sqrt(sample_size))


S7d_data2$Cycles <- factor(S7d_data2$Cycles)

S7d_f <- ggplot(S7d_data2,aes(x=Cycles,y=ratio_mean,fill=Cycles))+
  geom_bar(stat = "identity",width=0.5)+
  geom_errorbar(aes(ymin = ratio_mean-ratio_se, ymax = ratio_mean+ratio_se,color=Cycles), width = 0.3) +
  labs(x=NULL,y="Minor Perc")+ 
  theme_classic()+
  my_theme

ggsave(S7d_f,filename = paste0(outpath,"S7d_f.pdf"), 
       width = 8.5,
       height = 8.5,units = "cm")






