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



#S11b:##################


S11b_data <- read.table("Com_repeat_HLAC_1ng_3umi_umi5_types_3umi_1.txt",header = T)%>%
  filter(.,Geno_type_num <=5)

S11b_data$Cycles <- as.factor(S11b_data$Cycles)

S11b_f <- ggplot(S11b_data,aes(x=Cycles,y=Geno_type_num,fill=Cycles))+
  geom_violin(scale = "width",trim=FALSE)+
  labs(
    x=NULL,y="The diversity of 5' UMIs in each 3' UMI")+ 
  theme_classic()+
  my_theme

ggsave(S11b_f,filename = paste0(outpath,"S11b_f.pdf"),
       width = 3,
       height = 5,units = "cm")




#S11c:#################

S11c_data <- read.table("HLAC_1ng_All_OtoM_ratio_3umi_5umi_1.txt",header = T)


S11c_data2 <- S11c_data%>%
  dplyr::group_by(Cycles)%>%
  summarise(ratio_mean=mean(Ratio),
            ratio_sd=sd(Ratio),
            sample_size = n(),
            ratio_se=ratio_sd/sqrt(sample_size))

S11c_data2$Cycles <- factor(S11c_data2$Cycles)


S11c_f <- ggplot(S11c_data2,aes(x=Cycles,y=ratio_mean,fill=Cycles))+
  geom_bar(stat = "identity",width=0.5)+
  geom_errorbar(aes(ymin = ratio_mean-ratio_se, ymax = ratio_mean+ratio_se,color=Cycles), width = 0.3) +
  labs(x=NULL,y="MtoO Ratio")+ 
  theme_classic()+
  my_theme


ggsave(S11c_f,filename = paste0(outpath,"S11c_f.pdf"), 
       width = 3,
       height = 5,units = "cm")





BC <- S11b_f+S11c_f

ggsave(BC,filename = paste0(outpath,"S11bc_f.pdf"), 
       width = 9,
       height = 7.5,units = "cm")











