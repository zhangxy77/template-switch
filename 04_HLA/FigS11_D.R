#250702
#新HLA（1ng）数据、2reads
#S11D、E：利用muscle比对后确定每个umi3的major umi5和minor umi5来统计频率



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



#S10_D##################

cy20 <- read.table("HLAC_1ng_20.txt")%>%
  mutate(cycles=20)
cy25 <- read.table("HLAC_1ng_25.txt")%>%
  mutate(cycles=25)
cy30 <- read.table("HLAC_1ng_30.txt")%>%
  mutate(cycles=30)

combine_all <- rbind(cy20,cy25,cy30)


#改列名
combine_all <- combine_all %>%
  dplyr::rename(
    "zmw_hole" = "V1",
    "index" = "V2",
    "geno" = "V3",
    "umi5" = "V4",
    "umi3" = "V5"
  )%>%
  select(c(umi3,umi5,cycles))



all_reads_result <- combine_all %>%

  group_by(umi3, umi5,cycles) %>%
  summarise(
    umi5_count = n(),  
    .groups = "drop"
  ) %>%

  group_by(umi3, cycles) %>%
  mutate(
    umi3_count = sum(umi5_count)
  ) %>%
  ungroup() 



umi3_reads_filter <- all_reads_result%>%filter(.,umi3_count>=2)



load("umi_ts_20_F1.RData")
load("umi_ts_25_F1.RData")
load("umi_ts_30_F1.RData")


umi_ts_20 <- umi_ts_20%>%mutate(cycles=20)
umi_ts_25 <- umi_ts_25%>%mutate(cycles=25)
umi_ts_30 <- umi_ts_30%>%mutate(cycles=30)


umi_ts_all_data <- rbind(umi_ts_20,umi_ts_25,umi_ts_30)%>%
  select(c(umi3,umi5,type,cycles))




merged_df <- umi_ts_all_data %>%
  left_join(
    umi3_reads_filter,
    by = c("umi3", "umi5", "cycles") 
  )%>%
  mutate(.,major_freq=umi5_count/umi3_count)


#检查是否有na值
na_rows_all <- merged_df %>%
  filter(if_any(everything(), is.na))
rm(na_rows_all)


major_data <- merged_df%>%filter(.,type=="major")

major_data$cycles <- as.factor(major_data$cycles)

S11d_f <- ggplot(major_data,aes(x=cycles,y=major_freq,fill=cycles))+
  geom_violin(scale = "width",trim=FALSE)+
  scale_y_continuous(breaks = c(0,0.25,0.5,0.75,1.0,1.25))+
  labs(
    x=NULL,
    y="The frequency of HiFi reads of major 5' UMI in each 3' UMI")+ 
  theme_classic()+
  my_theme



ggsave(S11d_f,filename = paste0(outpath,"S11d_f.pdf"), 
       width = 8.5,
       height = 8.5,units = "cm")






