library(patchwork)
library(plyr)
library(dplyr)
library(reshape2, lib.loc = "/usr/local/lib64/R/library")
library(ggplot2)
library(tidyr)
library(purrr)


my_theme <- theme(legend.position = "bottom",
                  plot.title=element_text(hjust=0.5,size=5),
                  axis.text = element_text(size = 4),
                  strip.text = element_text(size = 5),
                  axis.title = element_text(size = 5),
                  legend.title = element_text(size = 5),
                  legend.text = element_text(size = 4)) 




simu_result_minor <- function(data_path){
  
  files <- list.files(data_path)
  
  all_result <- lapply(files,function(u){
    
    raw <- read.table(paste(data_path,u,sep=""),header = T)%>%
      dplyr::rename(cycle_20=V1,cycle_30=V2,cycle_40=V3)
    

    rownames(raw) <- c("OtoM_OtoO_Ratio", "OtoM", "OtoO","Minor_Ratio","Major_reads","Minor_reads")
    
    raw2 <- t(raw)%>%
      as.data.frame()%>%
      mutate(Cycles = as.numeric(gsub("cycle_", "", rownames(.))))
    

    rownames(raw2) <- NULL
    
    splitted <- unlist(strsplit(u, split = "_"))
    
    data <- data.frame(
      Len = rep(as.numeric(gsub("bp", "", splitted[1])),3),
      Template_molecules = rep(as.numeric(splitted[2]),3),
      wildtype_ratio = rep(as.numeric(gsub("wt", "", splitted[3])),3),
      error_rate = rep(as.numeric(splitted[4]),3),
      throughput_reads=rep(as.numeric(splitted[5]),3),
      ID=rep(as.numeric(gsub(".txt","",splitted[7])),3)
    )
    
    raw3 <- cbind(data,raw2)
    
    return(raw3)
    
  })%>%rbind.fill()
  
  return(all_result)
}

#Fig1_H#########################


sample_data <- simu_result_minor("sample_result/")%>%
  filter(.,Cycles=="20")%>%
  mutate(.,error_rate_per_site=error_rate,
         log_error_rate_per_site=log10(error_rate_per_site))

sample_data2 <- sample_data%>%
  subset(.,select = c(-OtoM,-OtoO,-Major_reads,-Minor_reads))%>%
  dplyr::group_by(Len,Template_molecules,wildtype_ratio,error_rate,throughput_reads,Cycles,error_rate_per_site,log_error_rate_per_site)%>%
  summarise(sample_MtoO_Ratio_mean=mean(OtoM_OtoO_Ratio),
            sample_Minor_Ratio_mean=mean(Minor_Ratio),
            .groups = 'drop')%>%
  mutate(.,log_MtoO_mean=log10(sample_MtoO_Ratio_mean),
         log_Minor_mean=log10(sample_Minor_Ratio_mean))




Minor_data <- read.table("/Umi3_geno_fig1_reads_2/Ura_diff_length_minor_pro_3umi_2.txt",header = T)%>%
  filter(.,Length=="1000" & Cycles=="20" &Template=="50" &Enzyme=="Q")


Minor_data2 <- Minor_data %>%
  summarise(sample_Minor_Ratio_mean = mean(Minor_proportion, na.rm = TRUE))%>%
  mutate(log_Minor_mean=log10(sample_Minor_Ratio_mean))


model_lm <- lm(sample_Minor_Ratio_mean~error_rate_per_site, data = sample_data2)
summary(model_lm)
sample_data2$prediction <- predict(model_lm,newdata =sample_data2)

f4 <- function(x, a, model) {
  predict(model, newdata = data.frame(error_rate_per_site = x)) - a
}

x_pred <- uniroot(f4, interval = c(-1,1), a = Minor_data2$sample_Minor_Ratio_mean, model = model_lm)$root

pic <- ggplot(sample_data,aes(x=error_rate_per_site,y=Minor_Ratio))+
  geom_point(size=0.5)+
  geom_hline(yintercept = Minor_data2$sample_Minor_Ratio_mean, color = "#F8766D",linetype = "dashed")
  geom_line(data=sample_data2,aes(x=error_rate_per_site,y=prediction),color="#00BFC4")+
  scale_x_continuous(
    breaks = seq(from = 3, to = 10, by = 1)*10^(-6),
  ) +
  labs(x="TSR",y="Minor ratio")+
  theme_classic()+
  my_theme


ggsave(pic,filename = paste0(outpath,"wt0.5_8TSR_data_Fig1_fitted_curve2.pdf"), 
       width = 9,
       height = 4.5,units = "cm")



data_fig1_H <- sample_data%>%
  select(c(Template_molecules,wildtype_ratio,Cycles,throughput_reads,error_rate_per_site,Minor_Ratio))


colnames(data_fig1_H) <- c("Templates","Wild-type proportion", "Cycles",
                           "Number of sampling sequences","TSR","Minor genotypes frequency" )

write.csv(data_fig1_H,file=paste0(outpath,"data_fig1_H.csv"),row.names = F)





#Fig1_I###################

reads_sample_data <- simu_result_minor("sample_result4/")%>%
  filter(.,Cycles=="20")%>%
  mutate(.,log_throughput_reads=log10(throughput_reads),
         error_rate_per_site=error_rate,
         log_error_rate_per_site=log10(error_rate_per_site))


reads_sample_data2 <- reads_sample_data%>%
  subset(.,select = c(-OtoM,-OtoO,-Major_reads,-Minor_reads))%>%
  dplyr::group_by(Len,Template_molecules,wildtype_ratio,error_rate,throughput_reads,log_throughput_reads,Cycles,error_rate_per_site,log_error_rate_per_site)%>%
  summarise(sample_OtoM_Ratio_mean=mean(OtoM_OtoO_Ratio),
            sample_Minor_Ratio_mean=mean(Minor_Ratio),
            .groups = 'drop')



diff_reads_minor_pic <- ggplot(reads_sample_data,aes(x=factor(throughput_reads),y=Minor_Ratio))+
  geom_point(size=0.5)+
  geom_smooth(aes(group = 1), method = "lm", span = 0.7,se = FALSE, color = "#00BFC4", linewidth = 0.5) +
  labs(x="the number of sampling sequences",y="Minor ratio")+
  scale_y_continuous(limits = c(0.01, 0.05))+
  theme_classic()+
  my_theme


ggsave(diff_reads_minor_pic,filename = paste0(outpath,"geom_smooth_wt0.5_diff_reads_minor_pic2_Fig1_fitted_curve.pdf"), 
       width = 9,
       height = 4.5,units = "cm")



data_fig1_I2 <- reads_sample_data%>%
  select(c(Template_molecules,wildtype_ratio,Cycles,error_rate_per_site,throughput_reads,Minor_Ratio))


colnames(data_fig1_I2) <- c("Templates","Wild-type proportion", "Cycles","TSR",
                           "The number of sampling sequences","Minor genotypes frequency" )

write.csv(data_fig1_I2,file=paste0(outpath,"data_fig1_I2.csv"),row.names = F)


