
library(dplyr)
library(ggbreak)
library(readxl)
library(ggplot2)


my_theme <- theme(legend.position = "bottom",
                  plot.title=element_text(hjust=0.5,size=20),
                  axis.text = element_text(size = 15),
                  strip.text = element_text(size = 18),
                  axis.title = element_text(size = 18),
                  legend.title = element_text(size = 16),
                  legend.text = element_text(size = 14))



all_data <- read.table("final_table_umi_2.txt",header = T)
         


#Fig1_C:#########################

c_data <- all_data%>%
  dplyr::select(Template,Length2,Enzyme,Cycles,Biorepeat,UMI_num_filter,OnetoOne_UMI,OnetoMore_UMI)%>%
  mutate(.,c_ratio=OnetoMore_UMI/OnetoOne_UMI)


c_data2 <- c_data%>%
  dplyr::group_by(Template,Length2,Enzyme,Cycles)%>%
  summarise(c_ratio_mean=mean(c_ratio),
            c_ratio_sd=sd(c_ratio),
            sample_size = n(),
            c_ratio_se=c_ratio_sd/sqrt(sample_size),
            .groups = "drop")  

c_data2$Template <- factor(c_data2$Template)

#Q5

c_Qdata <- c_data2%>%filter(.,Enzyme=="Q")%>%arrange(Template,Length2,Cycles)

c_Qdata_FC <- c_Qdata%>%
  mutate(.,FC=(c_ratio_mean)/(c_Qdata[1,5]%>%as.numeric()))


c_Q_his <- ggplot(data=c_Qdata_FC,aes(x=Template,y=c_ratio_mean,fill=Template))+
  geom_bar(stat = "identity")+
  geom_errorbar(aes(ymin = c_ratio_mean-c_ratio_se, ymax = c_ratio_mean+c_ratio_se,color=Template), width = 0.3) +
  labs(x=NULL,y="MtoO Ratio")+ 
  theme_classic()+
  facet_wrap(~Length2+Cycles,scales = "fixed",ncol = 9, nrow = 1,
             labeller = labeller(Length2 = function(x) paste0(x, " bp"),
                                 Cycles = function(x) paste0(x, " cycles")))+
  scale_fill_discrete(labels = c("50ng/100μl","500ng/100μl"))+
  scale_color_discrete(labels = c("50ng/100μl","500ng/100μl"))+
  scale_y_break(c(3, 6),scales = 1)+
  scale_y_break(c(0.4, 1),scales = 1)+
  scale_y_continuous(sec.axis=sec_axis(~./(c_Qdata[1,5]%>%as.numeric()),name="Fold Change",breaks = c(1,50,100,300,500,700,1000,1600,1900,2200)))+#双y轴
  geom_hline(yintercept = c_Qdata[1,5]%>%as.numeric(), col="#404040",lty = 2)+
  labs(title="fig1_c_Q5_se")+ 
  my_theme+
  theme(legend.position = "bottom")

ggsave(c_Q_his,filename = paste0(outpath,"fig1_c_Q5_se.pdf"),
       device = cairo_pdf,  
       width = 14,
       height = 6)



data_fig1_C <- c_data%>%
  filter(.,Enzyme=="Q")%>%
  mutate(.,Enzyme="Q5")%>%
  select(-UMI_num_filter)


colnames(data_fig1_C) <- c("Templates","Length", "Enzyme","Cycles","Repeat",
                           "The number of barcodes with single genotype",
                           "The number of barcodes with multiple genotypes",
                           "Ratio")

write.csv(data_fig1_C,file=paste0(outpath,"data_fig1_C.csv"),row.names = F)




#Phusion

c_Pdata <- c_data2%>%filter(.,Enzyme=="P")%>%arrange(Template,Length2,Cycles)

c_Pdata_FC <- c_Pdata%>%
  mutate(.,FC=(c_ratio_mean)/(c_Pdata[1,5]%>%as.numeric()))

c_P_his<- ggplot(data=c_Pdata_FC,aes(x=Template,y=c_ratio_mean,fill=Template))+
  geom_bar(stat = "identity")+
  geom_errorbar(aes(ymin = c_ratio_mean-c_ratio_se, ymax = c_ratio_mean+c_ratio_se,color=Template), width = 0.3) +
  labs(x=NULL,y="MtoO Ratio")+ 
  theme_classic()+
  facet_wrap(~Length2+Cycles,scales = "fixed",ncol = 9, nrow = 1,
             labeller = labeller(Length2 = function(x) paste0(x, " bp"),
                                 Cycles = function(x) paste0(x, " cycles")))+
  scale_fill_discrete(labels = c("50ng/100μl","500ng/100μl"))+
  scale_color_discrete(labels = c("50ng/100μl","500ng/100μl"))+
  scale_y_break(c(2.5, 6),scales = 1)+
  scale_y_break(c(0.3, 1),scales = 1)+
  scale_y_continuous(sec.axis=sec_axis(~./(c_Pdata[1,5]%>%as.numeric()),name="Fold Change",breaks = c(1,12,24,120,180,240,600,700,800)))+#双y轴
  geom_hline(yintercept = c_Pdata[1,5]%>%as.numeric(), col="#404040",lty = 2,linewidth=0.5)+
  my_theme+
  theme(legend.position = "bottom")


ggsave(c_P_his,filename = paste0(outpath,"figS2b_Phusion_se.pdf"), 
       device = cairo_pdf,
       width = 14,
       height = 6,units = "cm")





#Fig1_E:######################

e_data <- all_data%>%
  dplyr::select(Template,Length2,Enzyme,Cycles,Biorepeat,UMI_reads_sum_with_major,major_reads,minor_reads)%>%
  mutate(.,e_ratio=minor_reads/UMI_reads_sum_with_major)



e_data2 <- e_data%>%
  dplyr::group_by(Template,Length2,Enzyme,Cycles)%>%
  summarise(e_ratio_mean=mean(e_ratio),
            e_ratio_sd=sd(e_ratio),
            sample_size = n(),
            e_ratio_se=e_ratio_sd/sqrt(sample_size))

e_data2$Template <- factor(e_data2$Template)

#Q5

e_Qdata <- e_data2%>%filter(.,Enzyme=="Q")%>%arrange(Template,Length2,Cycles)

e_Qdata_FC <- e_Qdata%>%
  mutate(.,FC=(e_ratio_mean)/(e_Qdata[1,5]%>%as.numeric()))

e_Q_his <- ggplot(data=e_Qdata_FC,aes(x=Template,y=e_ratio_mean,fill=Template))+
  geom_bar(stat = "identity")+
  geom_errorbar(aes(ymin = e_ratio_mean-e_ratio_se, ymax = e_ratio_mean+e_ratio_se,color=Template), width = 0.3) +
  labs(x=NULL,y="minor/all")+ 
  theme_classic()+
  facet_wrap(~Length2+Cycles,scales = "fixed",ncol = 9, nrow = 1,
             labeller = labeller(Length2 = function(x) paste0(x, " bp"),
                                 Cycles = function(x) paste0(x, " cycles")))+
  scale_fill_discrete(labels = c("50ng/100μl","500ng/100μl"))+
  scale_color_discrete(labels = c("50ng/100μl","500ng/100μl"))+
  scale_y_break(c(0.07, 0.14),scales = 1)+
  scale_y_continuous(sec.axis=sec_axis(~./(e_Qdata[1,5]%>%as.numeric()),name="Fold Change",breaks = c(1,30,60,100,200,300,400)))+#双y轴
  geom_hline(yintercept = e_Qdata[1,5]%>%as.numeric(), col="#404040",lty = 2)+
  labs(title="fig1_e_Q5_se")+
  my_theme+
  theme(
    legend.position = "bottom"
  )


ggsave(e_Q_his,filename = paste0(outpath,"fig1_e_Q5_se.pdf"), 
       device = cairo_pdf, 
       width = 14,
       height = 6)



data_fig1_E <- e_data%>%
  filter(.,Enzyme=="Q")%>%
  mutate(.,Enzyme="Q5")%>%
  select(-major_reads)


colnames(data_fig1_E) <- c("Templates","Length", "Enzyme","Cycles","Repeat",
                           "The number of HiFi reads in all genotypes",
                           "The number of HiFi reads in all minor genotypes",
                           "Ratio")

write.csv(data_fig1_E,file=paste0(outpath,"data_fig1_E.csv"),row.names = F)


#Phusion

e_Pdata <- e_data2%>%filter(.,Enzyme=="P")%>%arrange(Template,Length2,Cycles)

e_Pdata_FC <- e_Pdata%>%
  mutate(.,FC=(e_ratio_mean)/(e_Pdata[1,5]%>%as.numeric()))



e_P_his <- ggplot(data=e_Pdata_FC,aes(x=Template,y=e_ratio_mean,fill=Template))+
  geom_bar(stat = "identity")+
  geom_errorbar(aes(ymin = e_ratio_mean-e_ratio_se, ymax = e_ratio_mean+e_ratio_se,color=Template), width = 0.3) +
  labs(x=NULL,y="minor/all")+ 
  theme_classic()+
  facet_wrap(~Length2+Cycles,scales = "fixed",ncol = 9, nrow = 1,
             labeller = labeller(Length2 = function(x) paste0(x, " bp"),
                                 Cycles = function(x) paste0(x, " cycles")))+
  scale_fill_discrete(labels = c("50ng/100μl","500ng/100μl"))+
  scale_color_discrete(labels = c("50ng/100μl","500ng/100μl"))+
  scale_y_break(c(0.06, 0.15),scales = 1)+
  #scale_y_break(c(0.4, 1),scales = 1)+
  scale_y_continuous(sec.axis=sec_axis(~./(e_Pdata[1,5]%>%as.numeric()),name="Fold Change",breaks = c(1,10,20,90,120,150)))+
  geom_hline(yintercept = e_Pdata[1,5]%>%as.numeric(), col="#404040",lty = 2)+
  labs(title="fig1_e_Phusion_se")+
  my_theme+
  theme(
    legend.position = "bottom"#图例位置
  )

ggsave(e_P_his,filename = paste0(outpath,"fig1_e_Phusion_se.pdf"), 
       device = cairo_pdf, 
       width = 14,
       height = 6)



#Fig1_B:#########################

b_data <- read.table("fig1_b_all_data.csv",header = T,sep = ",")

b_data$Template <- factor(b_data$Template)


##Q5

b_Qdata<-b_data%>% filter(.,Enzyme=="Q")%>%subset(.,select=c(-Sample))%>% filter(.,One_umiToGeno_num <=5)
expanded_Qdata <- b_Qdata[rep(1:nrow(b_Qdata),b_Qdata$Umi_num),]

b_Q_vio_pic <- ggplot(expanded_Qdata,aes(x=Template,y=One_umiToGeno_num,fill=Template))+
  geom_violin(scale = "width",trim=FALSE)+
  facet_wrap(~Length2+Cycles,scales = "fixed",ncol = 9, nrow = 1,
             labeller = labeller(Length2 = function(x) paste0(x, " bp"),
                                 Cycles = function(x) paste0(x, " cycles")))+
  scale_fill_discrete(labels = c("50ng/100μl","500ng/100μl"))+
  labs(title="Fig1_b_Combind_repeat_Q5_umi2",x=NULL,y="The diversity of genotypes in each barcode")+ 
  theme_classic()+
  my_theme+
  theme(legend.position = "bottom")

ggsave(b_Q_vio_pic,filename = paste0(outpath,"Fig1_b_Combind_repeat_Q5_umi2.pdf"), 
       device = cairo_pdf,  
       width = 13,
       height = 6)





data_fig1_B <- b_data%>%
  filter(.,Enzyme=="Q")%>%
  mutate(.,Enzyme="Q5")%>%
  select(Template,Length2,Enzyme,Cycles,One_umiToGeno_num,Umi_num)


colnames(data_fig1_B) <- c("Templates","Length", "Enzyme","Cycles",
                           "The diversity of genotypes in each barcode",
                           "The number of barcodes")

write.csv(data_fig1_B,file=paste0(outpath,"data_fig1_B.csv"),row.names = F)





##Phusion
b_Pdata<-b_data%>% filter(.,Enzyme=="P")%>%subset(.,select=c(-Sample))%>% filter(.,One_umiToGeno_num <=5)
expanded_Pdata <- b_Pdata[rep(1:nrow(b_Pdata),b_Pdata$Umi_num),]

b_P_vio_pic <- ggplot(expanded_Pdata,aes(x=Template,y=One_umiToGeno_num,fill=Template))+
  geom_violin(scale = "width",trim=FALSE)+
  facet_wrap(~Length2+Cycles,scales = "fixed",ncol = 9, nrow = 1,
             labeller = labeller(Length2 = function(x) paste0(x, " bp"),
                                 Cycles = function(x) paste0(x, " cycles")))+
  scale_fill_discrete(labels = c("50ng/100μl","500ng/100μl"))+
  labs(title="Fig1_b_Combind_repeat_Phusion_umi2",x=NULL,y="The diversity of genotypes in each barcode")+ 
  theme_classic()+
  my_theme+
  theme(legend.position = "bottom")#图例位置

ggsave(b_P_vio_pic,filename = paste0(outpath,"Fig1_b_Combind_repeat_Phusion_umi2.pdf"), 
       device = cairo_pdf,  
       width = 13,
       height = 6)


#Fig1_D:##########

d_data <- read.table("fig1_d_all_data.csv",header = T,sep = ",")

d_data$Template <- factor(d_data$Template)


#Q5

d_Qdata <- d_data%>%filter(.,Enzyme=="Q")

d_Q_vio_pic <- ggplot(d_Qdata,aes(x=Template,y=Major_proportion,fill=Template))+
  geom_violin(scale = "width",trim=FALSE)+
  facet_wrap(~Length2+Cycle,scales = "fixed",ncol = 9, nrow = 1,
             labeller = labeller(Length2 = function(x) paste0(x, " bp"),
                                 Cycle = function(x) paste0(x, " cycles")))+
  scale_fill_discrete(name="Template",labels = c("50ng/100μl","500ng/100μl"))+#
  scale_y_continuous(breaks = c(0,0.25,0.5,0.75,1.0,1.25))+
  labs(title="Fig1_d_Combind_repeat_Q5_umi2",x=NULL,y="The frequency of HiFi reads of major genotype in each barcode")+ 
  theme_classic()+
  my_theme+
  theme(legend.position = "bottom")


ggsave(d_Q_vio_pic,filename = paste0(outpath,"Fig1_d_Combind_repeat_Q5_umi2.pdf"), 
       device = cairo_pdf,  
       width = 13,
       height = 6)



data_fig1_D <- d_Qdata%>%
  mutate(.,Enzyme="Q5")%>%
  select(Template,Length2,Enzyme,Cycle,Umi,Major_proportion)


colnames(data_fig1_D) <- c("Templates","Length", "Enzyme","Cycles",
                           "Barcode",
                           "The frequency of HiFi reads of major genotype in each barcode")

write.csv(data_fig1_D,file=paste0(outpath,"data_fig1_D.csv"),row.names = F)








#Phusion

d_Pdata <- d_data%>%filter(.,Enzyme=="P")


d_P_vio_pic <- ggplot(d_Pdata,aes(x=Template,y=Major_proportion,fill=Template))+
  geom_violin(scale = "width",trim=FALSE)+
  facet_wrap(~Length2+Cycle,scales = "fixed",ncol = 9, nrow = 1,
             labeller = labeller(Length2 = function(x) paste0(x, " bp"),
                                 Cycle = function(x) paste0(x, " cycles")))+
  scale_fill_discrete(name="Template",labels = c("50ng/100μl","500ng/100μl"))+
  scale_y_continuous(breaks = c(0,0.25,0.5,0.75,1.0,1.25))+
  labs(title="Fig1_d_Combind_repeat_Phusion_umi2",x=NULL,y="The frequency of HiFi reads of major genotype in each barcode")+ 
  theme_classic()+
  my_theme+
  theme(legend.position = "bottom")


ggsave(d_P_vio_pic,filename = paste0(outpath,"Fig1_d_Combind_repeat_Phusion_umi2.pdf"), 
       device = cairo_pdf,  
       width = 13,
       height = 6)

#Fig1_G:####################


raw_1000 <- read.table("1000_ComRepeat_switch_error_v5_umi_2.txt",header = TRUE)
raw_1000<- raw_1000[!(raw_1000$switch_rate_and_error_rate==0 & raw_1000$error_rate==0 & raw_1000$other_type!=0),]


raw_1000_data <- raw_1000%>%
  subset(.,select=c(Template,Length,Enzyme,Cycles,major_loc,switch_rate_and_error_rate,error_rate,switch_rate))%>%
  mutate(Cycles2 = case_when(
    Cycles==20 ~ 18.56,
    Cycles==30 ~ 21.87,
    Cycles==40 ~ 22.35))%>%
  dplyr::mutate(.,s_e_rate_per_cycle=switch_rate_and_error_rate/Cycles2,
                error_rate_per_cycle=error_rate/Cycles2,
                switch_rate_per_cycle=switch_rate/Cycles2)




Q_1000_data <- raw_1000_data%>%
  filter(Length==1000,Enzyme=="Q")%>%
  subset(.,select=c(-switch_rate_and_error_rate,-error_rate,-switch_rate))

Q_1000_data2 <- Q_1000_data%>%
  melt(.,id.vars = c("Template","Length","Enzyme","Cycles","major_loc"),
       measure.vars = c("error_rate_per_cycle","switch_rate_per_cycle"))


Q_pic <- ggplot(Q_1000_data2,aes(x=variable,y=value,fill=variable))+
  geom_boxplot(outlier.alpha = 0)+
  coord_cartesian(ylim = c(-0.00015,0.00023))+
  facet_wrap(~Template+Cycles,scales = "fixed",ncol = 6,
             labeller = labeller(Template = function(x) paste0(x, " ng"),
                                 Cycles = function(x) paste0(x, " cycles")))+
  labs(x=NULL,y=NULL,title="Combind_repeat_Q5_1000bp_per_cycle(qPCR)_use_freq_v5")+
  scale_fill_discrete(labels=c("NMR","TSR"))+
  scale_x_discrete(labels=c("NMR","TSR"))+
  stat_summary(fun.data = function(x) data.frame(y=0.0002, label = paste("N =", length(x))), geom="text", size = 5)+
  stat_summary(fun.data = function(x) data.frame(y=0.00015, label = ifelse(median(x)==0,"0",format(median(x), digits = 2,scientific = TRUE))), geom="text", size = 5)+
  theme_classic()+
  my_theme

ggsave(Q_pic,filename = paste0(outpath,"Fig1_g_Combind_repeat_Q5_umi2_use_freq_v5_2.pdf"), 
       width = 14,
       height = 6)



data_fig1_G <- Q_1000_data%>%
  mutate(.,Enzyme="Q5")%>%
  select(Template,Length,Enzyme,Cycles,Cycles2,error_rate_per_cycle,switch_rate_per_cycle)


colnames(data_fig1_G) <- c("Templates","Length", "Enzyme","Cycles","Effective cycles",
                           "NMR","TSR")

write.csv(data_fig1_G,file=paste0(outpath,"data_fig1_G.csv"),row.names = F)






#Phusion

P_1000_data <- raw_1000_data%>%
  filter(Length==1000,Enzyme=="P")%>%
  subset(.,select=c(-switch_rate_and_error_rate,-error_rate,-switch_rate))

P_1000_data2 <- P_1000_data%>%
  melt(.,id.vars = c("Template","Length","Enzyme","Cycles","major_loc"),
       measure.vars = c("error_rate_per_cycle","switch_rate_per_cycle"))

P_pic <- ggplot(P_1000_data2,aes(x=variable,y=value,fill=variable))+
  geom_boxplot(outlier.alpha = 0)+
  coord_cartesian(ylim = c(-0.00015,0.00023))+
  facet_wrap(~Template+Cycles,scales = "fixed",ncol = 6,
             labeller = labeller(Template = function(x) paste0(x, " ng"),
                                 Cycles = function(x) paste0(x, " cycles")))+
  labs(x=NULL,y=NULL,title="Combind_repeat_Phusion_1000bp_per_cycle(qPCR)_use_freq_v5")+
  scale_fill_discrete(labels=c("NMR","TSR"))+
  scale_x_discrete(labels=c("NMR","TSR"))+
  stat_summary(fun.data = function(x) data.frame(y=0.0002, label = paste("N = ", length(x))), geom="text", size = 5)+
  stat_summary(fun.data = function(x) data.frame(y=0.00015, label = ifelse(median(x)==0,"0",format(median(x), digits = 2,scientific = TRUE))), geom="text", size = 5)+
  theme_classic()+
  my_theme


ggsave(P_pic,filename = paste0(outpath,"Fig1_g_Combind_repeat_Phusion_umi2_use_freq_v5_2.pdf"), 
       width = 14,
       height = 6)




