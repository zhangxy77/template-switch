library(ggplot2)
library(plyr)
library(dplyr)
library(reshape2)
library(scales)
library(ggsignif)

my_theme <- theme(
  plot.title=element_text(hjust=0.5,size=8),
  axis.text = element_text(size = 6),
  strip.text = element_text(size = 8),
  axis.title = element_text(size = 8),
  legend.title = element_text(size = 8),
  legend.text = element_text(size = 6))  


#method1数据####################
m1_raw <- read.table("Separate_new_m1_major_geno_reads_filter_2_v3.txt",header = TRUE,sep = "\t")

m1_raw2 <- m1_raw%>%
  mutate(.,Freq=Major_umi_types_num/Major_umi_types_num_sum,
         Method="Method1")%>%
  subset(.,select=c(Ku_type,Day,Template,Biorepeat,Method,Genotype,Syn_Type,Freq))


m1_raw3 <- m1_raw2%>%
  reshape2::dcast(.,Ku_type+Day+Template+Method+Genotype+Syn_Type~Biorepeat,value.var = "Freq")%>%
  dplyr::rename(c("r1"="1","r2"="2"))



#method2数据####################
m15_raw <- read.table("Separate_new_m2_major_geno_reads_filter_2_v3.txt",header = TRUE,sep = "\t")

m15_raw2 <- m15_raw%>%
  mutate(.,Freq=Major_umi_reads_num/Major_umi_reads_num_sum,
         Method="Method2")%>%
  subset(.,select=c(Ku_type,Day,Template,Biorepeat,Method,Genotype,Syn_Type,Freq))

m15_raw3 <- m15_raw2%>%
  reshape2::dcast(.,Ku_type+Day+Template+Method+Genotype+Syn_Type~Biorepeat,value.var = "Freq")%>%
  dplyr::rename(c("r1"="1","r2"="2"))


#method3数据################
m2_raw <- read.table("Separate_Fitness_m2_pair_reads_filter_2.txt",header = TRUE,sep = "\t")

m2_raw2 <- m2_raw%>%
  mutate(.,Freq=Major_umi_reads_num/Major_umi_reads_num_sum,
         Method="Method3")%>%
  subset(.,select=c(Ku_type,Day,Template,Biorepeat,Method,Genotype,Syn_Type,Freq))

m2_raw3 <- m2_raw2%>%
  reshape2::dcast(.,Ku_type+Day+Template+Method+Genotype+Syn_Type~Biorepeat,value.var = "Freq")%>%
  dplyr::rename(c("r1"="1","r2"="2"))




#Fig_S5:###########################

cor_pic <- function(rawdata,name){

  cor_results <- rawdata %>%
    drop_na()%>%
    group_by(Day,Template) %>%
    summarise(
      min_r1 = min(r1),  
      R = cor.test(r1, r2, method = "pearson")$estimate,
      p_value = cor.test(r1, r2, method = "pearson")$p.value,
      .groups = "drop") %>%
    mutate(label = paste("italic(R) == ", ifelse(R>0.99,"0.99",R),
                         ifelse(p_value == 0, paste0(" * ',' ~ italic(P) < '","1e-100"),
                                paste0(" * ',' ~ italic(P) == '",format(signif(p_value, digits = 3), scientific = TRUE))),"'"))
  
  
  p <- ggplot(data=rawdata,aes(x=r1,y=r2))+
    geom_point(size = 0.2)+
    geom_abline(slope = 1,color="grey",linetype="dashed")+
    
    facet_wrap(~Day,scales = "fixed",
               labeller = labeller(Day = function(x) paste0("Day ",x)))+
    scale_x_log10(labels = scientific,limits = c(1e-6, 1),breaks = 10^seq(-6, 1,by=2))+#坐标轴采用科学计数法
    scale_y_log10(labels = scientific,limits = c(1e-6, 1),breaks = 10^seq(-6, 1,by=2))+
    geom_text(data = cor_results, aes(x = 1e-6, y = 1, label = label), hjust = 0, vjust = 1,size = 2,parse = TRUE,inherit.aes = FALSE) +#，inherit.aes = FALSE 避免继承全局映射
    labs(title=paste0("Method",name),
         x=paste0("Genotype frequency",name," in technical replicate1"),
         y=paste0("Genotype frequency",name," in technical replicate2"))+
    theme_bw()+
    my_theme+
    theme(#axis.text = element_text(size = 15),
    legend.position = "",
    legend.title = element_blank(),
    panel.grid = element_blank(),
    strip.background = element_rect(fill = "transparent",   
                                    color = "black"))
  
  return(p)
  
}

m1_pic <- cor_pic(m1_raw3%>%filter(.,Template=="50"),"1")
m2_pic <- cor_pic(m15_raw3%>%filter(.,Template=="50"),"2")
m3_pic <- cor_pic(m2_raw3%>%filter(.,Template=="50"),"3")


all <- m1_pic/m2_pic/m3_pic

ggsave(all,filename = paste0(output,"figS5_point0.2.pdf"), 
       width = 16,
       height = 24,units = "cm")




#Fig2_I:##############


pow <- function(x, y) {
  return(x^y)
}


load("/TUY/TUY_generation.RData")
generation_data <- TUY_data_g_sum
g <- mean(TUY_data_g_sum$sum_generation)


fitness_fun_biorepeat <- function(data,template2,biorepeat){
  data1 <- data%>%filter(.,Template==template2 & Biorepeat==biorepeat)

  F0 <- data1%>%filter(.,Day=="0"&Genotype=="WT")%>%.$Freq
  F7 <- data1%>%filter(.,Day=="7"&Genotype=="WT")%>%.$Freq
  result <- data1%>%
    filter(.,Genotype!="WT")%>%
    reshape2::dcast(.,Ku_type+Template+Biorepeat+Method+Genotype+Syn_Type~Day,value.var = "Freq")%>%
    dplyr::rename(c("f0"="0","f7"="7"))%>%
    mutate(.,F0=F0,F7=F7,g=g,fitness=pow((f7/f0)/(F7/F0),1/g))
  result$Mutants_count <- ifelse(result$Genotype == "WT", 0, sapply(strsplit(result$Genotype, " "), length))
  return(result)
}


m1_fitness <- rbind(fitness_fun_biorepeat(m1_raw2,"50","1"),fitness_fun_biorepeat(m1_raw2,"50","2"))
m15_fitness <- rbind(fitness_fun_biorepeat(m15_raw2,"50","1"),fitness_fun_biorepeat(m15_raw2,"50","2"))
m2_fitness <- rbind(fitness_fun_biorepeat(m2_raw2,"50","1"),fitness_fun_biorepeat(m2_raw2,"50","2"))


All_fitness <- rbind(m1_fitness,m15_fitness,m2_fitness)%>%
  select(Ku_type,Template,Biorepeat,Method,Genotype,Syn_Type,fitness,Mutants_count)



fig2_i_data <- All_fitness%>%filter(.,Template==50 & Mutants_count==1 & fitness>0)

fig2_i_n_result <- fig2_i_data%>%
  group_by(Biorepeat,Method,Syn_Type) %>%
  summarise(n_rows = n(), 
            .groups = "drop")

fig2_i_n_result$Syn_Type <- factor(fig2_i_n_result$Syn_Type,levels=c("synonymous","non_synonymous"))
fig2_i_n_result$Method <- factor(fig2_i_n_result$Method)


fig2_i_n_result2 <- fig2_i_n_result%>%
  group_by(Method,Syn_Type)%>%
  summarise(SD=sd(n_rows),
            Mean=mean(n_rows),
            .groups = "drop")



ebtop <- function(x) {
  mean(x) + (sd(x) / sqrt(length(x)))
}


ebbottom <- function(x) {
  mean(x) - (sd(x) / sqrt(length(x)))
}




method_pairs <- combn(levels(fig2_i_n_result$Method), 2, simplify = FALSE)


fig2_i_syn_n <- fig2_i_n_result%>%filter(.,Syn_Type=="synonymous")

fig2_i_syn_pic <- ggplot(fig2_i_syn_n, aes(x = Method, y = n_rows, fill = Syn_Type)) +
  stat_summary(geom = "bar",
               fun = mean,
               width = 0.5)+
  stat_summary(geom = "errorbar",
               fun.min = ebbottom,
               fun.max = ebtop,
               width=0.1,color="black")+
  scale_y_continuous(limits = c(0, 700))+
  facet_wrap(~Syn_Type, scales = "free_x", nrow = 1,
             labeller = labeller(Syn_Type = c("non_synonymous" = "Non-synonymous",
                                              "synonymous" = "Synonymous"))) +
  scale_fill_manual(values = c("synonymous" = "#00BFC4", "non_synonymous" = "#F8766D"),
                    labels = c("synonymous" = "Synonymous", "non_synonymous" = "Non-synonymous")) +
  ggsignif::geom_signif(
    comparisons = method_pairs,  
    map_signif_level = function(p) sprintf("italic(P) == %.2g", p),
    test = "t.test",             
    test.args="less",           
    tip_length = 0.05,           
    textsize = 2,               
    y_position = c(450,550,500),
    parse = TRUE
  )  +
  labs(#title = "fig2_I_t.test_onetail",
    x = NULL,
    y = "Number of genotypes",
    fill = NULL) +
  theme_classic() +
  my_theme+
  theme(legend.position = "")



fig2_i_nonsyn_n <- fig2_i_n_result%>%filter(.,Syn_Type=="non_synonymous")

fig2_i_nonsyn_pic <- ggplot(fig2_i_nonsyn_n, aes(x = Method, y = n_rows, fill = Syn_Type)) +
  stat_summary(geom = "bar",
               fun = mean,
               width = 0.5)+
  stat_summary(geom = "errorbar",
               fun.min = ebbottom,
               fun.max = ebtop,
               width=0.1,color="black")+
  scale_y_continuous(limits = c(0, 2200))+
  facet_wrap(~Syn_Type, scales = "free_x", nrow = 1,
             labeller = labeller(Syn_Type = c("non_synonymous" = "Non-synonymous",
                                              "synonymous" = "Synonymous"))) +
  scale_fill_manual(values = c("synonymous" = "#00BFC4", "non_synonymous" = "#F8766D"),
                    labels = c("synonymous" = "Synonymous", "non_synonymous" = "Non-synonymous")) +
  ggsignif::geom_signif(
    comparisons = method_pairs,  
    map_signif_level = function(p) sprintf("italic(P) == %.2g", p),
    test = "t.test",             
    test.args="less",           
    tip_length = 0.05,          
    textsize = 2,               
    y_position = c(1400,1800,1600), 
    parse = TRUE
  )  +
  labs(#title = "fig2_I_t.test_onetail",
    x = NULL,
    y = "Number of genotypes",
    fill = NULL) +
  theme_classic() +
  my_theme+
  theme(legend.position = "")  


fig2_i_all_pic <- fig2_i_syn_pic+fig2_i_nonsyn_pic

ggsave(fig2_i_all_pic,filename = paste0(output,"fig2_i_all_pic.pdf"), 
       width = 8.5,
       height = 4,units = "cm")



ggsave(fig2_i_all_pic,filename = paste0(output,"v2_se_fig2_i_all_pic.pdf"),
       width = 8.5,
       height = 4,units = "cm")



data_fic2_I <- fig2_i_n_result%>%
  mutate(Syn_Type = ifelse(Syn_Type == "non_synonymous", "non-synonymous", "synonymous"))

write.csv(data_fic2_I,file=paste0(output,"data_fic2_I.csv"))



