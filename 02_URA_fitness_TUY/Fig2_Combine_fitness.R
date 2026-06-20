library(colorspace)
library(ggplot2)
library(scales)
library(plyr)
library(dplyr)
library(reshape2, lib.loc = "/usr/local/lib64/R/library")
library(ggpubr)
library(patchwork)
library(ggtext)
library(purrr)
library(stringr, lib.loc = "/usr/local/lib64/R/library")



my_theme <- theme(
                  plot.title=element_text(hjust=0.5,size=8),
                  axis.text = element_text(size = 6),
                  strip.text = element_text(size = 8),
                  axis.title = element_text(size = 8),
                  legend.title = element_text(size = 8),
                  legend.text = element_text(size = 6))  



#method1数据####################

m1_raw <- read.table("Combine_new_m1_major_geno_reads_filter_2_v3.txt",header = TRUE,sep = "\t")%>%subset(.,select=c(-Biorepeat))
m1_raw2 <- m1_raw%>%
  mutate(.,Freq=Major_umi_types_num/Major_umi_types_num_sum,
         Method="M1")%>%
  subset(.,select=c(Ku_type,Day,Template,Method,Genotype,Syn_Type,Freq))


#method2数据################
m15_raw <- read.table("Combine_new_m2_major_geno_reads_filter_2_v3.txt",header = TRUE,sep = "\t")%>%subset(.,select=c(-Biorepeat))

m15_raw2 <- m15_raw%>%
  mutate(.,Freq=Major_umi_reads_num/Major_umi_reads_num_sum,
         Method="M2")%>%
  subset(.,select=c(Ku_type,Day,Template,Method,Genotype,Syn_Type,Freq))

#method3数据################
m2_raw <- read.table("Combine_Fitness_m2_pair_reads_filter_2.txt",header = TRUE,sep = "\t")%>%subset(.,select=c(-Biorepeat))

m2_raw2 <- m2_raw%>%
  mutate(.,Freq=Major_umi_reads_num/Major_umi_reads_num_sum,
         Method="M3")%>%
  subset(.,select=c(Ku_type,Day,Template,Method,Genotype,Syn_Type,Freq))


All <- rbind(m1_raw2,m15_raw2,m2_raw2)


All2 <- All%>%
  reshape2::dcast(.,Ku_type+Day+Template+Genotype+Syn_Type~Method,value.var = "Freq")

All3 <- All2 %>% 
  mutate(.,M1_log=log10(M1),
         M2_log=log10(M2),
         M3_log=log10(M3),)

#Fig2_C:#############

cc <- function(fn,x,y,dday,ttemplate){
  zz <- fn%>%filter(.,Day==dday & Template==ttemplate)
  z <- cor.test(zz[[x]],zz[[y]])
  z2 <- ifelse(z$estimate > 0.99, "0.99", ifelse(z$estimate < 0.99, round(z$estimate, digits = 2), z$estimate))
  pp <- z$p.value
  return(paste("R == ", z2,
               ifelse(pp == 0, paste0(" * ',' ~ italic(P) < '","1e-100"),
                      paste0(" * ',' ~ italic(P) == '",format(signif(pp, digits = 3), scientific = TRUE))),"'"))
}


plot_fun2_log <- function(fn,day1,template1,x_var, y_var,cutoff){
  fn_f <- fn%>%filter(.,Day==day1 & Template==template1)
  fn_f2 <- fn_f%>%mutate(
    position = case_when(
      !!sym(y_var) > (!!sym(x_var) + cutoff) ~ "above",
      !!sym(y_var) < (!!sym(x_var) - cutoff) ~ "below",
      TRUE ~ "between"
    ))
  

  red_count <- sum(fn_f2$position == "above")
  blue_count <- sum(fn_f2$position == "below")
  

  count_label <- paste0(
    "<span style='color:#F8766D'>", red_count, "</span>",
    " vs. ",
    "<span style='color:#00BFC4'>", blue_count, "</span>"
  )
  

  p_label <- ""
  combined_label <- count_label  
  

  if (red_count + blue_count > 0) {
    binom_test <- binom.test(
      x = red_count, 
      n = red_count + blue_count,
      p = 0.5  
    )
    
    

    p_label <- ifelse(binom_test$p.value < 1e-100, paste0("<i>P</i> < ","1e-100"),
                      paste0("<i>P</i> = ",format(signif(binom_test$p.value, digits = 3), scientific = TRUE)))
    
    
    combined_label <- paste(count_label,p_label,sep = "\n") %>%
      str_replace_all("\n", "<br/>")  
  }
  
  p1 <- fn_f2%>%
    ggplot(aes(x=!!sym(x_var),y=!!sym(y_var)))+

    geom_abline(intercept = cutoff, slope = 1, linetype = "dashed", color = "gray50") +
    geom_abline(intercept = -cutoff, slope = 1, linetype = "dashed", color = "gray50") +
    geom_abline(slope = 1,color="black",linetype="solid")+
    geom_point(data = subset(fn_f2, position == "between"), 
               color = "black", size = 0.5) +
    geom_point(data = subset(fn_f2, position == "above"), 
               color = "#F8766D", size = 0.5) +
    geom_point(data = subset(fn_f2, position == "below"), 
               color = "#00BFC4", size = 0.5) +
    annotate("text",
           x = Inf, y = -Inf,
        label = p_label,
       hjust = -0.1, vjust = 1.5,
      size = 4, color = "black") +
    annotate("richtext",
             x = 0, y = -4,
             label = combined_label,
             hjust = 1, vjust = 1,
             size = 3,
             fill = NA, label.color = NA, 
             label.padding = grid::unit(rep(0, 4), "mm")) +
    scale_x_continuous(limits = c(-6, 0))+
    scale_y_continuous(limits = c(-6, 0))+
    facet_wrap(~Day,scales = "free",
               labeller = labeller(Day = function(x) paste0("Day ",x)))+
    annotate("text",x = -6, y = 0,
             label = cc(fn,strsplit(x_var, "_")[[1]][1],strsplit(y_var, "_")[[1]][1],day1,template1), 
             size = 3,hjust=0,vjust=1,parse = TRUE)+
    labs(
      x = paste0("log10(Genotype frequency", gsub("[^0-9]", "", x_var),")"),   
      y = paste0("log10(Genotype frequency", gsub("[^0-9]", "", y_var),")")    
    ) +
    theme_classic()+
    my_theme
  
  return(p1)
  
}

markers <- c("M1_log", "M2_log", "M3_log")


combinations <- combn(markers, 2, simplify = FALSE)


#Day0

plots_0 <- map(combinations, ~ plot_fun2_log(All3,"0","50", .x[1], .x[2],0.3))


final_plot_0 <- wrap_plots(plots_0, nrow = 1) + 
  plot_annotation(title = "Pairwise Marker Comparisons (Day 0, Template 50,intercept 0.3)")


#Day7

plots_7 <- map(combinations, ~ plot_fun2_log(All3,"7","50", .x[1], .x[2],0.3))


final_plot_7 <- wrap_plots(plots_7, ncol = 3) + 
  plot_annotation(title = "Pairwise Marker Comparisons (Day 7, Template 50,intercept 0.3)")

final_plot <- final_plot_0 / final_plot_7

ggsave(final_plot,filename = paste0(output,"fig2_C_log10_0.3.pdf"), 
       width = 19,
       height = 8,
       units = "cm")


data_fic2_C <- All3%>%filter(.,Template=="50")%>%
  mutate(Syn_Type = ifelse(Syn_Type == "non_synonymous", "non-synonymous", "synonymous"))

data_fic2_C[is.na(data_fic2_C)] <- "-"

write.csv(data_fic2_C,file=paste0(output,"data_fic2_C.csv"),row.names = F)



#fitness############


pow <- function(x, y) {
  return(x^y)
}


load("/TUY/TUY_generation.RData")
generation_data <- TUY_data_g_sum
g <- mean(TUY_data_g_sum$sum_generation)


fitness_fun <- function(data,template2){
  data1 <- data%>%filter(.,Template==template2)

  F0 <- data1%>%filter(.,Day=="0"&Genotype=="WT")%>%.$Freq
  F7 <- data1%>%filter(.,Day=="7"&Genotype=="WT")%>%.$Freq
  result <- data1%>%
    filter(.,Genotype!="WT")%>%
    reshape2::dcast(.,Ku_type+Template+Method+Genotype+Syn_Type~Day,value.var = "Freq")%>%
    dplyr::rename(c("f0"="0","f7"="7"))%>%
    mutate(.,F0=F0,F7=F7,g=g,fitness=pow((f7/f0)/(F7/F0),1/g))
  result$Mutants_count <- ifelse(result$Genotype == "WT", 0, sapply(strsplit(result$Genotype, " "), length))
  return(result)
}

m1_fitness <- rbind(fitness_fun(m1_raw2,"5"),fitness_fun(m1_raw2,"50"))
m15_fitness <- rbind(fitness_fun(m15_raw2,"5"),fitness_fun(m15_raw2,"50"))
m2_fitness <- rbind(fitness_fun(m2_raw2,"5"),fitness_fun(m2_raw2,"50"))


m1_snp_geno <- m1_fitness%>%
  filter(.,Mutants_count==1&fitness>0)%>%
  subset(select=c(Template,Genotype))

filtered_m2 <- inner_join(m2_fitness, m1_snp_geno, by = c("Genotype", "Template"))%>%
  mutate(.,Method="M3_overlap_genotype")


All_fitness <- rbind(m1_fitness,m15_fitness,filtered_m2,m2_fitness)%>%
  select(Ku_type,Template,Method,Genotype,Syn_Type,fitness,Mutants_count)
All_fitness$Method <- factor(All_fitness$Method,levels = c("M1","M2","M3_overlap_genotype","M3"))

#FigS6:######################

fig_s6_fitness <- All_fitness %>%
  filter(.,Template==50 & Mutants_count==1 & fitness>0 & Method != "M3_overlap_genotype")%>%
  mutate(
    Method2 = case_when(
      Method == "M1" ~ "Method1",
      Method == "M2" ~ "Method2",
      Method == "M3" ~ "Method3"
    )
  )

fig_s6_fitness_data <- fig_s6_fitness%>%
  group_by(Method2) %>%
  summarise(min_fitness=min(fitness),
            Median=median(fitness),
            Mean=mean(fitness),
            SD=sd(fitness),
            .groups = "drop")%>%
  mutate(label=paste0("Median = ",format(round(Median, 3), nsmall = 3),"\n",
                      "Mean = ",format(round(Mean, 3), nsmall = 3),"\n",
                      "s.d. = ",format(round(SD, 3), nsmall = 3)))


fig_s6 <- fig_s6_fitness%>%
  ggplot(aes(x=fitness,y=after_stat(count / tapply(count,PANEL,sum)[PANEL])))+
  geom_histogram(bins=30,fill="#FB6F66",alpha=0.8,color="#FB6F66",linewidth = 0.5)+
  facet_wrap(~Method2,scales = "free")+
  geom_vline(xintercept = 1, linetype = "dashed", color = "darkgrey") +
  geom_vline(data = fig_s6_fitness_data,aes(xintercept = Median),  
             color = "black", linetype = "dashed") +
  geom_text(data = fig_s6_fitness_data, aes(x = 1.01, y = 0.15, label = label), hjust = -0.05, vjust = 1,size = 2) +
  scale_x_continuous(labels = label_number(accuracy = 0.001)) +
  labs(
       x=expression(Fitness(italic(w))),y="Fraction")+
  theme_classic()+
  my_theme


ggsave(fig_s6,filename = paste0(output,"fig_s6.pdf"), #文件名称及其类型，一般通过改变后缀生成相应格式的图片
       width = 18,#宽
       height = 6,units = "cm")


#Fig2_K:####################


pairwise_wilcox_test_fun1 <- function(df) {
  

  methods <- unique(df$Method2)
  combn(methods, 2, simplify = FALSE) %>%
    map_df(function(pair) {

      x <- df %>% filter(Method2 == pair[1]) %>% pull(fitness)
      y <- df %>% filter(Method2 == pair[2]) %>% pull(fitness)
      

      test <- wilcox.test(x, y, alternative = "two.sided")
      

      tibble(
        Method1_type = pair[1],
        Method2_type = pair[2],
        p_value = test$p.value
      )
    })
}
density_wilcox_results <- pairwise_wilcox_test_fun1(fig_s6_fitness)%>%
  mutate(
    p_adj = p.adjust(p_value, method = "BH"),
    method_comparison = paste(Method1_type, "vs.", Method2_type),
    p_label_two_sided = paste0("italic(P) == ",format(signif(p_value, digits = 3), scientific = TRUE)),
    p_label_one_sided = paste0("italic(P) == ",format(signif(p_value/2, digits = 3), scientific = TRUE))
  )



write.csv(density_wilcox_results,file=paste0(output,"fig2_k_p_value_wilcox.csv"))


fig2_k_2_v2 <- ggplot(data = fig_s6_fitness, aes(x = fitness, linetype = Method2)) +

  stat_density(geom = "line", position = "identity", linewidth = 0.8, color = "black") +
  scale_linetype_manual(
    values = c("Method1" = "solid", 
               "Method2" = "dashed", 
               "Method3" = "dotted")) +
  geom_vline(
    data = fig_s6_fitness_data, 
    aes(xintercept = Median, linetype = Method2), 
    color = "grey",
    linewidth = 0.8  
  ) +
  coord_cartesian(xlim = c(0.96,1.01))+
  labs(
    x = expression(Fitness(italic(w))), 
    y = "Density") +
  theme_classic() +
  my_theme +
  theme(
    legend.position = "",
    legend.title = element_blank(),
    panel.grid = element_blank(),
    strip.background = element_rect(
      fill = "transparent",
      color = "black"))


ggsave(fig2_k_2_v2,filename = paste0(output,"v2_fig2_k2.pdf"), 
       width = 5,#宽
       height = 3,units = "cm")




data_fic2_K <- fig_s6_fitness%>%
  mutate(Syn_Type = ifelse(Syn_Type == "non_synonymous", "non-synonymous", "synonymous"))%>%
  select(Method2,Genotype,Syn_Type,fitness)

write.csv(data_fic2_K,file=paste0(output,"data_fic2_K.csv"),row.names = F)



#Fig2_DEF:##############################

M123_overlap_data <- All_fitness%>%
  filter(.,Method!="M3")%>%
  filter(.,Template==50 & Mutants_count==1 & fitness>0)%>%
  mutate(
    Method2 = case_when(
      Method == "M1" ~ "Method1",
      Method == "M2" ~ "Method2",
      Method == "M3_overlap_genotype" ~ "Method3"
    )
  )

M123_overlap_data$Syn_Type <- factor(M123_overlap_data$Syn_Type,levels = c("synonymous","non_synonymous"))
M123_overlap_data$Method2 <- factor(M123_overlap_data$Method2)


wilcox_results <- M123_overlap_data%>%
  group_by(Method2) %>%
  summarise(
    p_value=wilcox.test(fitness[Syn_Type=="synonymous"],fitness[Syn_Type=="non_synonymous"], alternative = "greater")$p.value,
    .groups = "drop"
  )%>%
  mutate(label=paste("italic(P) == ", format(round(p_value, 3), nsmall = 3)))


n_result <- M123_overlap_data%>%
  group_by(Method2) %>%
  summarise(n_rows = n(), 
            .groups = "drop")%>%
  mutate(label=paste("N == ", n_rows))


fig2_def <- ggplot(M123_overlap_data,aes(x=fitness,color=Syn_Type, linetype = Method2))+
  stat_ecdf(geom="smooth")+
  facet_wrap(~Method2,scales = "free")+
  coord_cartesian(xlim = c(0.96,1.01))+
  scale_color_manual(
    values = c("synonymous" = "#00BFC4","non_synonymous" = "#F8766D"),  
    labels = c("Synonymous","Non-synonymous")
  ) +
  scale_linetype_manual(
    values = c("Method1" = "solid", "Method2" = "dashed", "Method3" = "dotted") 
  ) +
  geom_text(data = wilcox_results, aes(x = 0.96, y = 1, label = label), hjust = 0, vjust = 1,size = 5,
            parse = TRUE,inherit.aes = FALSE) +
  geom_text(data = n_result, aes(x = 0.96, y = 0.9, label = label), hjust = 0, vjust = 1,size = 5,
            parse = TRUE,inherit.aes = FALSE) +
  labs(x=expression(Fitness(italic(w))),y="Cumulative frequency")+
  theme_bw()+
  my_theme+
  theme(#axis.text = element_text(size = 15),
    legend.position = "",
    legend.title = element_blank(),
    panel.grid = element_blank(),
    strip.background = element_rect(fill = "transparent",  
                                    color = "black"))



ggsave(fig2_def,filename = paste0(output,"fig2_def.pdf"),
       width = 19,#宽
       height = 4,
       units = "cm")


ggsave(fig2_def,filename = paste0(output,"v2_fig2_def.pdf"), 片
       width = 19,#宽
       height = 4,
       units = "cm")


data_fic2_DEFGH <- M123_overlap_data%>%
  mutate(Syn_Type = ifelse(Syn_Type == "non_synonymous", "non-synonymous", "synonymous"))%>%
  select(Method2,Genotype,Syn_Type,fitness)

write.csv(data_fic2_DEFGH,file=paste0(output,"data_fic2_DEFGH.csv"),row.names = F)




#Fig2_J:###############

M3_fitness_data <- All_fitness%>%
  filter(.,Method=="M3")%>%
  filter(.,Template==50 & Mutants_count==1 & fitness>0)%>%
  mutate(
    Method2 = case_when(
      Method == "M3" ~ "Method3"
    )
  )

M3_fitness_data$Syn_Type <- factor(M3_fitness_data$Syn_Type)
M3_fitness_data$Method2 <- factor(M3_fitness_data$Method2)


wilcox_results_M3 <- M3_fitness_data%>%
  group_by(Method2) %>%
  summarise(
    p_value=wilcox.test(fitness[Syn_Type=="synonymous"],fitness[Syn_Type=="non_synonymous"], alternative = "greater")$p.value,
    .groups = "drop"
  )%>%
  mutate(label=paste("italic(P) == ", format(round(p_value, 3), nsmall = 3)))


n_result_M3 <- M3_fitness_data%>%
  group_by(Method2) %>%
  summarise(n_rows = n(), 
            .groups = "drop")%>%
  mutate(label=paste("N == ", n_rows))


fig2_j <- ggplot(M3_fitness_data,aes(x=fitness,color=Syn_Type, linetype = Method2))+
  stat_ecdf(geom="smooth")+
  facet_wrap(~Method2,scales = "free")+
  coord_cartesian(xlim = c(0.96,1.01))+
  geom_text(data = wilcox_results_M3, aes(x = 0.96, y = 1, label = label), hjust = 0, vjust = 1,size = 5,
            parse = TRUE,inherit.aes = FALSE) +
  geom_text(data = n_result_M3, aes(x = 0.96, y = 0.9, label = label), hjust = 0, vjust = 1,size = 5,
            parse = TRUE,inherit.aes = FALSE) +
  scale_color_discrete(labels = c("Non-synonymous","Synonymous"))+
  scale_linetype_manual(
    values = c("Method1" = "solid", "Method2" = "dashed", "Method3" = "dotted")  
  ) +
  labs(x=expression(Fitness(italic(w))),y="Cumulative frequency")+
  theme_bw()+
  my_theme+
  theme(#axis.text = element_text(size = 15),
    legend.position = "",
    legend.title = element_blank(),#隐藏图例标题
    panel.grid = element_blank(),
    strip.background = element_rect(
      fill = "transparent",    # 透明底色
      color = "black"))


ggsave(fig2_j,filename = paste0(output,"v2_fig2_j.pdf"),
       width = 6.3,#宽
       height = 4,
       units = "cm")



data_fic2_J <- M3_fitness_data%>%
  mutate(Syn_Type = ifelse(Syn_Type == "non_synonymous", "non-synonymous", "synonymous"))%>%
  select(Method2,Genotype,Syn_Type,fitness)

write.csv(data_fic2_J,file=paste0(output,"data_fic2_J.csv"),row.names = F)




#Fig2_GH:##############


pairwise_wilcox_test <- function(data, syn_type) {

  df <- data %>% filter(Syn_Type == syn_type)
  

  methods <- unique(df$Method2)
  combn(methods, 2, simplify = FALSE) %>%
    map_df(function(pair) {

      x <- df %>% filter(Method2 == pair[1]) %>% pull(fitness)
      y <- df %>% filter(Method2 == pair[2]) %>% pull(fitness)
      

      test <- wilcox.test(x, y, alternative = "two.sided")
      

      tibble(
        Syn_Type = syn_type,
        Method1_type = pair[1],
        Method2_type = pair[2],
        p_value = test$p.value
      )
    })
}



syn_wilcox_results <- bind_rows(
  pairwise_wilcox_test(M123_overlap_data, "synonymous"),
  pairwise_wilcox_test(M123_overlap_data, "non_synonymous")
) %>%
  mutate(
    p_adj = p.adjust(p_value, method = "BH"), 
    method_comparison = paste(Method1_type, "vs.", Method2_type),
    p_label_two_sided = paste("italic(P) == ", format(round(p_value, 3), nsmall = 3)),
    p_label_one_sided = paste("italic(P) == ", format(round(p_value/2, 3), nsmall = 3))
  )

write.csv(syn_wilcox_results,file=paste0(output,"fig2_gh_p_value_wilcox.csv"))




fig2_gh <-ggplot(M123_overlap_data,aes(x=fitness,color=Syn_Type, linetype = Method2))+
  stat_ecdf(geom="smooth")+
  facet_wrap(~Syn_Type,scales = "free",
             labeller = labeller(Syn_Type = c(
               "non_synonymous" = "Non-synonymous",
               "synonymous" = "Synonymous"
             )))+
  coord_cartesian(xlim = c(0.96,1.01))+
  scale_color_manual(
    values = c("synonymous" = "#00BFC4","non_synonymous" = "#F8766D"), 
    labels = c("Synonymous","Non-synonymous")
  ) +
  scale_linetype_manual(
    values = c("Method1" = "solid", "Method2" = "dashed", "Method3" = "dotted")  
  ) +
  labs(x=expression(Fitness(italic(w))),y="Cumulative frequency")+
  theme_bw()+
  my_theme+
  theme(#axis.text = element_text(size = 15),
    legend.position = "right",
    legend.title = element_blank(),
    panel.grid = element_blank(),
    strip.background = element_rect(
      fill = "transparent",    
      color = "black"))


ggsave(fig2_gh,filename = paste0(output,"v2_fig2_gh.pdf"), 
       width = 19,#宽
       height = 4,units = "cm")





