library(dplyr)
library(ggplot2)

my_theme <- theme(legend.position = "bottom",
                  plot.title=element_text(hjust=0.5,size=8),
                  axis.text = element_text(size = 6),
                  strip.text = element_text(size = 8),
                  axis.title = element_text(size = 8),
                  legend.title = element_text(size = 8),
                  legend.text = element_text(size = 6))  



umi3_reads <- read.table("HLAC_umi3_reads.txt",header = T)



freq_data <- umi3_reads %>%
  mutate(
    category = factor(
      case_when(
        reads == 1 ~ "1",
        reads == 2 ~ "2",
        reads >= 3 ~ "≥3"
      ),
      levels = c("1", "2", "≥3")
    )
  ) %>%
  group_by(Cycles, Biorepeat,category) %>%  
  summarise(
    count = n(),
    .groups = "drop_last"
  ) %>%
  mutate(
    percentage = round(count / sum(count) * 100, 1)  
  ) %>%
  ungroup()



freq_data2 <- freq_data%>%
  dplyr::group_by(Cycles,category)%>%
  summarise(percentage_mean=mean(percentage),
            percentage_sd=sd(percentage),
            percentage_se = sd(percentage) / sqrt(n()),  
            sample_size = n(),                 
            .groups = "drop")  





S11a_pic <-ggplot(freq_data2, aes(x = category, y = percentage_mean)) +
  geom_col(width = 0.6,fill = "#00BA38") +
  geom_errorbar(aes(ymin = percentage_mean-percentage_se, ymax = percentage_mean+percentage_se), width = 0.2,color="black") +
  facet_wrap(~Cycles,scales = "fixed",ncol = 3,
             labeller = labeller(Cycles = function(x) paste0(x, " cycles")))+
  labs(
    x = "The number of HiFi reads",
    y = "The diversity of 3' UMIs(%)",
    title = NULL
  ) +
  theme_classic()+
  my_theme


ggsave(S11a_pic,filename = paste0(outpath,"S11a_pic_bio_se.pdf"), 
       device = cairo_pdf,  
       width = 19,
       height = 6.5,units = "cm")




