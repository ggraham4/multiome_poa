library(lme4)
library(tidyverse)
library(Seurat)
library(car)
library(emmeans)

obj <- readRDS("C:/Users/Gabe/Desktop/RNA Object.rds")

obj_data <- obj@assays$RNA$data[,obj$harmony.wnn_res0.4_clusters =='2']%>%
  t()%>%
  as.data.frame()%>%
  dplyr::select('LOC111577263')

aromatase_data <- data.frame(aromatase_expression = obj_data,
                             individual = obj$individual[obj$harmony.wnn_res0.4_clusters =='2'],
                             Sex = obj$Status[obj$harmony.wnn_res0.4_clusters =='2'])%>%
  group_by(individual, Sex)%>%
  summarize(mean_expression = mean(LOC111577263),
            se = sd(LOC111577263)/sqrt(n()))%>%
  subset(Sex != 'NRM')
library(ggsignif)

aromatase_data$Sex <- factor(aromatase_data$Sex, levels = c('M','D','E',"NF",'F'))
aromatase_plot <- ggplot(aromatase_data, aes(x = Sex, y = mean_expression))+
  geom_boxplot(alpha = 0.25, outlier.shape = NA, aes(color = Sex, fill = Sex))+
  geom_jitter(  shape = 1, color = 'black', size =2)+
  theme_classic()+
  theme(legend.position ='none',plot.title = element_text(hjust=0.5))+
  labs(title = 'Cluster 2', y = 'Mean cyp19a1b')+
  geom_signif(xmin = c(2.1), xmax = c(5), y_position = c(2.6), annotation =c("**"), color = "black", tip_length = c(0,0), textsize=6)+
  ylim(0.45,3)
aromatase_plot

ggsave(plot = aromatase_plot,
       file = "aromatase_expression.svg",
       device = "svg",
       units = "in",
       width = 1.6,
       height = 2,
       path = "Bachelors Thesis/Plots/Figure 3")

