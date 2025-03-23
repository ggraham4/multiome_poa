# sex differences in radial glia WGCNA
{
  library(parallel)
  library(clusterProfiler)
  library(blme)
  library(Seurat)
  library(tidyverse)
  library(tidyr)
  library(lme4)
  library(dplyr)
  library(MASS)
  library(SeuratObject)
  library(Signac)
  library(glmGamPoi)
  library(scran)
  library(parallel)
  library(factoextra)
  library(readxl)
  library(factoextra)
  library(forcats)
  library(ggrepel)
  library(biomaRt)
  library(openxlsx)
  library(emmeans)
  library(CytoTRACE)
  library(ggrepel)
  library(tidyverse)
  library(cowplot)
  library(patchwork)
  library(Polychrome)
  library(scCustomize)
  library(hdWGCNA)
  
  P40 <- createPalette(40, c("#FF0000", "#00FF00", "#0000FF"), range = c(30, 80))
  swatch(P40)
  names(P40) <- NULL
  
  mean_expression_cluster_plot<- readRDS('Functions/mean_expression_cluster_plot.rds')
  prop_cluster_plot<- readRDS( 'Functions/prop_cluster_plot.rds')
  define_degs_prop<- readRDS('Functions/define_degs_prop.rds')
  mean_expression_cluster_data<- readRDS('Functions/mean_expression_cluster_data.rds')
  prop_deg_function.rds<- readRDS('Functions/DEG_functions/prop_deg_function.rds')
  define_behavior_degs<- readRDS('Functions/define_behavior_degs')
  clown_go<- readRDS('Functions/clown_go')
  define_degs<- readRDS('Functions/define_degs')
  
}

radial_glia <- readRDS(file='A:/anemonefish_multiome_radial_glia_03_20_2025/radial_glia_object.rds')

# module eigengenes:
MEs <- GetMEs(radial_glia, harmonized=FALSE)

MEs$individual = radial_glia$individual
MEs$Status = radial_glia$Status
MEs$sub = radial_glia$sub

MEs_model <- MEs%>%
  subset(Status %in%c('M','D','F'))

###brown
brown_model = lmer(brown~Status*sub +(1|individual), data = MEs_model)
car::Anova(brown_model, type = 'III')
pairs(emmeans(brown_model, 'Status', by ='sub'), adjust ='none') #DF differ in 2_4, everyone in 2_1

### blue
blue_model = lmer(blue~Status*sub +(1|individual), data = MEs_model)
car::Anova(blue_model, type = 'III') # no difference

### yellow 
yellow_model = lmer(yellow~Status*sub +(1|individual), data = MEs_model)
car::Anova(yellow_model, type = 'III') 
pairs(emmeans(yellow_model, 'Status', by ='sub'), adjust ='none') # 2_2 df and dm

### turquoise
turquoise_model = lmer(turquoise~Status*sub +(1|individual), data = MEs_model)
car::Anova(turquoise_model, type = 'III') 
pairs(emmeans(turquoise_model, 'Status', by ='sub'), adjust ='none')## 2_7 but who cares
#  significant difference in 2_3


### Plot brown differences those are going to be interesting
brown_data <- MEs%>%
  group_by(individual, Status, sub)%>%
  summarize(brown = mean(brown))%>%
  filter(Status != 'NRM')
brown_data
brown_data$Status <- factor(brown_data$Status, levels = c('M','D','E', 'NF','F'))

library(ggsignif)
brown_plot_2_4 <- ggplot(subset(brown_data, sub =='2_4'), aes(x = Status, y = brown))+
  geom_boxplot(aes(color = Status, fill = Status), alpha = 0.25, outlier.shape=NA)+
  theme_classic()+
  geom_jitter(color= 'black', shape =1, size = 2)+
  labs(x = 'Sex', fill = 'Sex', color = 'Sex', title = '2_4')+
  theme(legend.position = 'none', plot.title = element_text(hjust = 0.5),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank())+
  geom_signif(xmin = c(2.1), xmax = c(5), y_position = c(5.2), annotation =c("**"), color = "black", tip_length = c(0,0), textsize=6)+
  geom_signif(xmin = c(1), xmax = c(1.9), y_position = c(5.2), annotation =c("*"), color = "black", tip_length = c(0,0), textsize=6)+
  ylim(0,6)
brown_plot_2_4

#ggsave(plot = brown_plot_2_4,
#       file = "brown_2_4.svg",
#       device = "svg",
#       units = "in",
#       width = 1.25,
#       height = 2,
#       path = "Bachelors Thesis/Plots/Figure 5")

pairs(emmeans(brown_model, 'Status', by ='sub'), adjust ='none') #DF differ in 2_4, everyone in 2_1

brown_plot_2_1 <- ggplot(subset(brown_data, sub =='2_1'), aes(x = Status, y = brown))+
  geom_boxplot(aes(color = Status, fill = Status), alpha = 0.25, outlier.shape=NA)+
  theme_classic()+
  geom_jitter(color= 'black', shape =1, size = 2)+
  labs(x = 'Sex', fill = 'Sex', color = 'Sex', y = 'Brown Expression', title = '2_1')+
  theme(legend.position = 'none', plot.title = element_text(hjust = 0.5))+
  geom_signif(xmin = c(2.1), xmax = c(5), y_position = c(5.2), annotation =c("***"), color = "black", tip_length = c(0,0), textsize=6)+
  geom_signif(xmin = c(1), xmax = c(1.9), y_position = c(5.2), annotation =c("**"), color = "black", tip_length = c(0,0), textsize=6)+
  ylim(0,6)
brown_plot_2_1

#ggsave(plot = brown_plot_2_1,
#       file = "brown_plot_2_1.svg",
#       device = "svg",
#       units = "in",
#       width = 1.6,
#       height = 2,
#       path = "Bachelors Thesis/Plots/Figure 5")

turquoise_data <- MEs%>%
  group_by(individual, Status, sub)%>%
  summarize(turquoise = mean(turquoise))%>%
  filter(Status != 'NRM')
turquoise_data
turquoise_data$Status <- factor(turquoise_data$Status, levels = c('M','D','E', 'NF','F'))

pairs(emmeans(turquoise_model, 'Status', by ='sub'), adjust ='none') 

turq_plot_2_3 <- ggplot(subset(turquoise_data, sub =='2_3'), aes(x = Status, y = turquoise))+
  geom_boxplot(aes(color = Status, fill = Status), alpha = 0.25, outlier.shape=NA)+
  theme_classic()+
  geom_jitter(color= 'black', shape =1, size = 2)+
  labs(x = 'Sex', fill = 'Sex', color = 'Sex', y = 'Turquoise Expression', title = '2_3')+
  theme(legend.position = 'none', plot.title = element_text(hjust = 0.5))+
  geom_signif(xmin = c(2), xmax = c(5), y_position = c(15), annotation =c("*"), color = "black", tip_length = c(0,0), textsize=6)+
  ylim(min(subset(turquoise_data, sub =='2_3')$turquoise),17 )
turq_plot_2_3

ggsave(plot = turq_plot_2_3,
       file = "turq_plot_2_3.svg",
       device = "svg",
       units = "in",
       width = 1.6,
       height = 2,
       path = "Bachelors Thesis/Plots/Figure 5")



