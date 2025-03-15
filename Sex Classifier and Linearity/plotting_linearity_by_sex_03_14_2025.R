
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
  library(CytoTRACE)
  library('glmGamPoi')
  library(scran)
  library(parallel)
  library(factoextra)
  library(readxl)
  library(factoextra)
  library(forcats)
  library(ggrepel)
  library(biomaRt)
  library(openxlsx)
  library(glmnet)  
  library(forcats)
  `%notin%` <- Negate(`%in%`)
  
}


linearity_data <- read_csv('Sex Classifier and Linearity/linearity_score_03_13_2025.csv')
sex_data <- read_csv('Sex Classifier and Linearity/sex_predictions_03_13_2025.csv')

joined_data <- linearity_data%>%
  right_join(sex_data, by = 'cluster')

continuity_sex_joined <- joined_data 

joined_data <- na.omit(joined_data)
#attach(continuity_sex_joined)
#cor.test(x =mean_continuum_score, y = mean_pred )
#detach(continuity_sex_joined)

continuity_sex_joined_plot <- ggplot(subset(continuity_sex_joined, cluster %notin% c(15,30)), aes(x = mean_continuum_score, y = mean_pred))+
  geom_point(aes(x = mean_continuum_score, y = mean_pred), size = 2, alpha =0.5)+
  geom_text_repel(aes(label = cluster), size =3, show.legend = F, color = 'darkred')+
  theme_classic()+
  theme(axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10), axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12))+
  labs(x = 'Mean Linearity', y = 'Mean Predicted Sex', color = 'Clusters', shape = 'Clusters')+
  xlim(min(joined_data$mean_continuum_score)- max(joined_data$se), max(joined_data$mean_continuum_score)+ max(joined_data$se))+
  ylim(min(joined_data$mean_pred)- max(joined_data$se_pred), max(joined_data$mean_pred)+ max(joined_data$se_pred))
continuity_sex_joined_plot


ggsave(plot = continuity_sex_joined_plot,
       file = "linearity_by_sex.svg",
       device = "svg",
       units = "in",
       width = 1.8,
       height = 2,
       path = "Bachelors Thesis/Plots/Figure 3")



