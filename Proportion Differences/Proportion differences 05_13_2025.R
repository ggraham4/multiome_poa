#Proportion DEGs 05 13 2025
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
  
  
  
}

`%notin%` <- Negate(`%in%`)

obj <- readRDS('A:/optimal_clustering_05_06_2025/nemo.orig_harmony.integration_all_testd_clusters.rds')

cells_in_d <- nrow(obj@meta.data[obj@meta.data$Status == 'D',])
cells_in_m <- nrow(obj@meta.data[obj@meta.data$Status == 'M',])
cells_in_f <- nrow(obj@meta.data[obj@meta.data$Status == 'F',])

#### GLMER APPROACH ####
library(emmeans)

glmer_data <- data.frame()
for(i in 0:26){
  if(FALSE){next}
  print(i)
  f_loop_data_cluster <- data.frame()
  for(f in unique(obj@meta.data$individual)){
    cells_in_fish <- nrow(obj@meta.data[obj@meta.data$individual==f,])
    cells_in_fish_cluster <- nrow(obj@meta.data[obj@meta.data$individual==f&
                                                  obj@meta.data$res0.8_50nn_40PC_45LSI==i,])
    Status <- unique(obj@meta.data$Status[obj@meta.data$individual==f])
    f_loop_data <- data.frame(cells_in_fish=cells_in_fish,
                              cells_in_fish_cluster=cells_in_fish_cluster,
                              individual =f,
                              Status=Status,
                              cluster = i)%>%
      subset(Status %in% c("M",'D','F'))
    f_loop_data_cluster <- rbind(f_loop_data, f_loop_data_cluster)
    
  }
  glmer_matrix_cluster <- matrix(NA, nrow(f_loop_data_cluster),2)
  glmer_matrix_cluster[,1] <- f_loop_data_cluster$cells_in_fish_cluster
  glmer_matrix_cluster[,2] <- f_loop_data_cluster$cells_in_fish - f_loop_data_cluster$cells_in_fish_cluster
  
  glmer_model <- glmer(glmer_matrix_cluster~f_loop_data_cluster$Status + (1|individual), family = binomial('logit'), data = f_loop_data_cluster)
  car::Anova(glmer_model, type = 'III')
  pairs <- as.data.frame(pairs(emmeans(glmer_model, 'Status'),adjust = 'tukey'))
  
  anova_p.value <- car::Anova(glmer_model, type = 'III')[2,3]
  
  d_m_estimate <- pairs$estimate[pairs$contrast=='D - M']
  d_f_estimate <- pairs$estimate[pairs$contrast=='D - F']
  f_m_estimate <- pairs$estimate[pairs$contrast=='F - M']
  
  
  d_m_p.value <- pairs$p.value[pairs$contrast=='D - M']
  d_f_p.value <- pairs$p.value[pairs$contrast=='D - F']
  f_m_p.value <- pairs$p.value[pairs$contrast=='F - M']
  
  singular <- isSingular(glmer_model)
  
  error <- ifelse(length(glmer_model@optinfo$conv$lme4$code) != 0, substr(glmer_model@optinfo$conv$lme4$messages, 1, 50), NA)
  
  glmer_output <- data.frame(cluster = i,
                             d_m_estimate=d_m_estimate,
                             d_f_estimate = d_f_estimate,
                             f_m_estimate = f_m_estimate,
                             anova_p.value = anova_p.value,
                             d_m_p.value = d_m_p.value,
                             d_f_p.value = d_f_p.value,
                             f_m_p.value = f_m_p.value,
                             singular = singular,
                             error = error
  )
  
  
  glmer_data <- rbind(glmer_data, glmer_output)                          
}
glmer_data$issignif = ifelse(glmer_data$anova_p.value <0.05, '*', NA)

view(glmer_data)


### proportion differences in clusters 
#>24 #### ok this is quickly becoming a very interesting cluster
#>22 # old cluster 27, same result as before
#> 9 #part of 0 and 1, interesting

DimPlot(obj, label = T, group.by = 'old_clusters_res0.4_wnn', reduction = 'harmony_wnn.umap')
DimPlot(obj, label = T, reduction = 'harmony_wnn.umap', group.by = 'res0.8_50nn_40PC_45LSI')
DimPlot(obj, label = T, reduction = 'harmony_wnn.umap', group.by = 'gabe_celltype_ids')

prop_cluster <- function(cluster){
  
  output_plot_data <- data.frame()
  for(i in(unique(obj@meta.data$individual))){
    cells_in_individual <- nrow(obj@meta.data[obj@meta.data$individual==i,])
    cells_in_individual_in_cluster <- nrow(obj@meta.data[obj@meta.data$individual==i & obj@meta.data$res0.8_50nn_40PC_45LSI == cluster,])
    proportion <- cells_in_individual_in_cluster/cells_in_individual
    status <- unique(obj@meta.data$Status[obj@meta.data$individual==i])
    individual <- i
    
    gen_data <- data.frame(cluster = cluster,
                           cells_in_individual= cells_in_individual,
                           cells_in_individual_in_cluster = cells_in_individual_in_cluster,
                           proportion = proportion,
                           Status = status,
                           individual = i)
    output_plot_data <- rbind(gen_data, output_plot_data)
  }
  output_plot_data$Status <- factor(output_plot_data$Status, levels = c('NRM','M','D','E','NF','F'))
  return(output_plot_data)
  
}

plot_prop_cluster <- function(cluster){
  prop_cluster <- function(cluster){
    
    output_plot_data <- data.frame()
    for(i in(unique(obj@meta.data$individual))){
      cells_in_individual <- nrow(obj@meta.data[obj@meta.data$individual==i,])
      cells_in_individual_in_cluster <- nrow(obj@meta.data[obj@meta.data$individual==i & obj@meta.data$res0.8_50nn_40PC_45LSI == cluster,])
      proportion <- cells_in_individual_in_cluster/cells_in_individual
      status <- unique(obj@meta.data$Status[obj@meta.data$individual==i])
      individual <- i
      
      gen_data <- data.frame(cluster = cluster,
                             cells_in_individual= cells_in_individual,
                             cells_in_individual_in_cluster = cells_in_individual_in_cluster,
                             proportion = proportion,
                             Status = status,
                             individual = i)
      output_plot_data <- rbind(gen_data, output_plot_data)
    }
    output_plot_data$Status <- factor(output_plot_data$Status, levels = c('NRM','M','D','E','NF','F'))
    return(output_plot_data)
    
  }
  
  output_plot_data <- prop_cluster(cluster)
  
  output_plot_data$factor <- ifelse(output_plot_data$Status == "NRM", 0, NA)
  output_plot_data$factor <- ifelse(output_plot_data$Status == "M", 1, output_plot_data$factor)
  output_plot_data$factor <- ifelse(output_plot_data$Status == "D", 2, output_plot_data$factor)
  output_plot_data$factor <- ifelse(output_plot_data$Status == "E", 3, output_plot_data$factor)
  output_plot_data$factor <- ifelse(output_plot_data$Status == "NF", 4, output_plot_data$factor)
  output_plot_data$factor <- ifelse(output_plot_data$Status == "F", 5, output_plot_data$factor)
  
  plot <- ggplot(output_plot_data, aes(x = fct_reorder(individual, factor), y = proportion, color = Status))+
    geom_boxplot(aes(color = Status, fill = Status, group = Status), alpha = 0.25, outlier.shape = NA)+
    geom_point()+
    theme_classic()+
    labs(x  ='FishID', y = 'Proportion', title = paste0('cluster_',cluster))+
    theme(axis.text.x = element_text(angle = -45, vjust = 1, hjust=0))
  return(plot)
  
}


plot_prop_cluster(24)

plot_prop_cluster(22)

plot_prop_cluster(9)

#i am very suspicious of these clusters as immature



