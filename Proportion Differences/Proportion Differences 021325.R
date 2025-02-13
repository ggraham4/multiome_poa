#Proportion DEGs 122724
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
  plot_prop_cluster <- function(cluster){
    
    output_plot_data <- data.frame()
    for(i in(unique(obj@meta.data$individual))){
      cells_in_individual <- nrow(obj@meta.data[obj@meta.data$individual==i,])
      cells_in_individual_in_cluster <- nrow(obj@meta.data[obj@meta.data$individual==i & obj@meta.data$harmony.wnn_res0.4_clusters == cluster,])
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
  
  prop_cluster <- function(cluster){
    
    output_plot_data <- data.frame()
    for(i in(unique(obj@meta.data$individual))){
      cells_in_individual <- nrow(obj@meta.data[obj@meta.data$individual==i,])
      cells_in_individual_in_cluster <- nrow(obj@meta.data[obj@meta.data$individual==i & obj@meta.data$harmony.wnn_res0.4_clusters == cluster,])
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
  
  
}

obj <- readRDS("C:/Users/Gabe/Desktop/RNA Object.rds")

cells_in_d <- nrow(obj@meta.data[obj@meta.data$Status == 'D',])
cells_in_m <- nrow(obj@meta.data[obj@meta.data$Status == 'M',])
cells_in_f <- nrow(obj@meta.data[obj@meta.data$Status == 'F',])

#### GLMER APPROACH ####
library(emmeans)

glmer_data <- data.frame()
for(i in 0:31){
  print(i)
  f_loop_data_cluster <- data.frame()
  for(f in unique(obj@meta.data$individual)){
    cells_in_fish <- nrow(obj@meta.data[obj@meta.data$individual==f,])
    cells_in_fish_cluster <- nrow(obj@meta.data[obj@meta.data$individual==f&
                                                  obj@meta.data$harmony.wnn_res0.4_clusters==i,])
    Status <- unique(obj@meta.data$Status[obj@meta.data$individual==f])
    f_loop_data <- data.frame(cells_in_fish=cells_in_fish,
                              cells_in_fish_cluster=cells_in_fish_cluster,
                              individual =f,
                              Status=Status,
                              cluster = i)
    f_loop_data_cluster <- rbind(f_loop_data, f_loop_data_cluster)
    
  }
  glmer_matrix_cluster <- matrix(NA, nrow(f_loop_data_cluster),2)
  glmer_matrix_cluster[,1] <- f_loop_data_cluster$cells_in_fish_cluster
  glmer_matrix_cluster[,2] <- f_loop_data_cluster$cells_in_fish - f_loop_data_cluster$cells_in_fish_cluster
  
  glmer_model <- glmer(glmer_matrix_cluster~f_loop_data_cluster$Status + (1|individual), family = binomial('logit'), data = f_loop_data_cluster)
  car::Anova(glmer_model, type = 'III')
  pairs <- as.data.frame(pairs(emmeans(glmer_model, 'Status'),adjust = 'none'))
  
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
glmer_data$anova_q.value <- p.adjust(glmer_data$anova_p.value, 'fdr',nrow(glmer_data))


glmer_data$d_m_q.value <- p.adjust(glmer_data$d_m_p.value, 'fdr',nrow(glmer_data))
glmer_data$d_f_q.value <- p.adjust(glmer_data$d_f_p.value, 'fdr',nrow(glmer_data))
glmer_data$f_m_q.value <- p.adjust(glmer_data$f_m_p.value, 'fdr',nrow(glmer_data))

### q value is very conservative
glmer_data$issignif <- ifelse(apply(glmer_data[, 12:14], 1, function(x) any(x < 0.05)), '*', NA)
glmer_signif <- glmer_data$issignif
names(glmer_signif) <- 0:31
glmer_signif

### lets try with p value
glmer_data_p <- glmer_data
glmer_data_p$issignif <- ifelse(apply(glmer_data_p[, 6:8], 1, function(x) any(x < 0.05)), '*', NA)
glmer_p_signif <- glmer_data_p$issignif
names(glmer_p_signif) <- 0:31
glmer_p_signif

### means and ses
meta = obj@meta.data
total_cells =meta%>%
  group_by(Status)%>%
  summarize(cells = n()
  )

cluster_cells =meta%>%
  group_by(harmony.wnn_res0.4_clusters, Status)%>%
  summarize(cells = n()
  )%>%
  full_join(total_cells, by = 'Status')%>%
  mutate(prop = cells.x/cells.y)%>%
  mutate(
         SD = prop*(1-prop)
  )%>%
  subset(Status %in% c('M','F','D'))

plot_prop_cluster(3)
plot_prop_cluster(9)
plot_prop_cluster(27)


define_degs_glm <- function(data, singular = TRUE) {
  if (!singular) {
    data <- data[data$singular == FALSE, ]
  }
  
  # Assign classes based on conditions
  data$class <- NA  # Initialize class column
  
  data$class[data$d_m_q.value < 0.05 & 
               data$d_m_estimate > 0 & 
               data$d_f_q.value > 0.05] <- 'Early Upregulated'
  
  
  data$class[data$d_m_q.value > 0.05 & 
               data$d_f_q.value < 0.05 & 
               data$d_f_estimate < 0] <- 'Late Upregulated'
  
  data$class[data$d_m_q.value < 0.05 & 
               data$d_m_estimate < 0 & 
               data$d_f_q.value > 0.05] <- 'Early Downregulated'
  
  
  data$class[data$d_m_q.value > 0.05 & 
               data$d_f_q.value < 0.05 & 
               data$d_f_estimate > 0] <- 'Late Downregulated'
  
  data$class[data$d_m_q.value < 0.05 & 
               data$d_f_q.value < 0.05 & 
               data$d_f_estimate > 0 & 
               data$d_m_estimate > 0] <- 'Transiently Upregulated'
  
  data$class[data$d_m_q.value < 0.05 & 
               data$d_f_q.value < 0.05 & 
               data$d_f_estimate < 0 & 
               data$d_m_estimate < 0] <- 'Transiently Downregulated'
  
  data$class[data$d_m_q.value < 0.05 & 
               data$d_f_q.value < 0.05 & 
               data$d_m_estimate < 0 & 
               data$d_f_estimate > 0] <- 'Progressively Downregulated'
  
  data$class[data$d_m_q.value < 0.05 & 
               data$d_f_q.value < 0.05 & 
               data$d_m_estimate > 0 & 
               data$d_f_estimate < 0] <- 'Progressively Upregulated'
  
  data$class[data$f_m_q.value < 0.05 & 
               data$d_f_q.value > 0.05 & 
               data$d_m_q.value > 0.05 & 
               data$f_m_estimate > 0 ] <- 'Terminally Upregulated'
  
  data$class[data$f_m_q.value < 0.05 & 
               data$d_f_q.value > 0.05 & 
               data$d_m_q.value > 0.05 & 
               data$f_m_estimate < 0  ] <- 'Terminally Downregulated'
  
  data$issignif <- ifelse(data$f_m_q.value<0.05|
                            data$d_m_q.value<0.05|
                            data$d_f_q.value<0.05,
                          '*',NA)
  
  return(data)
}

glm_data_class <- define_degs_glm(glmer_data)
glm_data_class[,c(1,ncol(glm_data_class)-1, ncol(glm_data_class))]

plot_prop_cluster(6)

