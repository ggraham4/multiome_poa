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
  library(emmeans)

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

  define_degs <- function(data, singular = TRUE) {
  if (!singular) {
    data <- data[data$singular == FALSE, ]
  }
  
  # Assign classes based on conditions
  data$class <- NA  # Initialize class column
  
  data$class[data$`q_value_D - M` < 0.05 & 
             data$`OR_D - M` > 1 & 
             data$`q_value_D - F` > 0.05] <- 'Early Upregulated'
  
  
  data$class[data$`q_value_D - M` > 0.05 & 
             data$`q_value_D - F` < 0.05 & 
             data$`OR_D - F` < 1] <- 'Late Upregulated'
  
  data$class[data$`q_value_D - M` < 0.05 & 
             data$`OR_D - M` < 1 & 
             data$`q_value_D - F` > 0.05] <- 'Early Downregulated'
  
  
  data$class[data$`q_value_D - M` > 0.05 & 
             data$`q_value_D - F` < 0.05 & 
             data$`OR_D - F` > 1] <- 'Late Downregulated'
  
  data$class[data$`q_value_D - M` < 0.05 & 
             data$`q_value_D - F` < 0.05 & 
             data$`OR_D - F` > 1 & 
             data$`OR_D - M` > 1] <- 'Transiently Upregulated'
  
  data$class[data$`q_value_D - M` < 0.05 & 
             data$`q_value_D - F` < 0.05 & 
             data$`OR_D - F` < 1 & 
             data$`OR_D - M` < 1] <- 'Transiently Downregulated'
  
    data$class[data$`q_value_D - M` < 0.05 & 
             data$`q_value_D - F` < 0.05 & 
             data$`OR_D - M` < 1 & 
             data$`OR_D - F` > 1] <- 'Progressively Downregulated'
    
      data$class[data$`q_value_D - M` < 0.05 & 
                data$`q_value_D - F` < 0.05 & 
             data$`OR_D - M` > 0 & 
             data$`OR_D - F` < 1] <- 'Progressively Upregulated'
      
      data$class[data$`q_value_F - M` < 0.05 & 
                data$`q_value_D - F` > 0.05 & 
                data$`q_value_D - M` > 0.05 & 
             data$`OR_F - M` > 1 ] <- 'Terminally Upregulated'
      
  data$class[data$`q_value_F - M` < 0.05 & 
                data$`q_value_D - F` > 0.05 & 
                data$`q_value_D - M` > 0.05 & 
             data$`OR_F - M` < 1  ] <- 'Terminally Downregulated'

     
      


  
  return(data)
}


}

obj <- readRDS('/Users/ggraham/Desktop/snRNA-seq R Files 122524/RNA Object.rds')

glm_data <- data.frame()
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
  glm_matrix_cluster <- matrix(NA, nrow(f_loop_data_cluster),2)
  glm_matrix_cluster[,1] <- f_loop_data_cluster$cells_in_fish_cluster
  glm_matrix_cluster[,2] <- f_loop_data_cluster$cells_in_fish - f_loop_data_cluster$cells_in_fish_cluster
  
  glm_model <- glm(glm_matrix_cluster~f_loop_data_cluster$Status, family = binomial('logit'), data = f_loop_data_cluster)
  car::Anova(glm_model, type = 'III')
  pairs <- as.data.frame(pairs(emmeans(glm_model, 'Status'),adjust = 'none'))
  
  anova_p.value <- anova(glm_model, test = 'Chisq')[2,5]
  
    d_m_estimate <- pairs$estimate[pairs$contrast=='D - M']
   d_f_estimate <- pairs$estimate[pairs$contrast=='D - F']
    f_m_estimate <- pairs$estimate[pairs$contrast=='F - M']

  
    d_m_p.value <- pairs$p.value[pairs$contrast=='D - M']
   d_f_p.value <- pairs$p.value[pairs$contrast=='D - F']
    f_m_p.value <- pairs$p.value[pairs$contrast=='F - M']
    

    glm_output <- data.frame(cluster = i,
                               d_m_estimate=d_m_estimate,
                               d_f_estimate = d_f_estimate,
                               f_m_estimate = f_m_estimate,
                               anova_p.value = anova_p.value,
                               d_m_p.value = d_m_p.value,
                               d_f_p.value = d_f_p.value,
                               f_m_p.value = f_m_p.value
    )
    
                               
    glm_data <- rbind(glm_data, glm_output)                          
}
glm_data$anova_q.value <- p.adjust(glm_data$anova_p.value, 'fdr',nrow(glm_data))


glm_data$d_m_q.value <- p.adjust(glm_data$d_m_p.value, 'fdr',nrow(glm_data))
glm_data$d_f_q.value <- p.adjust(glm_data$d_f_p.value, 'fdr',nrow(glm_data))
glm_data$f_m_q.value <- p.adjust(glm_data$f_m_p.value, 'fdr',nrow(glm_data))

glm_data_plot <- glm_data %>%
  dplyr::select(cluster,
                d_m_q.value,
                d_f_q.value,
                f_m_q.value,
                d_m_estimate,
                d_f_estimate,
                f_m_estimate) %>%
  pivot_longer(cols = d_m_estimate:f_m_estimate, 
               names_to = "contrast", 
               values_to = "estimate")
glm_data_plot$signif <- NA
glm_data_plot$signif <- ifelse(glm_data_plot$contrast=='d_m_estimat'e&glm_data_plot$d_m_q.value<0.05, '*', glm_data_plot$signif )
glm_data_plot$signif <- ifelse(glm_data_plot$contrast=='d_f_estimate'&glm_data_plot$d_f_q.value<0.05, '*', glm_data_plot$signif )
glm_data_plot$signif <- ifelse(glm_data_plot$contrast=='f_m_estimate'&glm_data_plot$f_m_q.value<0.05, '*', glm_data_plot$signif )

ggplot(glm_data_plot, aes(x = as.factor(cluster), y = estimate, color = contrast, fill = contrast, group = contrast)) +
  geom_bar(stat = 'identity', position = position_dodge(width = 1,preserve = 'single')) +
  geom_text(aes(label = signif), color = 'black', size = 7, position = position_dodge(width = 1))+
  geom_vline(xintercept = c(0:31)+0.5, linewidth = 0.25, linetype=2)

f_loop_data_cluster_2 <- data.frame()
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
    f_loop_data_cluster_2 <- rbind(f_loop_data, f_loop_data_cluster_2)
                              
  }
}

f_loop_data_cluster_summ <- f_loop_data_cluster_2%>%
  subset(Status =='F'|Status=='D'|Status=='M')%>%
  group_by(cluster, Status)%>%
  summarize(prop= sum(cells_in_fish_cluster)/sum(cells_in_fish))
f_loop_data_cluster_summ$Status<- factor(f_loop_data_cluster_summ$Status, levels = c('M','D',"F"))

tip_data <- data.frame(cluster = unique(f_loop_data_cluster_summ$cluster))

tip_data$tips_md <- lapply(tip_data$cluster, function(x) c(x - 0.33, x))
tip_data$tips_mf <- lapply(tip_data$cluster, function(x) c(x - 0.33, x+0.33))
tip_data$tips_fd <- lapply(tip_data$cluster, function(x) c(x, x+0.33))

tip_data$tip_height_md <- f_loop_data_cluster_summ$prop[f_loop_data_cluster_summ$Status=='D']*1.15
tip_data$tip_height_mf <- f_loop_data_cluster_summ$prop[f_loop_data_cluster_summ$Status=='F']*1.15
tip_data$tip_height_fd <- f_loop_data_cluster_summ$prop[f_loop_data_cluster_summ$Status=='M']*1.15

tip_data$signif_md <- glm_data_plot$signif[glm_data_plot$contrast=='d_m_estimate']
tip_data$signif_mf <- glm_data_plot$signif[glm_data_plot$contrast=='f_m_estimate']
tip_data$signif_fd <- glm_data_plot$signif[glm_data_plot$contrast=='d_f_estimate']

f_loop_data_cluster_max <- f_loop_data_cluster_summ%>%
  group_by(cluster)%>%
  summarize(max = max(prop))

tip_data$max <- f_loop_data_cluster_max$max 


library(ggsignif)

ggplot(f_loop_data_cluster_summ, aes(x = as.factor(cluster), y = prop, color = Status, fill = Status, group = Status)) +
  geom_bar(stat = 'identity', position = position_dodge(width = 1,preserve = 'single')) +
  geom_vline(xintercept = c(0:31)+0.5, linewidth = 1, linetype=1, color = 'white')+
    geom_signif(xmin = tip_data$cluster+0.66, xmax = tip_data$cluster+1.33,
                y_position = tip_data$max*1,
                annotation = tip_data$signif_mf , 
                color = "black",
                tip_length = rep(c(0,0),31), textsize = 5) +
    geom_signif(xmin = tip_data$cluster+0.66, xmax = tip_data$cluster+1,
                y_position = tip_data$max*.9,
                annotation = tip_data$signif_md , 
                color = "black",
                tip_length = rep(c(0,0),31), textsize = 5) +
    geom_signif(xmin = tip_data$cluster+1, xmax = tip_data$cluster+1.33,
                y_position = tip_data$max*.9,
                annotation = tip_data$signif_fd , 
                color = "black",
                tip_length = rep(c(0,0),31), textsize = 5) 


  Xgeom_signif(data =tip_data, aes(x = cluster, y = tip_height_mf, annotations =signif_mf), manual = T)



  glm_data_class <- define_degs_glm(glm_data)
glm_data_class[,c(1,ncol(glm_data_class)-1, ncol(glm_data_class))]


