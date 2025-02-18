# Eigengene sex prediction final
# here, I run PCA on all individuals, train a logistic on pcs 1 and 2 to predict male (0) or female (1), 
# and then use that trained model to predict the sex of dominants
#I then summarize across clusters and individuals

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
}

obj <- readRDS("C:/Users/Gabe/Desktop/RNA Object.rds")

mean_expression_cluster_data<- function(object, gene, cluster, clustering = 'harmony.wnn_res0.4_clusters'){
  options(dplyr.summarise.inform = FALSE)
  
  Counts_of_interest <- as.data.frame(counts[,gene])
  Counts_of_interest$expression <- Counts_of_interest[,1]
  Counts_of_interest$individual <- object@meta.data$individual[object@meta.data[[clustering]] == cluster]
  Counts_of_interest$Status <- object@meta.data$Status[object@meta.data[[clustering]] == cluster]
  
  results <-Counts_of_interest%>%
    group_by(individual, Status)%>%
    summarize(mean = mean(expression),
              se = sd(expression)/sqrt(n()))
  results$Sex <- results$Status
  
  
  results$Sex <- str_sub(results$individual, -1)
  results$Sex[results$individual == 'T17D'] = 'NF'
  results$Sex[results$individual == 'A12D'] = 'E'
  results$Sex[results$individual == 'T11D'] = 'E'
  results$Sex[results$individual == 'GH'] = 'NRM'
  return(results)
}


### Read in DEG data ###
degs_list<- list()
for(i in 0:31){
  print(i)
  data <- read.csv(
    paste0('DEG Outputs/012425 Neg Bin w Doms Lower Stringency/cluster_',i,'.csv')
  )
  degs<- data$gene[(data$f_m_q.value<0.1)]
  degs <- degs[!is.na(degs)]
  degs <- c(degs)
  
  degs_list[[paste0(i)]] = c(degs)
  
}

### Summarize Gene expression of DEGs ### 
expression_data <- list()
for(i in 0:31){
  print(i)
  cluster = paste0(i)
  cluster_data <- data.frame()
  
  counts <- t(obj@assays$RNA$data[,obj@meta.data$harmony.wnn_res0.4_clusters == i]) 
  #remove counts from mean expression function to make it faster
  
  for(g in degs_list[[cluster]]){
    gene_expression_data <- mean_expression_cluster_data(object = obj,
                                                         gene = g,
                                                         cluster = i)
    gene_expression_data$gene = g
    cluster_data <- rbind(cluster_data, gene_expression_data)
    
  }
  expression_data[[paste0(i)]] = cluster_data
}

#run principal component analysis on DEGs ####
pca_output <- list()
for(i in 0:31){
  print(i)
  
  #extract data frame
  pca_data <- expression_data[[paste0(i)]]
  
  #pivot and subset
  pca_data_pivoted <- pca_data %>%
    dplyr::select(individual, mean, Sex, gene)%>%
    dplyr::filter(Sex == 'D'| Sex =='F' | Sex == 'M')%>%
    pivot_wider(names_from = gene, 
                values_from = mean)
  
  #run PCA on individuals
  cluster_prcomp <- prcomp(pca_data_pivoted[,3:ncol(pca_data_pivoted)], scale. = T)
  
  pca_output[[paste0(i)]]$prcomp = cluster_prcomp
}

### Train models and test individuals ####
pca_predictions <- list()
pca_model_output <- list()
for(i in 0:31){
  print(i)
  
  #make training data
  pca_data <- expression_data[[paste0(i)]]
  
  pca_data_pivoted <- pca_data %>%
    dplyr::select(individual, mean, Sex, gene)%>%
    dplyr::filter(Sex == 'D'| Sex =='F' | Sex == 'M')%>%
    pivot_wider(names_from = gene, 
                values_from = mean)
  
  
  data <- as.data.frame(pca_output[[paste0(i)]]$prcomp$x)
  if(ncol(data)<=2){next}
  
  data$Sex <- pca_data_pivoted$Sex
  data$individual <- pca_data_pivoted$individual
  
  
  data$binary_sex <- ifelse(data$Sex == 'F', 1,NA)
  data$binary_sex <- ifelse(data$Sex == 'M', 0, data$binary_sex)
  
  training_data <- data[!is.na(data$binary_sex),]
  
  ###train model
  pc1_2_model <- glm(binary_sex ~ PC1+PC2, data = training_data, family = binomial('logit'))
  
  pc1_2_model_av = anova(pc1_2_model, test ='Chisq')
  
  pca_model_output[[paste0(i)]]$pc1_2_model = pc1_2_model
  pca_model_output[[paste0(i)]]$pc1_2_av = pc1_2_model_av
  
  
  #### test model ####
  test_data <- data[is.na(data$binary_sex),]
  
  
  #pc1_2
  predictions_pc1_2_model <- as.data.frame(predict(pc1_2_model, test_data, type = 'response'))
  predictions_pc1_2_model$individual = test_data$individual
  colnames(predictions_pc1_2_model) <- c('predicted','individual')
  predictions_pc1_2_model$model = 'pc1_2'
  
  
  predictions = rbind(predictions_pc1_2_model)
  predictions$cluster =i
  predictions$n_degs = ncol(pca_data_pivoted)-2
  pca_predictions[[paste0(i)]] = predictions
  
}

pca_predictions_merged = do.call(rbind, pca_predictions)

library(forcats)
`%notin%` <- Negate(`%in%`)

pca_predictions_merged_cluster <- pca_predictions_merged%>% 
  group_by(cluster)%>%
  filter(cluster %notin% c(15,30))%>%
  summarize(mean_pred = mean(predicted),
            se_pred = sd(predicted)/sqrt(n()))

p_val_data <- data.frame()
for(i in 0:31){
  model <- pca_model_output[[paste0(i)]][["pc1_2_av"]]
  if(is.null(model)){next}
  model_summary <- as.data.frame(model)
  
  pc1_p.val <- model_summary$`Pr(>Chi)`[2]
  pc2_p.val <- model_summary$`Pr(>Chi)`[3]
  
  new_data <- data.frame(cluster = i,
                         pc1_p.val = pc1_p.val,
                         pc2_p.val = pc2_p.val)
  p_val_data <- rbind(p_val_data, new_data)
  
  }
p_val_data$asterisk <- ifelse((p_val_data$pc1_p.val)<0.05 |(p_val_data$pc2_p.val)<0.05, '*', NA)



pca_predictions_merged_cluster$Progression <- ifelse(pca_predictions_merged_cluster$mean_pred<0.33,'Male',NA)
pca_predictions_merged_cluster$Progression <- ifelse(pca_predictions_merged_cluster$mean_pred>0.66,'Female',pca_predictions_merged_cluster$Progression)
pca_predictions_merged_cluster$Progression <- ifelse(is.na(pca_predictions_merged_cluster$Progression), 'Intermediate',pca_predictions_merged_cluster$Progression) 

pca_predictions_merged_cluster_plot <- pca_predictions_merged_cluster%>%
  right_join(p_val_data, by = 'cluster')

pca_pred_plot <- ggplot(subset(pca_predictions_merged_cluster_plot, cluster %notin% c(15,30)), aes(x = fct_reorder(as.factor(cluster), mean_pred), y = mean_pred, color = Progression, shape = Progression))+
  geom_pointrange(aes(x = fct_reorder(as.factor(cluster), mean_pred), y = mean_pred, ymin = mean_pred-se_pred, ymax=mean_pred+se_pred))+
  geom_text(aes(x =  fct_reorder(as.factor(cluster), mean_pred), y = 1.05*(mean_pred+se_pred+0.01), label = asterisk), color = 'black', size =10)+
  labs(x = 'Cluster', y = 'Mean Prediction +/- SE')+
  theme_classic()+
  theme(legend.position = c(0.15,.85))+
  theme(axis.text.x = element_text(size = 10,angle = -45, vjust = 1, hjust=0), axis.text.y = element_text(size = 10), axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12))

pca_pred_plot
length(pca_predictions_merged_cluster_plot$cluster[pca_predictions_merged_cluster_plot$cluster %notin% c(15,30)])
ggsave(plot = pca_pred_plot,
       file = "pca_pred_plot.svg",
       device = "svg",
       units = "in",
       width = 5,
       height = 5,
       path = "Bachelors Thesis/Plots/Figure 2")

