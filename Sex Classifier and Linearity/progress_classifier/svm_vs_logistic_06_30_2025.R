# eigengene sex prediction with new data

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

obj <- readRDS("~/Desktop/optimal_clustering_rna_only.rds")


### Read in DEG data ###
for(i in 0:25){
  print(i)
  if(i == 21){next}
  data <- read.csv(
    paste0('~/Desktop/multiome_poa/Sex Classifier and Linearity/progress_classifier/degs_0.1_f_m/cluster_',i,'.csv'))
  assign(paste0('cluster_',i,'_degs'), data)
}

#run principal component analysis on DEGs ####
pca_output <- list()
for(i in 0:25){
  print(i)
  if(i == 21){next}
  #extract data frame
  pca_data <- get(paste0('cluster_',i,'_degs'))
  
  #subset
  pca_data_pivoted <- pca_data %>%
    dplyr::filter(Status == 'D'| Status =='F' | Status == 'M')
  
  #run PCA on individuals
  cluster_prcomp <- prcomp(pca_data_pivoted[,4:ncol(pca_data_pivoted)], scale. = T)
  
  pca_output[[paste0(i)]]$prcomp = cluster_prcomp
}

### Train models and test individuals ####
pca_predictions <- list()
pca_model_output <- list()
for(i in 0:25){
  print(i)
    if(i == 21){next}

  #extract data frame
  pca_data <- get(paste0('cluster_',i,'_degs'))
  
  #subset
  pca_data_pivoted <- pca_data %>%
    dplyr::filter(Status == 'D'| Status =='F' | Status == 'M')
  
  
  
  data <- as.data.frame(pca_output[[paste0(i)]]$prcomp$x)
  if(ncol(data)<2){next}
  
  data$Status <- pca_data_pivoted$Status
  data$individual <- pca_data_pivoted$individual
  
  
  data$binary_sex <- ifelse(data$Status == 'F', 1,NA)
  data$binary_sex <- ifelse(data$Status == 'M', 0, data$binary_sex)
  
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
  summarize(mean_pred = mean(predicted),
            se_pred = sd(predicted)/sqrt(n()))


### read in SVM data
svm_data = read.csv('/Users/ggraham/Desktop/multiome_poa/Sex Classifier and Linearity/progress_classifier/svm_predictions_dominants.csv')


joined = pca_predictions_merged_cluster%>%
  right_join(svm_data, join_by('cluster'))

ggplot(joined, aes(x = as.factor(cluster), y = mean_pred, color = 'Logistic PC1 + PC2'))+
  geom_point(aes(shape ='Logistic PC1 + PC2'), shape = 17, size = 4)+
  geom_point(aes(y = X0, color = 'SVM', shape = 'SVM'), shape = 16, size = 4)+
  labs(x = 'Prediction', y = 'Cluster')


ggplot(joined, aes(x = X0, y = mean_pred), size = 4)+
  geom_point()+
  geom_smooth(method = 'lm')+
  labs(x = 'SVM', y = 'Logistic', title = 'Prediction Correlations', subtitle = 'Pearsons r = 0.298')

  cor.test(joined$X0, joined$mean_pred)
  
  
  
  
  
