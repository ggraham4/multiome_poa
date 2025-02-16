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
mean_expression_cluster_data<- readRDS('Functions/mean_expression_cluster_data.rds')


### Read in DEG data ###
degs_list<- list()
for(i in 0:31){
  print(i)
  data <- read.csv(
    paste0('DEG Outputs/012425 Neg Bin w Doms Lower Stringency/cluster_',i,'.csv')
    )
  degs<- data$gene[(data$f_m_q.value<0.1|
                   data$d_m_q.value<0.1|
                   data$d_f_q.value<0.1)]
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
  cluster_prcomp <- prcomp(pca_data_pivoted[,3:ncol(pca_data_pivoted)])
  
  pca_output[[paste0(i)]] = cluster_prcomp
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
  
  
  data <- as.data.frame(pca_output[[paste0(i)]]$x)
  data$Sex <- pca_data_pivoted$Sex
  data$individual <- pca_data_pivoted$individual
  
  
  data$binary_sex <- ifelse(data$Sex == 'F', 1,NA)
  data$binary_sex <- ifelse(data$Sex == 'M', 0, data$binary_sex)
  
  training_data <- data[!is.na(data$binary_sex),]
  
  ###train model
  pc1_model <- glm(binary_sex ~ PC1, data = training_data, family = binomial('logit'))
  
  pc2_model <- glm(binary_sex ~ PC2, data = training_data, family = binomial('logit'))
  
  pc1_2_model <- glm(binary_sex ~ PC1+PC2, data = training_data, family = binomial('logit'))
  
  ### collect p values
  pc1_2_model_av = anova(pc1_2_model, test ='Chisq')
  
  pc1_p.value = summary(pc1_model)$coefficients[1,4]
  pc2_p.value = summary(pc2_model)$coefficients[1,4]
  
  pc_1_2_pc1_p.value = pc1_2_model_av$`Pr(>Chi)`[2]
  pc_1_2_pc2_p.value = pc1_2_model_av$`Pr(>Chi)`[3]
  
  cluster_outputs <- data.frame(pc1_p.value = pc1_p.value,
                                pc2_p.value = pc2_p.value,
                                pc_1_2_pc1_p.value = pc_1_2_pc1_p.value,
                                pc_1_2_pc1_p.value= pc_1_2_pc2_p.value)
  pca_model_output[[paste0(i)]]=cluster_outputs
  
  #### test model ####
  test_data <- data[is.na(data$binary_sex),]
  
  #pc1
  predictions_pc1_model <- as.data.frame(predict(pc1_model, test_data, type = 'response'))
  predictions_pc1_model$individual = test_data$individual
  colnames(predictions_pc1_model) <- c('predicted','individual')
  predictions_pc1_model$model = 'pc1'
  
  #pc2
  predictions_pc2_model <- as.data.frame(predict(pc2_model, test_data, type = 'response'))
  predictions_pc2_model$individual = test_data$individual
  colnames(predictions_pc2_model) <- c('predicted','individual')
  predictions_pc2_model$model = 'pc2'
  
  #pc1_2
  predictions_pc1_2_model <- as.data.frame(predict(pc1_2_model, test_data, type = 'response'))
  predictions_pc1_2_model$individual = test_data$individual
  colnames(predictions_pc1_2_model) <- c('predicted','individual')
  predictions_pc1_2_model$model = 'pc1_2'
  
  predictions = rbind(predictions_pc1_model, predictions_pc2_model, predictions_pc1_2_model)
  predictions$cluster =i
  pca_predictions[[paste0(i)]] = predictions
  
}

pca_predictions_merged = do.call(rbind, pca_predictions)

###pc1
pca_predictions_pc1 <- pca_predictions_merged%>%
  subset(model == 'pc1')%>%
  group_by(cluster)%>%
  summarize(mean_pred = mean(predicted),
            se_pred = sd(predicted)/sqrt(n()))

library(forcats)
pca_predictions_pc1$guess <- ifelse(pca_predictions_pc1$mean_pred>0.66, 'Female',NA)
pca_predictions_pc1$guess <- ifelse(pca_predictions_pc1$mean_pred<0.33, 'Male',pca_predictions_pc1$guess)
pca_predictions_pc1$guess <- ifelse(is.na(pca_predictions_pc1$guess), 'Intermediate',pca_predictions_pc1$guess)


ggplot(pca_predictions_pc1, aes(x = fct_reorder(as.factor(cluster), mean_pred), y = mean_pred, color = guess))+
  geom_pointrange(aes(x = fct_reorder(as.factor(cluster), mean_pred), y = mean_pred, ymin = mean_pred-se_pred, ymax = mean_pred+se_pred))

### pc2

pca_predictions_pc2 <- pca_predictions_merged%>%
  subset(model == 'pc2')%>%
  group_by(cluster)%>%
  summarize(mean_pred = mean(predicted),
            se_pred = sd(predicted)/sqrt(n()))

library(forcats)
pca_predictions_pc2$guess <- ifelse(pca_predictions_pc2$mean_pred>0.66, 'Female',NA)
pca_predictions_pc2$guess <- ifelse(pca_predictions_pc2$mean_pred<0.33, 'Male',pca_predictions_pc2$guess)
pca_predictions_pc2$guess <- ifelse(is.na(pca_predictions_pc2$guess), 'Intermediate',pca_predictions_pc2$guess)


ggplot(pca_predictions_pc2, aes(x = fct_reorder(as.factor(cluster), mean_pred), y = mean_pred, color = guess))+
  geom_pointrange(aes(x = fct_reorder(as.factor(cluster), mean_pred), y = mean_pred, ymin = mean_pred-se_pred, ymax = mean_pred+se_pred))

#pc1_2

pca_predictions_pc1_2 <- pca_predictions_merged%>%
  subset(model == 'pc1_2')%>%
  group_by(cluster)%>%
  summarize(mean_pred = mean(predicted),
            se_pred = sd(predicted)/sqrt(n()))

library(forcats)
pca_predictions_pc1_2$guess <- ifelse(pca_predictions_pc1_2$mean_pred>0.66, 'Female',NA)
pca_predictions_pc1_2$guess <- ifelse(pca_predictions_pc1_2$mean_pred<0.33, 'Male',pca_predictions_pc1_2$guess)
pca_predictions_pc1_2$guess <- ifelse(is.na(pca_predictions_pc1_2$guess), 'Intermediate',pca_predictions_pc1_2$guess)


ggplot(pca_predictions_pc1_2, aes(x = fct_reorder(as.factor(cluster), mean_pred), y = mean_pred, color = guess))+
  geom_pointrange(aes(x = fct_reorder(as.factor(cluster), mean_pred), y = mean_pred, ymin = mean_pred-se_pred, ymax = mean_pred+se_pred))

#### how well are they correlated

pca_predictions_together = pca_predictions_pc1%>%
  full_join(pca_predictions_pc2, by = 'cluster')%>%
  full_join(pca_predictions_pc1_2, by = 'cluster')

ggplot(pca_predictions_together, aes(x = mean_pred.x, y  = mean_pred.y))+
  geom_smooth(method = 'lm')

cor.test(pca_predictions_together$mean_pred.x, pca_predictions_together$mean_pred.y)
#-0.37

ggplot(pca_predictions_together, aes(x = mean_pred.x, y  = mean_pred))+
  geom_smooth(method = 'lm')
cor.test(pca_predictions_together$mean_pred.x, pca_predictions_together$mean_pred)
#0.47

ggplot(pca_predictions_together, aes(x = mean_pred.y, y  = mean_pred))+
  geom_smooth(method = 'lm')
cor.test(pca_predictions_together$mean_pred, pca_predictions_together$mean_pred.y)
#0.33


