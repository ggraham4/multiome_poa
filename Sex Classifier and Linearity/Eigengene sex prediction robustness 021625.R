#Eigengene sex prediction leave-one-out analysis

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


obj <- readRDS("C:/Users/Gabe/Desktop/RNA Object.rds")

### Read in DEG data ###
degs_list<- list()
for(i in 0:31){
  print(i)
  data <- read.csv(
    paste0('DEG Outputs/012425 Neg Bin w Doms Lower Stringency/cluster_',i,'.csv')
  )
  degs<- data$gene[(data$f_m_q.value<0.15)]
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

#### leave one out analysis

fish_to_assess <- unique(pca_data_pivoted$individual[pca_data_pivoted$Sex %in% c('F','M')])

leave_one_out <- list()
for(c in 0:31){
  print(c)
  if(length(degs_list[[paste0(c)]])<4){next}
  cluster_expression_data <- expression_data[[paste0(c)]]
  cluster_expression_data_pivoted <- cluster_expression_data%>%
    dplyr::select(individual, mean, Sex, gene)%>%
    dplyr::filter( Sex =='F' | Sex == 'M')%>%
    pivot_wider(names_from = gene, 
                values_from = mean)
  new_data_frame = data.frame()
  for(i in fish_to_assess){
    print(i)
    fish_real_sex = cluster_expression_data_pivoted$Sex[cluster_expression_data_pivoted$individual==i]
    if(length(fish_real_sex)==0){next}
    cluster_expression_data_pivoted$Sex[cluster_expression_data_pivoted$individual==i] <- 'Test'
    
    #run PCA
    cluster_prcomp <- prcomp(cluster_expression_data_pivoted[,3:ncol(cluster_expression_data_pivoted)], scale. = T)
    cluster_prcomp_x <- as.data.frame(cluster_prcomp$x)
    cluster_prcomp_x$Sex =cluster_expression_data_pivoted$Sex
    cluster_prcomp_x$individual= cluster_expression_data_pivoted$individual
    cluster_prcomp_x$binary_sex = ifelse(cluster_prcomp_x$Sex == 'F',1, 0)
    cluster_prcomp_x$binary_sex = ifelse(cluster_prcomp_x$individual == i,NA, cluster_prcomp_x$binary_sex)
    
    #train model 1_2 and 1_4 on the known fish
    training_data <- cluster_prcomp_x[!is.na(cluster_prcomp_x$binary_sex),]
    pca_model_1_2 <- glm(binary_sex~PC1+PC2, data=cluster_prcomp_x, family = binomial('logit'))
    pca_model_1_4 <- glm(binary_sex~PC1+PC2+PC3+PC4, data=cluster_prcomp_x, family = binomial('logit'))
    
    data_to_predict = cluster_prcomp_x[cluster_prcomp_x$individual==i,]
    
    predictions_pc1_2_model <- as.data.frame(predict(pca_model_1_2, data_to_predict, type = 'response'))
    predictions_pc1_4_model <- as.data.frame(predict(pca_model_1_4, data_to_predict, type = 'response'))
    
    
    new_data <- data.frame(cluster = c,
                           fish = i,
                           fish_real_sex = fish_real_sex,
                           prediction_1_2 = predictions_pc1_2_model,
                           prediction_1_4 = predictions_pc1_4_model)
    
    colnames(new_data) <- c('cluster','fish','fish_real_sex','prediction_1_2', 'prediction_1_4')
    new_data_frame <- rbind(new_data, new_data_frame)
    
    cluster_expression_data_pivoted$Sex[cluster_expression_data_pivoted$individual==i] <- fish_real_sex
    
  
  }
  leave_one_out[[paste0(c)]] <- new_data_frame
  
}

leave_one_out[['19']]%>%
  pivot_longer(cols = c('prediction_1_2','prediction_1_4'),names_to = 'prediction_type',
               values_to = 'prediction')%>%
  ggplot(aes(x = fish_real_sex, y = prediction, shape = prediction_type))+
  geom_point(position = position_dodge2(0.5))+
  geom_boxplot(aes(fill= prediction_type,alpha = 0.25))

leave_one_out_performance <- data.frame()
for(i in  names(leave_one_out)){
  print(i)
  data = leave_one_out[[i]]
  
  data_pivoted <-   data%>%
    group_by(cluster, fish_real_sex)%>%
    summarize(mean_prediction_1_2 = mean(prediction_1_2),
              se_1_2 = sd(prediction_1_2)/sqrt(n()),
              mean_prediction_1_4 = mean(prediction_1_4),
              se_1_4 = sd(prediction_1_4)/sqrt(n()))
  
  leave_one_out_performance <- rbind(leave_one_out_performance, data_pivoted)
    
  
  
}

leave_one_out_performance_summarized = leave_one_out_performance%>%
  pivot_longer(cols = c('mean_prediction_1_2','mean_prediction_1_4'),
               names_to = 'prediction_type',
               values_to = 'prediction')
leave_one_out_performance_summarized$se_1_2 <- ifelse(leave_one_out_performance_summarized$prediction_type =='mean_prediction_1_2', NA,leave_one_out_performance_summarized$se_1_2)
leave_one_out_performance_summarized$se_1_4 <- ifelse(leave_one_out_performance_summarized$prediction_type =='mean_prediction_1_4', NA,leave_one_out_performance_summarized$se_1_4)
leave_one_out_performance_summarized$se = leave_one_out_performance_summarized$se_1_2
leave_one_out_performance_summarized$se <- ifelse(is.na(leave_one_out_performance_summarized$se),leave_one_out_performance_summarized$se_1_4,leave_one_out_performance_summarized$se)

ggplot(leave_one_out_performance_summarized, aes(x = fct_reorder(as.factor(cluster),mean(prediction)), y = prediction, color =fish_real_sex, shape = prediction_type))+
  geom_pointrange(aes(x = fct_reorder(as.factor(cluster),prediction), y = prediction, ymin = prediction-se, ymax = prediction+se), position = position_dodge2(1))


##### compare to lasso scoring #####
lasso_robust <- data.frame()
for(c in 0:31){
  print(c)
  lasso_exp_data <- expression_data[[paste0(c)]]
  lasso_exp_data_pivoted <- lasso_exp_data%>%
    dplyr::select(individual, mean, Sex, gene)%>%
    dplyr::filter( Sex =='F' | Sex == 'M')%>%
    pivot_wider(names_from = gene, 
                values_from = mean)
  
  
  for(i in fish_to_assess){
    print(i)
    
    fish_real_sex = lasso_exp_data_pivoted$Sex[lasso_exp_data_pivoted$individual==i]
    
    if(length(fish_real_sex)==0){next}
    
    lasso_exp_data_pivoted$Sex[lasso_exp_data_pivoted$individual==i] <- 'Test'
    
    lasso_exp_data_pivoted$fish.class <- ifelse(lasso_exp_data_pivoted$Sex == 'F',1,NA)
    lasso_exp_data_pivoted$fish.class <- ifelse(lasso_exp_data_pivoted$Sex == 'M',0,lasso_exp_data_pivoted$fish.class)
    
    training_data <- lasso_exp_data_pivoted[!is.na(lasso_exp_data_pivoted$fish.class),]
    
    
    training_data_model <- training_data[, c(4:ncol(training_data)-1)]
    
    if(ncol(training_data)<4){next}
    
    x.training <- as.matrix(training_data_model)
    
    if(ncol(x.training)<2){return(next)}
    
    y.training <- training_data$fish.class
    
    #calculate lambda
    lasso <-cv.glmnet(y = y.training, x = x.training, family = "binomial", alpha = 1, lambda = NULL)
    
    #use as lambda in final trainer
    min <- lasso$lambda.min
    
    #train
    lasso.final <- glmnet(x=x.training, y=y.training, alpha = 1, family = "binomial",
                          lambda = min)
    
    
    
    test_data <- subset(lasso_exp_data_pivoted, is.na(fish.class))
    
    x.test <- as.matrix(na.omit(test_data[,(4:ncol(test_data)-1)]))
    if(ncol(x.test)<2){return(NULL)}
    
    probabilities <- lasso.final %>% predict(newx = x.test, type = 'response')
    
    predicted_data <- as.data.frame(probabilities)
    
    predicted_data$individual = i
    predicted_data$fish_real_sex = fish_real_sex
    predicted_data$cluster = c
    
    lasso_robust <- rbind(predicted_data, lasso_robust)
    
    #reset the sex
    lasso_exp_data_pivoted$Sex[lasso_exp_data_pivoted$individual==i] <- fish_real_sex
    
    }
  
  
}
colnames(lasso_robust) <- c('s0','individual','fish_real_sex','cluster')

deg_data <- data.frame()
for(i in 0:31){
  n_degs = length(degs_list[[paste0(i)]])
  
  deg_data2 <- data.frame(cluster = i,
                         n_degs=n_degs)
  deg_data<- rbind(deg_data, deg_data2)
  
  }

lasso_robust%>%
  group_by(cluster, fish_real_sex)%>%
  summarize(mean_pred = mean(s0),
            se = sd(s0)/sqrt(n()))%>%
  right_join(deg_data, by = 'cluster')%>%
  ggplot(aes(x = fct_reorder(as.factor(cluster), mean_pred), y = mean_pred, color = fish_real_sex))+
  geom_pointrange(aes(x = fct_reorder(as.factor(cluster),mean_pred),
                      y = mean_pred,
                      ymin = mean_pred-se,
                      ymax = mean_pred+se,
                      size = n_degs), position = position_dodge2(0.5))+
  scale_size_continuous(range = c(0,1.5))


