#Dr rhodes wants me to do PCs 1 - 4 and also doesnt believe pc1 isnt a significant predictor of sex most of the time
#so i need to prove that via plotting
#I want to save in my lists 
#* Models
#* Plots
#* scores
#* maybe even leave one out
#* only look at genes that differ than male and female
#* ignore clusters with only 1 deg, add a size metric when plotting clusters in the dotplot

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
  data$Sex <- pca_data_pivoted$Sex
  data$individual <- pca_data_pivoted$individual
  
  if(ncol(data)<=5){next}
  
  data$binary_sex <- ifelse(data$Sex == 'F', 1,NA)
  data$binary_sex <- ifelse(data$Sex == 'M', 0, data$binary_sex)
  
  training_data <- data[!is.na(data$binary_sex),]
  
  ###train model
  pc1_model <- glm(binary_sex ~ PC1, data = training_data, family = binomial('logit'))
  
  pc2_model <- glm(binary_sex ~ PC2, data = training_data, family = binomial('logit'))
  
  pc1_2_model <- glm(binary_sex ~ PC1+PC2, data = training_data, family = binomial('logit'))
  
  pc1_4_model <- glm(binary_sex ~ PC1+PC2+PC3+PC4, data = training_data, family = binomial('logit'))
  
  pc1_2_model_av = anova(pc1_2_model, test ='Chisq')
  pc1_4_model_av= anova(pc1_4_model, test ='Chisq')
  
  pca_model_output[[paste0(i)]]$pc1_model = pc1_model
  pca_model_output[[paste0(i)]]$pc1_model = pc1_model
  pca_model_output[[paste0(i)]]$pc1_2_model = pc1_2_model
  pca_model_output[[paste0(i)]]$pc1_2_av = pc1_2_model_av
  pca_model_output[[paste0(i)]]$pc1_4_model = pc1_4_model
  pca_model_output[[paste0(i)]]$pc1_4_av = pc1_4_model_av
  
  pca_model_output[[paste0(i)]]$pc1_4_anova = pc1_4_model_av
  
  
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
  
  #pc1_4
  predictions_pc1_4_model <- as.data.frame(predict(pc1_4_model, test_data, type = 'response'))
  predictions_pc1_4_model$individual = test_data$individual
  colnames(predictions_pc1_4_model) <- c('predicted','individual')
  predictions_pc1_4_model$model = 'pc1_4'
  
  
  predictions = rbind(predictions_pc1_model, predictions_pc2_model, predictions_pc1_2_model, predictions_pc1_4_model)
  predictions$cluster =i
  predictions$n_degs = ncol(pca_data_pivoted)-2
  pca_predictions[[paste0(i)]] = predictions
  
}

pca_predictions_merged = do.call(rbind, pca_predictions)

library(forcats)
###pc1
pca_predictions_pc1 <- pca_predictions_merged%>%
  subset(model == 'pc1')%>%
  group_by(cluster, n_degs)%>%
  summarize(mean_pred = mean(predicted),
            se_pred = sd(predicted)/sqrt(n()))

pca_predictions_pc1$guess <- ifelse(pca_predictions_pc1$mean_pred>0.66, 'Female',NA)
pca_predictions_pc1$guess <- ifelse(pca_predictions_pc1$mean_pred<0.33, 'Male',pca_predictions_pc1$guess)
pca_predictions_pc1$guess <- ifelse(is.na(pca_predictions_pc1$guess), 'Intermediate',pca_predictions_pc1$guess)


ggplot(pca_predictions_pc1, aes(x = fct_reorder(as.factor(cluster), mean_pred), y = mean_pred, color = guess))+
  geom_pointrange(aes(x = fct_reorder(as.factor(cluster), mean_pred), y = mean_pred, ymin = mean_pred-se_pred, ymax = mean_pred+se_pred, size = n_degs))+
  scale_size_continuous(range = c(0,2))

### pc2
pca_predictions_pc2 <- pca_predictions_merged%>%
  subset(model == 'pc2')%>%
  group_by(cluster, n_degs)%>%
  summarize(mean_pred = mean(predicted),
            se_pred = sd(predicted)/sqrt(n()))

pca_predictions_pc2$guess <- ifelse(pca_predictions_pc2$mean_pred>0.66, 'Female',NA)
pca_predictions_pc2$guess <- ifelse(pca_predictions_pc2$mean_pred<0.33, 'Male',pca_predictions_pc2$guess)
pca_predictions_pc2$guess <- ifelse(is.na(pca_predictions_pc2$guess), 'Intermediate',pca_predictions_pc2$guess)


ggplot(pca_predictions_pc2, aes(x = fct_reorder(as.factor(cluster), mean_pred), y = mean_pred, color = guess))+
  geom_pointrange(aes(x = fct_reorder(as.factor(cluster), mean_pred), y = mean_pred, ymin = mean_pred-se_pred, ymax = mean_pred+se_pred, size = n_degs))+
  scale_size_continuous(range = c(0,2))

#pc1_2
pca_predictions_pc1_2 <- pca_predictions_merged%>%
  subset(model == 'pc1_2')%>%
  group_by(cluster, n_degs)%>%
  summarize(mean_pred = mean(predicted),
            se_pred = sd(predicted)/sqrt(n()))

pca_predictions_pc1_2$guess <- ifelse(pca_predictions_pc1_2$mean_pred>0.66, 'Female',NA)
pca_predictions_pc1_2$guess <- ifelse(pca_predictions_pc1_2$mean_pred<0.33, 'Male',pca_predictions_pc1_2$guess)
pca_predictions_pc1_2$guess <- ifelse(is.na(pca_predictions_pc1_2$guess), 'Intermediate',pca_predictions_pc1_2$guess)


ggplot(pca_predictions_pc1_2, aes(x = fct_reorder(as.factor(cluster), mean_pred), y = mean_pred, color = guess))+
  geom_pointrange(aes(x = fct_reorder(as.factor(cluster), mean_pred), y = mean_pred, ymin = mean_pred-se_pred, ymax = mean_pred+se_pred, size = n_degs))+
  scale_size_continuous(range = c(0,2))

#pc1_4
pca_predictions_pc1_4 <- pca_predictions_merged%>%
  subset(model == 'pc1_4')%>%
  group_by(cluster, n_degs)%>%
  summarize(mean_pred = mean(predicted),
            se_pred = sd(predicted)/sqrt(n()))

pca_predictions_pc1_4$guess <- ifelse(pca_predictions_pc1_4$mean_pred>0.66, 'Female',NA)
pca_predictions_pc1_4$guess <- ifelse(pca_predictions_pc1_4$mean_pred<0.33, 'Male',pca_predictions_pc1_4$guess)
pca_predictions_pc1_4$guess <- ifelse(is.na(pca_predictions_pc1_4$guess), 'Intermediate',pca_predictions_pc1_4$guess)


ggplot(pca_predictions_pc1_4, aes(x = fct_reorder(as.factor(cluster), mean_pred), y = mean_pred, color = guess))+
  geom_pointrange(aes(x = fct_reorder(as.factor(cluster), mean_pred), y = mean_pred, ymin = mean_pred-se_pred, ymax = mean_pred+se_pred, size = n_degs))+
  scale_size_continuous(range = c(0,2))

#### how well are they correlated

pca_predictions_together = pca_predictions_pc1
pca_predictions_together$pred_1 =pca_predictions_pc1$mean_pred
pca_predictions_together$pred_2 =pca_predictions_pc2$mean_pred
pca_predictions_together$pred1_2 =pca_predictions_pc1_2$mean_pred
pca_predictions_together$pred1_4 =pca_predictions_pc1_4$mean_pred

#1, 1_2
attach(pca_predictions_together)
cor.test(pred_1, pred1_2)
# 0.85

cor.test(pred_1, pred1_4)
#0.73

cor.test(pred1_2, pred1_4)
#0.74
detach(pca_predictions_together)


### add asterisks to model
pca_predictions_together$pc1_asterisk <- NA
pca_predictions_together$pc1_2_asterisk <- NA
pca_predictions_together$pc1_4_asterisk <- NA

for (i in unique(pca_predictions_together$cluster)) {  
  print(i)
  
  pc1_model <- pca_model_output[[paste0(i)]]$pc1_model
  p_val_1 <- summary(pc1_model)$coefficients[2,4]
  pca_predictions_together$pc1_asterisk[which(pca_predictions_together$cluster == i)] <- ifelse(p_val_1 < 0.05, '*', NA)
  
  pc1_2_av <- as.data.frame(pca_model_output[[paste0(i)]]$pc1_2_av)[-1,]
  
  pca_predictions_together$pc1_2_asterisk[which(pca_predictions_together$cluster == i)] <- ifelse(any(pc1_2_av$`Pr(>Chi)` < 0.05), '*', NA)

  
  pc1_4_av <- as.data.frame(pca_model_output[[paste0(i)]]$pc1_4_av)[-1,]
  
  pca_predictions_together$pc1_4_asterisk[which(pca_predictions_together$cluster == i)] <- ifelse(any(pc1_4_av$`Pr(>Chi)`< 0.05) , '*', NA)
  
}


### how do these compare to the lasso scores####
mean_expression_cluster_data_2 <- function(object, gene, cluster, clustering = 'harmony.wnn_res0.4_clusters'){
  counts <- t(obj@assays$RNA$data[,obj@meta.data$harmony.wnn_res0.4_clusters == cluster])
  Counts_of_interest <- as.data.frame(counts[,gene])
  Counts_of_interest[[gene]] <- Counts_of_interest[,1]
  Counts_of_interest$individual <- object@meta.data$individual[object@meta.data[[clustering]] == cluster]
  results <- data.frame()
  
  for (i in unique(object@meta.data$individual)) {
    Counts <- Counts_of_interest[[gene]][Counts_of_interest$individual==i]
    mean <- mean(Counts)
    mean_se <- sd(Counts) / sqrt(length(Counts))
    df <- data.frame(
      individual = i,
      mean = mean,
      se = mean_se
    )
    results <- rbind(results, df)
  }
  results$Sex <- str_sub(results$individual, -1)
  results$Sex[results$individual == 'T17D'] = 'NF'
  results$Sex[results$individual == 'A12D'] = 'E'
  results$Sex[results$individual == 'T11D'] = 'E'
  results$Sex[results$individual == 'GH'] = 'NRM'
  return(results)
}


lasso_scorer <- function(deg_data, cluster){
  #subset data
  data<- subset(deg_data, cluster == cluster & !is.na(gene))
  #remove clusters with no DEGs
  if(nrow(data)<1){return(NULL)}
  #coerce to numeric
  data$d_f_q.value <- as.numeric(data$d_f_q.value)
  data$d_m_q.value <- as.numeric(data$d_m_q.value)
  data$f_m_q.value <- as.numeric(data$f_m_q.value)
  
  #list degs
  
  degs <- data$gene[data$f_m_q.value<0.15| #chagning these to 0.15 so the lasso has access to the same genes the pca does
                      data$d_f_q.value<0.15|
                      data$d_m_q.value<0.15]
  degs <- degs[!is.na(degs)]
  
  if(length(degs) ==0){return(NULL)}
  
  degs_expression <- data.frame()
  for(b in degs){
    gene_expression <-mean_expression_cluster_data_2(
      obj,
      b,
      cluster
    )
    
    new_data <-data.frame(
      cluster= cluster,
      gene = b,
      mean = gene_expression$mean,
      sex = gene_expression$Sex,
      individual = gene_expression$individual
    )
    degs_expression <- rbind(new_data, degs_expression)
  }
  
  #Define sexes as binomial  
  degs_expression$fish.class <- ifelse(degs_expression$sex == 'F',1,NA)
  degs_expression$fish.class <- ifelse(degs_expression$sex == 'M',0,degs_expression$fish.class)
  #remove the other weird sexes
  degs_expression <- subset(degs_expression, sex == 'F'|
                              sex=='M'|
                              sex == 'D')
  
  #pivot to make matrix
  pivoted_data<- degs_expression%>%
    pivot_wider(names_from = gene, 
                values_from = mean)
  
  #training should only be males and females
  training_data <- pivoted_data[!is.na(pivoted_data$fish.class),]
  training_data <- training_data[complete.cases(training_data[, 5:ncol(training_data)]), ]
  
  
  
  x.training <- as.matrix(na.omit(training_data[,5:ncol(training_data)]))
  if(ncol(x.training)<2){return(NULL)}
  
  y.training <- training_data$fish.class
  
  #calculate lambda
  lasso <-cv.glmnet(y = y.training, x = x.training, family = "binomial", alpha = 1, lambda = NULL)
  
  #use as lambda in final trainer
  min <- lasso$lambda.min
  
  #train
  lasso.final <- glmnet(x=x.training, y=y.training, alpha = 1, family = "binomial",
                        lambda = min)
  
  #now predict dominants
  test_data <- pivoted_data[is.na(pivoted_data$fish.class),]
  test_data <- test_data[complete.cases(test_data[, 5:ncol(test_data)]), ]
  
  
  x.test <- as.matrix(na.omit(test_data[,5:ncol(test_data)]))
  if(ncol(x.test)<2){return(NULL)}
  
  
  #predict probabilities
  probabilities <- lasso.final %>% predict(newx = x.test, type = 'response')
  
  
  #make a dataframe with results
  predicted.classes <- ifelse(probabilities > 0.5, 'f', "m")
  
  predicted_data <- as.data.frame(probabilities)
  made.data <- as.data.frame(probabilities)
  made.data$predicted <- predicted.classes
  made.data$fish <- test_data$individual
  made.data$probabilities <- probabilities
  
  return(made.data)
  
}

probability_data <- data.frame()
for(i in 0:31){
  print(i)
  file.name <-  paste0('C:/Users/Gabe/Desktop/multiome_poa/DEG Outputs/012425 Neg Bin w Doms Lower Stringency/cluster_',
                       i,'.csv')
  
  data<- read.csv(file.name)
 
   prediction <-#tryCatch({
    lasso_scorer(data, i)
  #)}, error=function(e){
  #message('error')
  #return(NULL)})
  if (is.null(prediction)) next
  prediction$cluster = i
  probability_data <- rbind(probability_data,prediction)
  
}
data_for_plot <- probability_data%>%
  group_by(cluster)%>%
  summarize(mean_pred = mean(s0),
            se = sd(s0)/sqrt(n())
  )
data_for_plot$sex = ifelse(data_for_plot$mean_pred>0.66, 'Female', NA)
data_for_plot$sex = ifelse(data_for_plot$mean_pred<0.33, 'Male', data_for_plot$sex)
data_for_plot$sex = ifelse(is.na(data_for_plot$sex ), 'intermediate', data_for_plot$sex )

ggplot(data_for_plot, aes(x = fct_reorder(as.factor(cluster), mean_pred), y = mean_pred, color = sex))+
  geom_pointrange(aes(x = fct_reorder(as.factor(cluster), mean_pred), y = mean_pred, ymin = mean_pred - se, ymax= mean_pred+se))


### how do they correlate
pca_predictions_together_lasso <- pca_predictions_together%>%
  right_join(data_for_plot, by = 'cluster')

cor.test(pca_predictions_together_lasso$pred_1, pca_predictions_together_lasso$mean_pred.y)
#0.28

cor.test(pca_predictions_together_lasso$pred_2, pca_predictions_together_lasso$mean_pred.y)
#0.16

cor.test(pca_predictions_together_lasso$pred1_2, pca_predictions_together_lasso$mean_pred.y)
#0.322

cor.test(pca_predictions_together_lasso$pred1_4, pca_predictions_together_lasso$mean_pred.y)
#0.47

#none of them correlate very well

ggplot(pca_predictions_together_lasso, aes(x = fct_reorder(as.factor(cluster), mean_pred.y), y = mean_pred.y, color = sex))+
  geom_pointrange(data=pca_predictions_together_lasso, aes(x = fct_reorder(as.factor(cluster), mean_pred.y),
                      y = mean_pred.y,
                      ymin = mean_pred.y - se, 
                      ymax= mean_pred.y+se))



