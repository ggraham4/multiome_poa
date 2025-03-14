#Redoing lasso predictions, continuity predictions with threshold 0.1
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
  
  mean_expression_cluster_data <- function(object, gene, cluster, clustering = 'harmony.wnn_res0.4_clusters'){
    counts <- t(object@assays$RNA$counts[,object@meta.data[[clustering]] == cluster])
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

  degs <- data$gene[data$f_m_q.value<0.1|
                      data$d_f_q.value<0.1|
                      data$d_m_q.value<0.1]
  degs <- degs[!is.na(degs)]
  
  if(length(degs) ==0){return(NULL)}

degs_expression <- data.frame()
for(i in degs){
  gene_expression <-mean_expression_cluster_data(
    obj,
    i,
    cluster
    )

  new_data <-data.frame(
    cluster= cluster,
    gene = i,
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

  continuity <- function(degs, cluster, clustering = 'harmony.wnn_res0.4_clusters', object = obj){
       
   mean_expression_cluster_data <- function(object = object, gene, cluster, clustering = 'harmony.wnn_res0.4_clusters'){

      counts <- t(object@assays$RNA$counts[,object@meta.data[[clustering]] == cluster])
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
  results$mean <- ifelse(is.na(results$mean)| is.nan(results$mean), 0, results$mean)
  return(results)
   }


  
    pca_data <- data.frame()
for(i in degs){


  data <- mean_expression_cluster_data(object,
                                       i,
                                       cluster)


  data$gene <- i
  pca_data <- rbind(pca_data, data)
}

pca_data_pivoted <- pca_data %>%
  dplyr::select(individual, mean, Sex, gene)%>%
  dplyr::filter(Sex == 'D'| Sex =='F' | Sex == 'M')%>%
  pivot_wider(names_from = gene, 
              values_from = mean)

cluster_prcomp <- prcomp(pca_data_pivoted[,3:ncol(pca_data_pivoted)])

cluster_pca_loadings <- cluster_prcomp$x
       
cluster_pca_loadings <- as.data.frame(cluster_pca_loadings[,1:2])
cluster_pca_loadings$Sex <- pca_data_pivoted$Sex


grouped_means <- cluster_pca_loadings %>%
  group_by(Sex) %>%
  summarize(across(starts_with("PC"), base::mean))


mean_m <- grouped_means[grouped_means$Sex=='M',2:3]

mean_f <- grouped_means[grouped_means$Sex=='F',2:3]

mean_d <- grouped_means[grouped_means$Sex=='D',2:3]

f_m_distance <- stats::dist(rbind(as.numeric(mean_m), as.numeric(mean_f)))

d_m_distance <- stats::dist(rbind(as.numeric(mean_m), as.numeric(mean_d)))

d_f_distance <- stats::dist(rbind(as.numeric(mean_d), as.numeric(mean_f)))

continuum_score <- f_m_distance/(d_m_distance+d_f_distance)

return(continuum_score)
}
  
}



obj <- readRDS('/Volumes/jrhodes/Fish Lab/Experiments/sex change single nuc POA/R/RNA Object.rds')

directory <- dir('/Users/ggraham/Desktop/snRNA-seq R Files 122524/Seurat Outputs/122324 Neg Bin with Doms')

#read in resuls
for(i in directory){
  
  d <- read.csv(
    paste0('/Users/ggraham/Desktop/snRNA-seq R Files 122524/Seurat Outputs/122324 Neg Bin with Doms/',
           i
    )
  )
  
  assign(i, d, envir = .GlobalEnv)
}

for(i in 0:31){
  
  file <-  get(paste0('cluster_',
                       i,
                       '.csv'
  ))
  new_name <- paste0('results_cluster',i)
   assign(new_name, file,envir=.GlobalEnv)
}


#make list of dataframe names
file.list <- c()
for(i in 0:31){
  
  file.name <-  paste0('results_cluster',i)
  
  file.list <- append(file.list, file.name)
}


#### lasso analysis ####
probability_data <- data.frame()
for(i in 0:31){
  print(i)
  file.name <-  paste0('results_cluster',
                       i
  )
  data<- get(file.name)
  prediction <-#tryCatch({
    lasso_scorer(data, i)
                 #)}, error=function(e){
    #message('error')
    #return(NULL)})
if (is.null(prediction)) next
  prediction$cluster = i
  probability_data <- rbind(probability_data,prediction)

}


###Continuum analysis ####
continuity_data <- data.frame()
for(i in 0:31){
  print(i)
  cluster_data <- read.csv(paste0('/Users/ggraham/Desktop/snRNA-seq R Files 122524/Seurat Outputs/122324 Neg Bin with Doms/cluster_',i,'.csv'))

  cluster_degs<- cluster_data$gene[cluster_data$f_m_q.value<0.1|
                                         cluster_data$d_m_q.value<0.1|
                                          cluster_data$d_f_q.value<0.1]
  
cluster_degs <- cluster_degs[!is.na(cluster_degs)]
if(length(cluster_degs)<2){next}

continuity_score <- continuity(degs = cluster_degs, cluster = i)

output <- data.frame(cluster = i, 
                     continuity_score = continuity_score)

continuity_data<- rbind(output, continuity_data)

  
}
degs <- data.frame()
for(i in 0:31){
  print(i)
  cluster_data <- read.csv(paste0('/Users/ggraham/Desktop/snRNA-seq R Files 122524/Seurat Outputs/122324 Neg Bin with Doms/cluster_',i,'.csv'))

  cluster_degs<- cluster_data$gene[cluster_data$f_m_q.value<0.1|
                                         cluster_data$d_m_q.value<0.1|
                                          cluster_data$d_f_q.value<0.1]
  
cluster_degs <- cluster_degs[!is.na(cluster_degs)]
if(length(cluster_degs)<2){next}


output <- data.frame(cluster = i, 
                     n_degs = length(cluster_degs))
degs<- rbind(output, degs)

  
}


                  
##### Combining the two ####
colnames(probability_data) <- c('s0','predicted','fish','probabilities','cluster')
probability_data_grouped <- probability_data%>%
  group_by(cluster)%>%
  summarize(mean_pred = mean(probabilities),
            se = sd(probabilities)/sqrt(n()))

joined_data <- right_join(probability_data_grouped,continuity_data, by = 'cluster' )
joined_data <- right_join(joined_data,degs, by = 'cluster' )


probability_data <- read.csv('/Volumes/jrhodes/Fish Lab/Experiments/sex change single nuc POA/Seurat Outputs/122524 Dominant sex predictions lasso degs.csv')



library(forcats)
joined_data$clusters_of_interest <- ifelse(joined_data$cluster%in% c(2,3,6,7,8,14,19,29), 'interest',NA)

ggplot(joined_data, aes(x = continuity_score, y = mean_pred, color = clusters_of_interest))+
  geom_point(aes(size = n_degs))+
  geom_pointrange(data = joined_data, aes(x =continuity_score, ymin = mean_pred-se, ymax = mean_pred+se))+
  geom_text(label = joined_data$cluster, vjust = 1.5, color = 'black')+
  geom_hline(yintercept = 0.5, linetype = 2)+
    geom_vline(xintercept = 0.5, linetype = 2)+
  geom_text(x = 0.1, y = 0.05, label = 'nonlinear masculine', color = 'black')+
  geom_text(x = 0.90, y = 0.05, label = 'linear masculine', color = 'black')+
    geom_text(x = 0.1, y = 0.95, label = 'nonlinear feminine', color = 'black')+
      geom_text(x = 0.9, y = 0.95, label = 'linear feminine', color = 'black')+
  xlim(0,1)+
  ylim(0,1)+
  labs(x = 'Continuity Score', y = 'Lasso Prediction')
  



  
  
  joined_data$color <- ifelse(joined_data$mean_pred>0.6,'female',NA)
joined_data$color <- ifelse(joined_data$mean_pred<0.4,'male',joined_data$color)
joined_data$color <- ifelse(is.na(joined_data$color), 'intermediate', joined_data$color)

   ggplot(joined_data, aes(x = fct_reorder(as.factor(cluster),mean_pred), y = mean_pred, color = color)) +
  geom_pointrange(aes(ymin = mean_pred - se, ymax = mean_pred + se, size = n_degs)) +
  geom_hline(yintercept = 0.5, linetype = 2) +
     scale_size(range =c(0,1))

   
     joined_data$color_cont <- ifelse(joined_data$continuity_score>0.6,'linear',NA)
joined_data$color_cont <- ifelse(joined_data$continuity_score<0.4,'nonlinear',joined_data$color_cont)
joined_data$color_cont <- ifelse(is.na(joined_data$color_cont), 'intermediate', joined_data$color_cont)

      ggplot(joined_data, aes(x = fct_reorder(as.factor(cluster),continuity_score), y = continuity_score, color = color_cont)) +
  geom_point(data =joined_data, aes(size = mean_pred))+
  geom_hline(yintercept = 0.5, linetype = 2) 
      
      ggplot(joined_data, aes(x =cluster, y = continuity_score))+
  geom_bar(stat = 'identity', aes(fill = color_cont))+
  scale_x_continuous(breaks = 0:31)+
  geom_hline(yintercept = 0.5, linetype =2)
      
      
  #### Do number of DEGs predict either of these measures
  
      continuity_lm <- lm(continuity_score ~n_degs, data = joined_data)
      summary(continuity_lm)
      
            lasso_lm <- lm(mean_pred ~n_degs, data = joined_data)
                  summary(lasso_lm)
                  
                  #does lasso predict continuity or vice versa
                  
      corr_lm <- lm(continuity_score~mean_pred, data =joined_data)
summary(corr_lm)
#nope

#### How do they look on the UMAP ####
cluster_to_predicted_sex <- setNames(joined_data$mean_pred, joined_data$cluster)

cell_clusters <- Idents(obj)
metadata <- cluster_to_predicted_sex[as.character(cell_clusters)]
names(metadata) <- names(cell_clusters)
obj <- AddMetaData(obj, metadata, col.name = "mean_pred")

FeaturePlot(obj, 'mean_pred')+
    scale_color_gradientn(colors = c('blue','lightgreen','red'))

 cluster_to_corr <- setNames(joined_data$continuity_score, joined_data$cluster)
metadata2 <- cluster_to_corr[as.character(cell_clusters)]
  names(metadata2) <- names(cell_clusters)
  obj <- AddMetaData(obj, metadata2, col.name = "continuity_score")
 
  FeaturePlot(obj, 'continuity_score')+
    scale_color_gradientn(colors = c('blue','lightgreen','red'))

  ggplot(joined_data, aes(x = cluster, y = n_degs))+
    stat_summary(geom = 'bar', fun = 'identity')+
    scale_x_continuous(breaks = c(0:31))

  
prob_data_fish <- probability_data%>%
  group_by(fish)%>%
  summarize(mean_pred = mean(probabilities),
            se_pred = sd(probabilities)/sqrt(n()))

prob_data_fish$predicted_sex <- ifelse(prob_data_fish$mean_pred>0.66, 'F',NA)
prob_data_fish$predicted_sex <- ifelse(prob_data_fish$mean_pred<0.33, 'M',prob_data_fish$predicted_sex)
prob_data_fish$predicted_sex <- ifelse(is.na(prob_data_fish$predicted_sex), 'Intermediate',prob_data_fish$predicted_sex)

ggplot(prob_data_fish, aes(x = fct_reorder(fish, mean_pred), y = mean_pred, color = predicted_sex))+
  geom_pointrange(aes(fct_reorder(fish, mean_pred), y = mean_pred, ymin = mean_pred-se_pred, ymax =mean_pred+se_pred))
      
joined_data$continuity_class <- ifelse(joined_data$continuity_score>0.66, 'linear', NA)
joined_data$continuity_class <- ifelse(joined_data$continuity_score<0.33, 'nonlinear', joined_data$continuity_class)
joined_data$continuity_class <- ifelse(is.na(joined_data$continuity_class), 'intermediate', joined_data$continuity_class)


  ggplot(joined_data, aes(x = fct_reorder(as.factor(cluster),continuity_score), y = continuity_score, color =continuity_class))+
           geom_point(aes(size = n_degs))
  