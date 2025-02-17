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
  `%notin%` <- Negate(`%in%`)
  
  
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
obj <- readRDS("C:/Users/Gabe/Desktop/RNA Object.rds")


continuity_data <- data.frame()
for(i in 0:31){
  print(i)
  cluster_data <- read.csv(
    paste0('DEG Outputs/012425 Neg Bin w Doms Lower Stringency/cluster_',i,'.csv')
  )
  
  cluster_degs<- cluster_data$gene[cluster_data$f_m_q.value<0.1|
                                     cluster_data$d_m_q.value<0.1|
                                     cluster_data$d_f_q.value<0.1]
  
  cluster_degs <- cluster_degs[!is.na(cluster_degs)]
  if(length(cluster_degs)<2){next}
  
  counts <- t(obj@assays$RNA$counts[,obj@meta.data$harmony.wnn_res0.4_clusters == i])
  
  continuity_score <- continuity(degs = cluster_degs, cluster = i)
  
  output <- data.frame(cluster = i, 
                       continuity_score = continuity_score)
  
  continuity_data<- rbind(output, continuity_data)
  
  
}

continuity_data$Linearity <- ifelse(continuity_data$continuity_score<0.33, 'Nonlinear', NA)
continuity_data$Linearity <- ifelse(continuity_data$continuity_score>0.66, 'Linear', continuity_data$Linearity)
continuity_data$Linearity <- ifelse(is.na(continuity_data$Linearity ), 'Intermediate', continuity_data$Linearity)

continuity_plot <- ggplot(subset(continuity_data, cluster %notin% c(15,30)), aes(x =fct_reorder(as.factor(cluster), continuity_score), y = continuity_score, color = Linearity, shape = Linearity))+
  geom_point(size = 3)+
  theme_classic()+
  labs(x = 'Cluster', y = 'Continuity Score')+
  theme(legend.position = c(0.15,.85))+
  theme(axis.text.x = element_text(size = 8,angle = -45, vjust = 1, hjust=0), axis.text.y = element_text(size = 10), axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12))

ggsave(plot = continuity_plot,
       file = "continuity_plot.svg",
       device = "svg",
       units = "in",
       width = 5,
       height = 5,
       path = "Bachelors Thesis/Plots/Figure 2")
 
#the more I think about this, the more I think I could get a standard error by routing through 
#* each individual dominant as opposed to just the mean dominant
#* 
#* also, i'm convinced I can get a p value from this some other way than our permutation
#* analysis, I think you should be able to connect the distribution of males to females 
#* and then show that dominants lie outside of that distribution
#* 
#* 
####### with an individual term ######
continuity_individual <- function(degs, cluster, clustering = 'harmony.wnn_res0.4_clusters', object = obj){
  
  mean_expression_cluster_data <- function(object = object, gene, cluster, clustering = 'harmony.wnn_res0.4_clusters'){
    
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
  
  cluster_pca_loadings$individual = pca_data_pivoted$individual
  
  individual_data <- data.frame()
  for(u in unique(cluster_pca_loadings$individual[cluster_pca_loadings$Sex == 'D'])){
    
    dominant_data <- cluster_pca_loadings[cluster_pca_loadings$individual == u,1:2]
    
    
    dominant_m_distance <- stats::dist(rbind(as.numeric(mean_m), as.numeric(dominant_data)))
    
    dominant_f_distance <- stats::dist(rbind(as.numeric(dominant_data), as.numeric(mean_f)))
    
    continuum_score_individual <- f_m_distance/(dominant_m_distance+dominant_f_distance)
    
    new_data <- data.frame(
      cluster = cluster,
      individual = u,
      continuum_score_individual =continuum_score_individual
    )
    individual_data <- rbind(individual_data, new_data)
    
  }  
  return(individual_data)
}

continuity_data_2 <- data.frame()
for(i in 0:31){
  print(i)
  cluster_data <- read.csv(
    paste0('DEG Outputs/012425 Neg Bin w Doms Lower Stringency/cluster_',i,'.csv')
  )
  
  cluster_degs<- cluster_data$gene[cluster_data$f_m_q.value<0.1|
                                     cluster_data$d_m_q.value<0.1|
                                     cluster_data$d_f_q.value<0.1]
  
  cluster_degs <- cluster_degs[!is.na(cluster_degs)]
  if(length(cluster_degs)<2){next}
  
  counts <- t(obj@assays$RNA$counts[,obj@meta.data$harmony.wnn_res0.4_clusters == i])
  
  continuity_score <- continuity_individual(degs = cluster_degs, cluster = i)
  
  output <- data.frame(cluster = i, 
                       continuity_score = continuity_score)
  
  continuity_data_2<- rbind(output, continuity_data_2)
  
  
}

colnames(continuity_data_2)

continuity_data_2_grouped <- continuity_data_2%>%
  group_by(cluster)%>%
  summarize(mean_pred = mean(continuity_score.continuum_score_individual),
            se_pred = sd(continuity_score.continuum_score_individual)/sqrt(n())
  )

continuity_data_2_grouped$Linearity <- ifelse(continuity_data_2_grouped$mean_pred<0.33, 'Nonlinear', NA)
continuity_data_2_grouped$Linearity <- ifelse(continuity_data_2_grouped$mean_pred>0.66, 'Linear', continuity_data_2_grouped$Linearity)
continuity_data_2_grouped$Linearity <- ifelse(is.na(continuity_data_2_grouped$Linearity ), 'Intermediate', continuity_data_2_grouped$Linearity)

continuity_plot2 <- ggplot(subset(continuity_data_2_grouped, cluster %notin% c(15,30)), aes(x =fct_reorder(as.factor(cluster), mean_pred), y = mean_pred, color = Linearity, shape = Linearity))+
  geom_pointrange(aes(x =fct_reorder(as.factor(cluster), mean_pred), y = mean_pred, ymin = mean_pred-se_pred, ymax = mean_pred+se_pred))+
  theme_classic()+
  labs(x = 'Cluster', y = 'Continuity Score')+
  theme(legend.position = c(0.15,.85))+
  theme(axis.text.x = element_text(size = 8,angle = -45, vjust = 1, hjust=0), axis.text.y = element_text(size = 10), axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12))

joined_cont <- continuity_data_2_grouped%>%right_join(continuity_data, by = 'cluster')

cor.test(joined_cont$mean_pred, joined_cont$continuity_score)
ggplot(joined_cont,aes(x=mean_pred, y=continuity_score))+
  geom_point()+
  geom_text(aes(label = cluster))
#i guess this actually makes sense? the mean coordinates of dominants is not the same as the mean dominant score I dont think
#* we are doing this transformation to bound it between 0 and 1 so I think it doesnt translate 1:1? I need to look more closely at it

ggsave(plot = continuity_plot2,
       file = "continuity_plot_individual.svg",
       device = "svg",
       units = "in",
       width = 5,
       height = 5,
       path = "Bachelors Thesis/Plots/Figure 2")
