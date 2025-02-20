### correlate continuity score and eigengene score

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
  library(forcats)
  `%notin%` <- Negate(`%in%`)
  
}


obj <- readRDS("C:/Users/Gabe/Desktop/RNA Object.rds")

##### continuity score ####
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
  
  ### run ANCOVA on pc loadings ###
  #initialize list
  ancovas <- list()
  cluster_pca_loadings_ancova <- cluster_pca_loadings
  cluster_pca_loadings_ancova$dummy_sex <- ifelse(cluster_pca_loadings_ancova$Sex == 'D', "Test", "Known")
  
  pc1_by_pc2 <- lm(PC1 ~ PC2+dummy_sex, data =cluster_pca_loadings_ancova)
  pc2_by_pc1 <- lm(PC2 ~ PC1+dummy_sex, data =cluster_pca_loadings_ancova)
  ancovas[['pc1_by_pc2']] =pc1_by_pc2
  ancovas[['pc2_by_pc1']] =pc2_by_pc1
  
  #assign to global env
  assign(paste0('ancovas_',cluster), ancovas, envir = .GlobalEnv)
  
  #continue with continuity score
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
###coalesce ancova_data
ancova_data <- data.frame()
for(i in 0:31){
  ancovas_list <- get(paste0('ancovas_',i))
  model_1_2 <- ancovas_list[['pc1_by_pc2']]
  model_2_1 <- ancovas_list[['pc2_by_pc1']]
  
  av_1_2 <- as.data.frame(anova(model_1_2, test = 'Chisq'))[-3,]
  av_2_1 <- as.data.frame(anova(model_2_1, test = 'Chisq'))[-3,]
  
  av_1_2$model = '1_2'
  av_1_2$vars = c('PC2','dummy_sex')
  av_2_1$model = '2_1'
  av_2_1$vars = c('PC1','dummy_sex')
  
  new_data <- rbind(av_1_2,av_2_1)
  new_data$cluster = i
  ancova_data <- rbind(ancova_data, new_data)
}

ancova_data$issignif <- ifelse(ancova_data$`Pr(>F)`<0.05, '*', NA)

ancova_data[ancova_data$vars!='dummy_sex',]

ancova_data_dummy_sex <- ancova_data[ancova_data$vars=='dummy_sex',]

continuity_data_2_plot <- continuity_data_2%>%
  group_by(cluster)%>%
  summarize(mean_pred = mean(continuity_score.continuum_score_individual),
            se_pred = sd(continuity_score.continuum_score_individual)/sqrt(n())
  )%>%
  right_join(ancova_data_dummy_sex, by = 'cluster')

continuity_data_2_plot$Linearity <- ifelse(continuity_data_2_plot$mean_pred<0.33, 'Nonlinear', NA)
continuity_data_2_plot$Linearity <- ifelse(continuity_data_2_plot$mean_pred>0.66, 'Linear', continuity_data_2_plot$Linearity)
continuity_data_2_plot$Linearity <- ifelse(is.na(continuity_data_2_plot$Linearity ), 'Intermediate', continuity_data_2_plot$Linearity)

continuity_plot2 <- ggplot(subset(continuity_data_2_plot, cluster %notin% c(15,30)), aes(x =fct_reorder(as.factor(cluster), mean_pred), y = mean_pred, color = Linearity, shape = Linearity))+
  geom_pointrange(aes(x =fct_reorder(as.factor(cluster), mean_pred), y = mean_pred, ymin = mean_pred-se_pred, ymax = mean_pred+se_pred))+
  theme_classic()+
  labs(x = 'Cluster', y = 'Continuity Score +/- SE')+
  theme(legend.position = c(0.15,.85))+
  theme(axis.text.x = element_text(size = 10,angle = -45, vjust = 1, hjust=0), axis.text.y = element_text(size = 10), axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12))+
  geom_text(aes(label = issignif,x =fct_reorder(as.factor(cluster), mean_pred), y = 1.05*(mean_pred+se_pred)), color = 'black', size=10, show.legend = F)
continuity_plot2

data_to_count <- subset(continuity_data_2_plot, cluster %notin% c(15,30))
data_to_count$look <- ifelse(data_to_count$cluster %in% c(19, 2, 8,7, 6, 3, 14), '*', NA)
data_to_count_look <-data_to_count[!is.na(data_to_count$look),]
length(unique(data_to_count$cluster))


#ggsave(plot = continuity_plot2,
#       file = "continuity_plot_individual.svg",
#       device = "svg",
#       units = "in",
#       width = 5,
#       height = 5,
#       path = "Bachelors Thesis/Plots/Figure 2")



##### eigengene score ######
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


##### summarize and join ####
pca_predictions_merged_summarized <- pca_predictions_merged%>%
  group_by(cluster)%>%
  summarize(mean_sex = mean(predicted),
            se_sex = sd(predicted)/sqrt(n())
  )
continuity_data_2_grouped <- continuity_data_2%>%
  group_by(cluster)%>%
  summarize(mean_continuity = mean(continuity_score.continuum_score_individual),
            se_continuity = sd(continuity_score.continuum_score_individual)/sqrt(n())
  )

continuity_sex_joined <- pca_predictions_merged_summarized%>%
  right_join(continuity_data_2_grouped, by = 'cluster')

cor.test(continuity_sex_joined$mean_continuity, continuity_sex_joined$mean_sex)
#no correlation 
continuity_sex_joined$cluster_label <- ifelse((continuity_sex_joined$mean_sex <0.66) &(continuity_sex_joined$mean_continuity <0.66),continuity_sex_joined$cluster, NA) 
continuity_sex_joined$color <- ifelse(continuity_sex_joined$cluster %in% c(2, 3, 6,7,8,14,19), 'Interest', 'Other')

continuity_sex_joined_plot <- ggplot(subset(continuity_sex_joined, cluster %notin% c(15,30)), aes(x = mean_continuity, y = mean_sex, color = color, shape = color))+
  geom_pointrange(aes(x = mean_continuity, y = mean_sex, ymin = mean_sex - se_sex, ymax = mean_sex+se_sex), linetype =2)+
  geom_pointrange(aes(x = mean_continuity, y = mean_sex, xmin = mean_continuity - se_continuity, xmax = mean_continuity+se_continuity), linetype = 2)+
  geom_text_repel(aes(label = cluster), size =5, show.legend = F)+
  theme_classic()+
  theme(axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10), axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12))+
  scale_color_manual(values = c('red','black'))+
  labs(x = 'Mean Continuity Score +/- SE', y = 'Mean Sex Prediction +/- SE', color = 'Clusters', shape = 'Clusters')+
  theme(legend.position = c(0.9,.85))+
  xlim(0,1)+
  ylim(0,1)

#ggsave(plot = continuity_sex_joined_plot,
#       file = "continuity_sex.svg",
#       device = "svg",
#       units = "in",
#       width = 5,
#       height = 5,
#       path = "Bachelors Thesis/Plots/Figure 2")

continuity_sex_joined$both <- ifelse((continuity_sex_joined$mean_sex<0.66 & continuity_sex_joined$mean_continuity<0.66), '*', NA)


