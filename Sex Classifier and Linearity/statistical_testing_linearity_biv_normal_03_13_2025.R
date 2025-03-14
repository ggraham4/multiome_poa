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

### define a function to test the linearity of each individual
continuity_individual <- function(degs, cluster, clustering = 'harmony.wnn_res0.4_clusters', object = obj){
  
  #define a function to summarize the expression of genes for each individual
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
    #for each Deg, summarize expression for each individual
    data <- mean_expression_cluster_data(object,
                                         i,
                                         cluster)
    
    data$gene <- i
    pca_data <- rbind(pca_data, data)
  }
  
  #pivot such that columns ar indiviidual, sex, and the mean expression for each gene
  pca_data_pivoted <- pca_data %>%
    dplyr::select(individual, mean, Sex, gene)%>%
    dplyr::filter(Sex == 'D'| Sex =='F' | Sex == 'M')%>% #might be interesting to interrogate where Es and NFs fall on this 
    pivot_wider(names_from = gene, 
                values_from = mean)
  
  #Run PCA
  cluster_prcomp <- prcomp(pca_data_pivoted[,3:ncol(pca_data_pivoted)])
  
  #Extract their coordinates
  cluster_pca_loadings <- cluster_prcomp$x
  cluster_pca_loadings <- as.data.frame(cluster_pca_loadings[,1:2])
  cluster_pca_loadings$Sex <- pca_data_pivoted$Sex
  
  ### run Bivariate normal on this data ###
  male_female_data <- cluster_pca_loadings%>%
    subset(Sex != 'D')
  
  #Extract mean PC1 and PC2
  mu_hat <- colMeans(male_female_data[,-3])
  #Extract Covariance from data
  Sigma_hat <- cov(male_female_data[,-3])
  
  #simulate data from males and females
  simulated_data <- mvrnorm(n = 10000, mu = mu_hat, Sigma = Sigma_hat)
  
  #calculate p-value for mean dominant
  dom_data <- cluster_pca_loadings%>%
    subset(Sex == 'D')
  
  dom_PC1 <- mean(dom_data[,1])
  dom_PC2 <- mean(dom_data[,2])
  point <- c(dom_PC1, dom_PC2)
  
  #calculate mahalanobis distance
  mahalanobis_dist <- mahalanobis(point, mu_hat, Sigma_hat)
  
  # Compute p-value using chi-square distribution with 2 degrees of freedom
  p_value <- 1 - pchisq(mahalanobis_dist, df = 2)
  p_value
  
  ### continue with continuity score ###
  grouped_means <- cluster_pca_loadings %>%
    group_by(Sex) %>%
    summarize(across(starts_with("PC"), base::mean))
  
  #summarize mean male and female
  mean_m <- grouped_means[grouped_means$Sex=='M',2:3]
  
  mean_f <- grouped_means[grouped_means$Sex=='F',2:3]

  #calculate the distance between the mean male and female
  f_m_distance <- stats::dist(rbind(as.numeric(mean_m), as.numeric(mean_f)))
  
  # add individual into this cause I'm gonna query each dominant individually
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
      continuum_score_individual =continuum_score_individual,
      p_value = p_value
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

continuity_data_2_summarized <- continuity_data_2%>%
  group_by(cluster)%>%
  summarize(mean_continuum_score = mean(continuity_score.continuum_score_individual),
            se =  sd(continuity_score.continuum_score_individual)/sqrt(n()),
            p_value = mean(continuity_score.p_value))


continuity_data_2_summarized$issignif <- ifelse(continuity_data_2_summarized$p_value<0.05, '*', NA)
continuity_data_2_summarized$Linearity <- ifelse(continuity_data_2_summarized$mean_continuum_score<0.33, 'Nonlinear', NA)
continuity_data_2_summarized$Linearity <- ifelse(continuity_data_2_summarized$mean_continuum_score>0.66, 'Linear', continuity_data_2_summarized$Linearity)
continuity_data_2_summarized$Linearity <- ifelse(is.na(continuity_data_2_summarized$Linearity ), 'Intermediate', continuity_data_2_summarized$Linearity)

###plot this result
continuity_plot2 <- ggplot(subset(continuity_data_2_summarized, cluster %notin% c(15,30)), aes(x =fct_reorder(as.factor(cluster), mean_continuum_score), y = mean_continuum_score, color = Linearity, shape = Linearity))+
  geom_pointrange(aes(x =fct_reorder(as.factor(cluster), mean_continuum_score), y = mean_continuum_score, ymin = mean_continuum_score-se, ymax = mean_continuum_score+se))+
  theme_minimal()+
  labs(x = 'Cluster', y = 'Linearity Score +/- SE')+
  theme(legend.position = c(0.17,.77))+
  theme(axis.text.x = element_text(size = 10,angle = -90, vjust = 1, hjust=0), axis.text.y = element_text(size = 10), axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12))+
  geom_text(aes(label = issignif,x =fct_reorder(as.factor(cluster), mean_continuum_score), y = 1.05*(mean_continuum_score+se)), color = 'black', size=10, show.legend = F)
continuity_plot2

write.csv(continuity_data_2_summarized, 'Sex Classifier and Linearity/linearity_score_03_13_2025.csv')

#ggsave(plot = continuity_plot2,
#       file = "continuity_plot_individual.svg",
#       device = "svg",
 #      units = "in",
  #     width = 4.5,
   #    height = 2.5,
    #   path = "Bachelors Thesis/Plots/Figure 3")




