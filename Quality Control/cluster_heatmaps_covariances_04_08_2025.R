library(Seurat)
library(tidyverse)
library(pheatmap)
library(ggplot2)

obj <- readRDS('C:/Users/Gabe/Desktop/RNA Object.rds')

### zack email 04_04_2025

#matrix 1 - proportion of each cluster within each individual
#> another way to state this, what % of cluster i is within a given individual j

individuals <- unique(obj$individual)

proportion_of_each_cluster_within_each_individual <- matrix(0, length(individuals), 32)
rownames(proportion_of_each_cluster_within_each_individual) <- individuals
colnames(proportion_of_each_cluster_within_each_individual) <- 0:31

for(c in 0:31){
  print(c)
  #calculate total cells in cluster
  total_cells <- nrow(obj@meta.data[obj$harmony.wnn_res0.4_clusters==c,])
  
  cluster_index <- which(colnames(proportion_of_each_cluster_within_each_individual)==c)
  for(i in individuals){
    #calculate total cells in individual in cluster
    total_cells_individual <- nrow(obj@meta.data[obj$harmony.wnn_res0.4_clusters==c & obj$individual==i,])
    
    prop <- total_cells_individual/ total_cells
    
    individual_index <- which(rownames(proportion_of_each_cluster_within_each_individual)==i)
    
    proportion_of_each_cluster_within_each_individual[individual_index,cluster_index] = prop
  }
}

pheatmap(log(proportion_of_each_cluster_within_each_individual+0.001),
         main='Log Proportion of Each Cluster Within Each Individual')

#0 and 1 covary pretty well
set.seed(123)
pheatmap(log(proportion_of_each_cluster_within_each_individual+0.001),
         k  = 5,
         main='Log Proportion of Each Cluster Within Each Individual K = 5') #here, clusters are individuals

poecwei_hierarchical <- hclust(dist(proportion_of_each_cluster_within_each_individual))
plot(poecwei_hierarchical, xlab ='Individual', sub = 'Complete Clustering')


# matrix 2 - proportion of each individual within each cluster
#what % of fish in cluster
proportion_of_each_individual_within_each_cluster <- matrix(0, length(individuals), 32)
rownames(proportion_of_each_individual_within_each_cluster) <- individuals
colnames(proportion_of_each_individual_within_each_cluster) <- 0:31

for(i in individuals){
  print(i)
  individual_index <- which(rownames(proportion_of_each_individual_within_each_cluster)==i)
  
  total_cells_in_fish <- nrow(obj@meta.data[obj$individual==i,])
  for(c in 0:31){
    cluster_index <- which(colnames(proportion_of_each_individual_within_each_cluster)==c)
    
    cells_in_fish_in_cluster <- nrow(obj@meta.data[obj$individual==i & obj$harmony.wnn_res0.4_clusters ==c,])
    
    prop = cells_in_fish_in_cluster / total_cells_in_fish
    
    proportion_of_each_individual_within_each_cluster[individual_index,cluster_index] = prop
    
  }
  
}

pheatmap(t(proportion_of_each_individual_within_each_cluster),
         main='Proportion of Each Individual Within Each Cluster')

pheatmap(t(proportion_of_each_individual_within_each_cluster),k = 5,
         main='Proportion of Each Individual Within Each Cluster K = 5')


poeiwec_hierarchical <- hclust(dist(t(proportion_of_each_individual_within_each_cluster)))
plot(poeiwec_hierarchical, xlab ='Cluster', sub = 'Complete Clustering')

#again, 0 and 1 covary

#mean and SD of relative proportions for each cluster, compare to known data

mean_and_sd = as.data.frame(t(proportion_of_each_individual_within_each_cluster))
mean_and_sd$cluster = rownames(mean_and_sd)

mean_and_sd_summarized <- mean_and_sd %>%
  pivot_longer(cols = -cluster, names_to = "individual", values_to = "value") %>%
  group_by(cluster) %>%
  summarize(
    mean_value = mean(value, na.rm = TRUE),
    sd_value = sd(value, na.rm = TRUE)
  )

#plot of means and sds
ggplot(mean_and_sd_summarized, aes(x = fct_reorder(cluster,mean_value, .desc = T), y = mean_value))+
  geom_bar(stat = 'identity')+
  labs(x = 'Cluster', y = 'Mean Proportion of Individual Within Cluster')

ggplot(mean_and_sd_summarized, aes(x = fct_reorder(cluster,mean_value, .desc = T), y = sd_value))+
  geom_bar(stat = 'identity')+
  labs(x = 'Cluster', y = 'SD Proportion of Individual Within Cluster')

