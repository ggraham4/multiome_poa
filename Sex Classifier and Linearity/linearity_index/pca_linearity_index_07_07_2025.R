#### linearity index ###
#libaries
library(tidyverse)
library(Seurat)
library(factoextra)
mecd = readRDS("Functions/mean_expression_cluster_data.rds")

obj = readRDS("C:/Users/Gabe/Desktop/nemo.orig_harmony.integration_all_testd_clusters.rds")

### read in DEGs
deg_data = data.frame()
base_path = "DEG Outputs/05_11_2025 Neg Bin w Doms New_clustering/"
base_dir = dir("DEG Outputs/05_11_2025 Neg Bin w Doms New_clustering")
for(i in base_dir){
  data = read.csv(paste0(base_path, i))
  
  subset_data = subset(data, av_q.value < 0.1 & singular == F)
  
  data_to_append = subset_data%>%
    dplyr::select(gene, cluster)
  deg_data = rbind(deg_data, data_to_append)
  
}

#summarize deg data
expression_data_all_clusters = list()
for(i in unique(deg_data$cluster)){
  print(i)
  deg_data_subset = subset(deg_data, cluster ==i)
  expression_data = data.frame()
    for(g in unique(deg_data_subset$gene)){
      gex = mecd(object = obj,gene = g,cluster = i, clustering='res0.8_50nn_40PC_45LSI')
      
      gex$gene = g
      expression_data = rbind(expression_data, gex)
    }
  
  expression_data_pivoted = expression_data%>%
    dplyr::select(individual, Status, mean, gene)%>%
    pivot_wider(names_from = gene, values_from=mean)
  
  expression_data_all_clusters[[paste0('cluster_',i)]] = expression_data_pivoted
  
}

## run PCA
pca_data = data.frame()
for(i in unique(deg_data$cluster)){
  print(i)
  
  cluster_data = expression_data_all_clusters[[paste0('cluster_',i)]]
  
  pca_matrix = as.matrix(cluster_data[,3:ncol(cluster_data)])
  
  if(ncol(pca_matrix)<2){next}
  pca = prcomp(pca_matrix, scale = T)
  
  pc_1 = pca$x[,1]
  pc_2 = pca$x[,2]
  
  loading_data = data.frame('individual' = cluster_data$individual,
                            'status' = cluster_data$Status,
                            'PC1' = pc_1,
                            'PC2' = pc_2,
                            'cluster' = i)
  pca_data = rbind(pca_data, loading_data)
  
}

#fviz_pca_ind(pca, habillage = cluster_data$Status)

## do linearity calculation
mean_data =data.frame()
for(i in unique(pca_data$cluster)){
  subset_data = subset(pca_data, cluster == i)
  
  #endpoints
  male_data = subset(subset_data, status =='M')
  female_data = subset(subset_data, status =='F')
  
  male_centroid = c(mean(male_data$PC1), mean(male_data$PC2))
  female_centroid = c(mean(female_data$PC1), mean(female_data$PC2))
  
  #query animals
  dominant_data = subset(subset_data, status =='D')
  dominant_centroid = c(mean(dominant_data$PC1), mean(dominant_data$PC2))
  
  expanded_data =subset(subset_data, status =='E')
  expanded_centroid = c(mean(expanded_data$PC1), mean(expanded_data$PC2))
  
  nf_data =subset(subset_data, status =='NF')
  nf_centroid = c(mean(nf_data$PC1), mean(nf_data$PC2)) 
  
  #calculate linearity scores
  m_f_d_matrix = matrix(c(male_centroid,
                        female_centroid,
                        dominant_centroid,
                        expanded_centroid,
                        nf_centroid),5,2,
                        byrow =T)
  
  
  m_f_distance = stats::dist(m_f_d_matrix[c(1,2),], method = 'euclidean')[1]
  m_d_distance = stats::dist(m_f_d_matrix[c(1,3),], method = 'euclidean')[1]
  f_d_distance = stats::dist(m_f_d_matrix[c(2,3),], method = 'euclidean')[1]
  
  m_e_distance= stats::dist(m_f_d_matrix[c(1,4),], method = 'euclidean')[1]
  f_e_distance = stats::dist(m_f_d_matrix[c(2,4),], method = 'euclidean')[1]
  
  m_nf_distance= stats::dist(m_f_d_matrix[c(1,5),], method = 'euclidean')[1]
  f_nf_distance = stats::dist(m_f_d_matrix[c(2,5),], method = 'euclidean')[1]
  
  linearity_index_dominant = m_f_distance/(m_d_distance+f_d_distance)
  linearity_index_expanded = m_f_distance/(m_e_distance+f_e_distance)
  linearity_index_nf = m_f_distance/(m_nf_distance+f_nf_distance)
  
  newd = data.frame(cluster =i,
                    m_f_distance = m_f_distance,
                    m_d_distance = m_d_distance,
                    f_d_distance = f_d_distance,
                    m_e_distance = m_e_distance,
                    f_e_distance=f_e_distance,
                    m_nf_distance=m_nf_distance,
                    f_nf_distance=f_nf_distance,
                    lin_dom=linearity_index_dominant,
                    lin_expanded =linearity_index_expanded,
                    lin_nf = linearity_index_nf)
  
  mean_data = rbind(mean_data, newd)
}

#plot for fun 
library(ggplot2)
ggplot(mean_data, aes(x = as.factor(cluster), y = lin_dom))+
  geom_point(size =4, shape=15, aes(colour = "D"))+
  geom_point(aes(y = lin_expanded, colour = "E"), shape = 17, size = 4)+
  geom_point(aes(y = lin_nf, colour = "NF"), shape = 16, size = 4)+
  ylim(0,1)

  
ggplot(subset(mean_data, cluster %in% c(0,1,6,11)), aes(x = as.factor(cluster), y = lin_dom))+
  geom_point(size =4, shape=15, aes(colour = "D"))+
  geom_point(aes(y = lin_expanded, colour = "E"), shape = 17, size = 4)+
  geom_point(aes(y = lin_nf, colour = "NF"), shape = 16, size = 4)+
  ylim(0,1)

ggplot(mean_data, aes(x = as.factor(cluster), y = lin_dom))+
  geom_point(size =4, shape=15, aes(colour = "D"))+
  ylim(0,1)




### calculate se of doms and expandeds

### also, i think it makes sense in this file to see if I can calculate a progress percent 
## by seeing what % of the mf distance the df distance is where a smaller number indicates closer
#to female so maybe 1/that number?
mean_data$dom_progress = mean_data$f_d_distance/mean_data$m_f_distance

ggplot(mean_data, aes(x = as.factor(cluster), y = dom_progress))+
  geom_point(size =4, shape=15, aes(colour = "D"))
# so now this says all the clusters are male, I think that is probably not true and there
# is so much noise in this because of the whole nonlinearity dimension, if dominants
#> are equidistant to males and females, they will get a score of 1 but that is not true
#> so you could probably also do something with the female distance, and im realizing
#> i am just recreating the linearity score arent i thats pretty funny
#> Maybe if you did d_f / d_m you could see which one it is closer to?

mean_data$dom_progress2 = mean_data$f_d_distance/mean_data$m_d_distance

ggplot(mean_data, aes(x = as.factor(cluster), y = dom_progress2))+
  geom_point(size =4, shape=15, aes(colour = "D"))

