#linearity_null_07_28_2025
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

###summarize deg data
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

### run PCA
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

### get real results
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


###create fake labels and permute
whole_cluster_null = data.frame()
for(clust in unique(pca_data$cluster)){
  message(clust)
  
  #subset data
  clust_data = subset(pca_data, cluster == clust & status != 'NRM')
  
  #extract statuses
  statuses= clust_data$status
  
  #make data frame to append null data into
  null_data = data.frame()
  for(iter in 1:1000){ #1000 times shuffle the statuses 
    
    #new column with fake stataus
  clust_data$fake_status = sample(statuses, length(statuses))
  
  ##perform calculation using fake status -- I'm only doing centroids
  #endpoints
  male_data = subset(clust_data, fake_status =='M')
  female_data = subset(clust_data, fake_status =='F')
  
  male_centroid = c(mean(male_data$PC1), mean(male_data$PC2))
  female_centroid = c(mean(female_data$PC1), mean(female_data$PC2))
  
  #query animals
  dominant_data = subset(clust_data, fake_status =='D')
  dominant_centroid = c(mean(dominant_data$PC1), mean(dominant_data$PC2))
  
  expanded_data =subset(clust_data, fake_status =='E')
  expanded_centroid = c(mean(expanded_data$PC1), mean(expanded_data$PC2))
  
  nf_data =subset(clust_data, fake_status =='NF')
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
  
  
  newd = data.frame(cluster =clust,
                    lin_dom=linearity_index_dominant,
                    lin_expanded =linearity_index_expanded,
                    lin_nf = linearity_index_nf,
                    iter = iter)
  
  null_data = rbind(null_data, newd)
  
  }
  whole_cluster_null = rbind(whole_cluster_null, null_data)
  
}

#still get the pileup at 1 why does this happen?
hist(whole_cluster_null$lin_dom)

hist(whole_cluster_null$lin_dom[whole_cluster_null$clust==6])

hist(log(whole_cluster_null$lin_dom)+0.01) 

#it does not produce a normal distribution because it is bounded, so I think it is actually fair

### perform one tailed p value calculation 

empirical_pvalue_dom = function(result){
  p_value = sum(whole_cluster_null$lin_dom <= result)/nrow(whole_cluster_null)
  return(p_value)
}


mean_data$dom_p.value <- lapply(X=mean_data$lin_dom, FUN = empirical_pvalue_dom)
mean_data$dom_p.value <- as.numeric(mean_data$dom_p.value)

hist(mean_data$dom_p.value)
mean_data$dom_signif = ifelse(mean_data$dom_p.value < 0.05, '*', NA)

hist(mean_data$lin_dom)

# well, I cant massage this into being significant in any way, is this really worth including in 
# the paper


empirical_pvalue_exp = function(result){
  p_value = sum(whole_cluster_null$lin_expanded <= result)/nrow(whole_cluster_null)
  return(p_value)
}
mean_data$exp_p.value <- lapply(X=mean_data$lin_expanded, FUN = empirical_pvalue_exp)
mean_data$exp_signif = ifelse(mean_data$exp_p.value < 0.05, '*', NA)

