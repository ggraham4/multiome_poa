library(Seurat)
library(dplyr)
library(tidyverse)
### read in data
obj = readRDS("~/Desktop/optimal_clustering_rna_only.rds")
deg_csvs <- dir('DEG Outputs/05_11_2025 Neg Bin w Doms New_clustering')
deg_path <- 'DEG Outputs/05_11_2025 Neg Bin w Doms New_clustering/'

# read in DEGs list, I think first I will do all singular and nonsingular
all_degs <- data.frame()
for(i in deg_csvs){
  data <- read.csv(paste0(deg_path,i))
  all_degs <- rbind(all_degs, data)
}

all_degs_subset <- subset(all_degs, av_q.value <0.1 & f_m_p.value <0.05 & singular ==F)

mecd = readRDS("Functions/mean_expression_cluster_data.rds")

for(clusterx in unique(all_degs_subset$cluster)){
  deg_mean_data = data.frame()
  print(clusterx)
  subset_cluster = subset(all_degs_subset, cluster == clusterx)
  degs = unique(subset_cluster$gene)
  for(deg in degs){
    deg_data = mecd(obj,
                    deg,
                    clusterx, 
                    clustering = 'final_clusters')
    deg_data$gene = deg
    
    deg_mean_data = rbind(deg_mean_data, deg_data)
  }
  pivoted_data = deg_mean_data%>%
    dplyr::select(-c(se, Sex))%>%
    pivot_wider(names_from = gene, values_from = mean)%>%
    filter(Status != 'NRM')
  path = 'Sex Classifier and Linearity/progress_classifier/degs_0.1_f_m/'
  file_name = paste0('cluster_',clusterx,'.csv')
  
  write.csv(pivoted_data, paste0(path, file_name))
  
}

