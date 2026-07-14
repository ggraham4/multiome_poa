library(Seurat)
library(dplyr)
library(tidyverse)

mecd = function(object, gene, cluster, clustering = 'final_clusters'){
  library(stringr)    
  options(dplyr.summarise.inform = FALSE)

    counts <- t(object@assays$RNA$data[,object@meta.data[[clustering]] == cluster]%>%as.matrix())
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

mecd_ATAC = function(object, gene, cluster, clustering = 'final_clusters'){
  library(stringr)    
  options(dplyr.summarise.inform = FALSE)

    counts <- t(object@assays$ATAC$data[,object@meta.data[[clustering]] == cluster]%>%as.matrix())
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


### read in data
obj = readRDS("~/Desktop/nemo.orig_harmony.integration_all_testd_clusters.rds")
deg_csvs <- dir('DEG Outputs/05_11_2025 Neg Bin w Doms New_clustering')
deg_path <- 'DEG Outputs/05_11_2025 Neg Bin w Doms New_clustering/'
dars = read.csv("Collaboration/all_clusters_DARs_peak_level_classified_with_support.csv")

# fix fragment file
# Individual-to-library mapping used for fragment assignment.

# read in DEGs list, I think first I will do all singular and nonsingular
all_degs <- data.frame()
for(i in deg_csvs){
  data <- read.csv(paste0(deg_path,i))
  all_degs <- rbind(all_degs, data)
}

all_degs_subset <- subset(all_degs, av_q.value <0.1 & f_m_p.value <0.05) 

all_dars = subset(dars, M_vs_F_FDR<0.05)

for(clusterx in unique(all_degs_subset$cluster)){
  deg_mean_data = data.frame()
  print(clusterx)
  
  # DEGs (RNA)
  subset_cluster = subset(all_degs_subset, cluster == clusterx)
  degs = unique(subset_cluster$gene)
  for(deg in degs){
    deg_data = mecd(obj,
                    deg,
                    clusterx, 
                    clustering = 'res0.8_50nn_40PC_45LSI')
    deg_data$gene = deg
    deg_mean_data = rbind(deg_mean_data, deg_data)
  }
  
  # DARs (ATAC) for this cluster
  subset_dars = subset(all_dars, cluster_id == clusterx)
  peaks = unique(subset_dars$peak)
  for(peak in peaks){
    peak_data = mecd_ATAC(obj,
                          peak,
                          clusterx,
                          clustering = 'res0.8_50nn_40PC_45LSI')
    peak_data$gene = peak
    deg_mean_data = rbind(deg_mean_data, peak_data)
  }
  
  pivoted_data = deg_mean_data%>%
    dplyr::select(-c(se, Sex))%>%
    pivot_wider(names_from = gene, values_from = mean)%>%
    filter(Status != 'NRM')
  
  path = 'Sex Classifier and Linearity/progress_classifier/Including DARs/degs_0.1_f_m_dars_0.05_f_m/'
  file_name = paste0('cluster_',clusterx,'.csv')
  write.csv(pivoted_data, paste0(path, file_name))
}
