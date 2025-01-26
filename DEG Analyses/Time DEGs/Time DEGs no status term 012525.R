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
  library(emmeans)
  library(CytoTRACE)
  library(ggrepel)
  
}

obj <- readRDS('C:/Users/Gabe/Desktop/RNA Object.rds')
Time_deg_function <- function(object = obj, 
                              cluster = 19,
                              clustering = 'harmony.wnn_res0.4_clusters'){
  start <- Sys.time()
  library(lme4)
  library(dplyr)
  library(parallelsugar)
  
  options(dplyr.summarise.inform = FALSE)
  
  #extract counts matrix 
  counts <- t(as.matrix(object@assays$RNA$data[, object@meta.data[[clustering]] == cluster]))
  
  
  
  filtered_cols_matrix <- counts[,which(colSums(counts) != 0)]
  
  meta_data <- data.frame(
    cells = rownames(object@meta.data[object@meta.data[[clustering]] == cluster,]),
    Status = object@meta.data$Status[object@meta.data[[clustering]] == cluster],
    individual =object@meta.data$individual[object@meta.data[[clustering]] == cluster],
    Time_Day_2 = object@meta.data$Time_Day_2[object@meta.data[[clustering]] == cluster],
    Time_Day_2 = object@meta.data$Time_Day_2[object@meta.data[[clustering]] == cluster]
  )
  
  genes <- colnames(filtered_cols_matrix)
  genes <- genes[!is.na(genes)]
  
  n_genes <- length(genes)
  
  
  deg_function <- function(gene){
    
    message(paste0('gene ', which(genes == gene), ' of ', n_genes))
    
    gene_expression <- filtered_cols_matrix[, gene]
    meta_data$gene <- gene_expression
    
    joined_data_restrictions <- meta_data%>%
      group_by(Status)%>%
      summarize(expression = sum(gene))
    
    if (any(joined_data_restrictions$expression < 1)) return(NULL)
    
    joined_data_for_model <- meta_data%>%
      group_by(Status, individual)%>%
      summarize(expression = mean(gene),
                se_expression = sd(gene)/sqrt(n()),
                Time_Day_2 = mean(Time_Day_2),
      )
    if (length(unique(joined_data_for_model$expression)) == 1) return(NULL)
    
    ###Time model
    Time_model <- lm(expression~Time_Day_2, data=joined_data_for_model)
    
    #summary
    Time_summary <- summary(Time_model)
    Time_coefs <- as.data.frame(Time_summary$coefficients)
    
    #values
    Time_estimate <- Time_coefs[2,1]
    Time_p.value <- Time_coefs[2,4]
    

    data_for_output <- data.frame(cluster = cluster,
                                  gene = gene,
                                  Time_estimate=Time_estimate,
                                  Time_p.value= Time_p.value)
    return(data_for_output)
  }
  
  #deg_output <- parallelsugar::mclapply(X=genes, FUN=deg_function,mc.cores= detectCores()-1)
  deg_output <-lapply(X=genes, FUN=deg_function)
  
  deg_output2 <- do.call(rbind, deg_output)
  
  deg_output2$Time_q.value <- p.adjust(deg_output2$Time_p.value,'fdr', nrow(deg_output2)) 

  return(deg_output2)
  
  
}


for(i in 0:31){
  start <- Sys.time()
  print(i)
  Time_deg_data <- Time_deg_function(object = obj, 
                                     cluster = i,
                                     clustering = 'harmony.wnn_res0.4_clusters')
  assign(paste0('time_results_cluster_',i), Time_deg_data, envir = .GlobalEnv)
  
  write.csv(Time_deg_data,paste0('X:/Fish Lab/Experiments/sex change single nuc POA/Seurat Outputs/012525 Time DEGs no Status/Time_results_cluster_',i,'.csv'))
  end <- Sys.time()
  print(end-start)
}















