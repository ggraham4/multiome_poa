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

}

obj <-  readRDS("~/Desktop/snRNA-seq R Files 122524/RNA Object.rds")


### For ANY implementation of mean_expression_cluster_data, this line needs to precede it
# counts <- t(object@assays$RNA$data[,object@meta.data[[clustering]] == cluster])

change_Mass_function <- function(object = obj, 
                              cluster = 19,
                              clustering = 'harmony.wnn_res0.4_clusters'){
  start <- Sys.time()
  library(lme4)
  library(dplyr)
  library(parallel)
  
  options(dplyr.summarise.inform = FALSE)
  
  #extract counts matrix 
  counts <- t(as.matrix(object@assays$RNA$data[, object@meta.data[[clustering]] == cluster & object@meta.data$Status == 'D']))
  
  filtered_cols_matrix <- counts[,which(colSums(counts) != 0)]
  
  meta_data <- data.frame(
    cells = rownames(object@meta.data[object@meta.data[[clustering]] == cluster & object@meta.data$Status == 'D',]),
    Status = object@meta.data$Status[object@meta.data[[clustering]] == cluster & object@meta.data$Status == 'D'],
    individual =object@meta.data$individual[object@meta.data[[clustering]] == cluster & object@meta.data$Status == 'D'],
    change_Mass = object@meta.data$Change_Mass[object@meta.data[[clustering]] == cluster &object@meta.data$Status == 'D']
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
                change_Mass = mean(change_Mass)
      )%>%na.omit()
    
    if (length(unique(joined_data_for_model$expression)) == 1) return(NULL)
    
    ###Mass model
    Mass_model <- lm(expression~change_Mass, data=joined_data_for_model)
    
    #summary
    Mass_summary <- summary(Mass_model)
    Mass_coefs <- as.data.frame(Mass_summary$coefficients)
    
    #values 
    Mass_estimate <- Mass_coefs[2,1]
    Mass_summary_p.value <- Mass_coefs[2,4]

    
    data_for_output <- data.frame(cluster = cluster,
                                  gene = gene,
                                  Mass_estimate=Mass_estimate,
                                  Mass_summary_p.value=Mass_summary_p.value
    )
    return(data_for_output)
  }
  
  #deg_output <- parallelsugar::mclapply(X=genes, FUN=deg_function,mc.cores= detectCores()-1)
  deg_output <- mclapply(X=genes, FUN=deg_function
                       , mc.cores= detectCores()-1
                       )
  
  deg_output2 <- do.call(rbind, deg_output)
  
  deg_output2$Mass_summary_q.value <- p.adjust(deg_output2$Mass_summary_p.value, 'fdr', nrow(deg_output2))

  
  
  return(deg_output2)
  
}

for(i in 31:0){
  start <- Sys.time()
  print(i)
  
  change_Mass_data <- change_Mass_function(object = obj, 
                                     cluster = i,
                                     clustering = 'harmony.wnn_res0.4_clusters')
  assign(paste0('change_Mass_results_cluster_',i), change_Mass_data, envir = .GlobalEnv)
  
  write.csv(change_Mass_data ,paste0('DEG Outputs/022625 Change Mass DEGs/results_cluster_',i,'.csv'))
  end <- Sys.time()
  print(end-start)
}

for(i in 0:31){
  print(i)
  data <- suppressMessages(read_csv(paste0('DEG Outputs/022625 Change Mass DEGs/results_cluster_',i,'.csv')))
  print(data$gene[data$Mass_summary_q.value<0.05])
  #none, what happens if i include expandeds
  }

