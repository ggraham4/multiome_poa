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
behavior_deg_function <- function(object = obj, 
                                  cluster = 19,
                                  clustering = 'harmony.wnn_res0.4_clusters'){
  start <- Sys.time()
  library(lme4)
  library(dplyr)
  library(parallelsugar)
  
  options(dplyr.summarise.inform = FALSE)
  
  #extract counts matrix 
  
  ###here, i removed the m, d, f restriction because status is no longer a term in the model
  counts <- t(as.matrix(object@assays$RNA$data[, object@meta.data[[clustering]] == cluster]))
  
  
  
  filtered_cols_matrix <- counts[,which(colSums(counts) != 0)]
  
  meta_data <- data.frame(
    cells = rownames(object@meta.data[object@meta.data[[clustering]] == cluster ,]),
    Status = object@meta.data$Status[object@meta.data[[clustering]] == cluster ],
    individual =object@meta.data$individual[object@meta.data[[clustering]] == cluster ],
    Behaviors_Day_2 = object@meta.data$Behaviors_Day_2[object@meta.data[[clustering]] == cluster ],
    Time_Day_2 = object@meta.data$Time_Day_2[object@meta.data[[clustering]] == cluster ]
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
      subset(Status != 'NRM'&
               Status != 'NF'&
               Status != 'E')%>%
      summarize(expression = sum(gene))
    #I only want to restrict genes that are 0 in m, d, or f, the other sexes are too small to be useful
    
    if (any(joined_data_restrictions$expression < 1)) return(NULL)
    
    joined_data_for_model <- meta_data%>%
      group_by(individual)%>%
      summarize(expression = mean(gene),
                se_expression = sd(gene)/sqrt(n()),
                Behaviors_Day_2 = mean(Behaviors_Day_2),
                Time_Day_2 = mean(Time_Day_2))%>%
                  na.omit()
      
    if (length(unique(joined_data_for_model$expression)) == 1) return(NULL)
    
    ###Behavior model
    behavior_model <- lm(expression~Behaviors_Day_2, data=joined_data_for_model)
    
    #summary
    behavior_summary <- summary(behavior_model)
    behavior_coefs <- as.data.frame(behavior_summary$coefficients)
    
    #values
    behavior_estimate <- behavior_coefs[2,1]
    behavior_p.value <- behavior_coefs[2,4]

    data_for_output <- data.frame(cluster = cluster,
                                  gene = gene,
                                  behavior_p.value = behavior_p.value,
                                  behavior_estimate = behavior_estimate,
                                  behavior_p.value = behavior_p.value)
    return(data_for_output)
  }
  
  #deg_output <- parallelsugar::mclapply(X=genes, FUN=deg_function,mc.cores= detectCores()-1)
  deg_output <-lapply(X=genes, FUN=deg_function)
  
  deg_output2 <- do.call(rbind, deg_output)
  
  deg_output2$behavior_q.value <- p.adjust(deg_output2$behavior_p.value,'fdr', nrow(deg_output2)) 
  return(deg_output2)
  
  
}


for(i in 0:31){
  start <- Sys.time()
  print(i)
  behavior_deg_data <- behavior_deg_function(object = obj, 
                                             cluster = i,
                                             clustering = 'harmony.wnn_res0.4_clusters')
  assign(paste0('results_cluster_',i), behavior_deg_data, envir = .GlobalEnv)
  
  write.csv(behavior_deg_data,paste0('X:/Fish Lab/Experiments/sex change single nuc POA/Seurat Outputs/012525 Behavior DEGs no Status/behavior_results_cluster_',i,'.csv'))
  end <- Sys.time()
  print(end-start)
}















