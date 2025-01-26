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
Test_deg_function <- function(object = obj, 
                              cluster = 19,
                              clustering = 'harmony.wnn_res0.4_clusters'){
  start <- Sys.time()
  library(lme4)
  library(dplyr)
  library(parallelsugar)
  
  options(dplyr.summarise.inform = FALSE)
  
  #extract counts matrix 
  counts <- t(as.matrix(object@assays$RNA$data[, object@meta.data[[clustering]] == cluster ]))
  
  
  
  filtered_cols_matrix <- counts[,which(colSums(counts) != 0)]
  
  meta_data <- data.frame(
    cells = rownames(object@meta.data[object@meta.data[[clustering]] == cluster ,]),
    Status = object@meta.data$Status[object@meta.data[[clustering]] == cluster ],
    individual =object@meta.data$individual[object@meta.data[[clustering]] == cluster ],
    Percent_Testicular = object@meta.data$Percent_[object@meta.data[[clustering]] == cluster ],
    Percent_Testicular = object@meta.data$Percent_Testicular[object@meta.data[[clustering]] == cluster ]
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
                Percent_Testicular = mean(Percent_Testicular),
                Percent_Testicular = mean(Percent_Testicular)
      )%>%na.omit()
    
    if (length(unique(joined_data_for_model$expression)) == 1) return(NULL)
    
    ###Test model
    Test_model <- lm(expression~Percent_Testicular, data=joined_data_for_model)
    
    #summary
    Test_summary <- summary(Test_model)
    Test_coefs <- as.data.frame(Test_summary$coefficients)

    #values 
    percent_testicular_estimate <- Test_coefs[2,1]
    percent_testicular_summary_p.value <- Test_coefs[2,4]


    
    data_for_output <- data.frame(cluster = cluster,
                                  gene = gene,
                                  percent_testicular_estimate=percent_testicular_estimate,
                                  percent_testicular_summary_p.value=percent_testicular_summary_p.value
    )
    return(data_for_output)
  }
  
  #deg_output <- parallelsugar::mclapply(X=genes, FUN=deg_function,mc.cores= detectCores()-1)
  deg_output <- lapply(X=genes, FUN=deg_function)
  
  deg_output2 <- do.call(rbind, deg_output)
  
  deg_output2$percent_testicular_summary_q.value <- p.adjust(deg_output2$percent_testicular_summary_p.value, 'fdr', nrow(deg_output2))

  
  
  return(deg_output2)
  
  
}


for(i in 0:31){
  start <- Sys.time()
  print(i)
  
  test_deg_data <- Test_deg_function(object = obj, 
                                     cluster = i,
                                     clustering = 'harmony.wnn_res0.4_clusters')
  assign(paste0('test_results_cluster_',i), test_deg_data, envir = .GlobalEnv)
  
  write.csv(test_deg_data,paste0('X:/Fish Lab/Experiments/sex change single nuc POA/Seurat Outputs/012625 Testicular DEGs no Status/Test_results_cluster_',i,'.csv'))
  end <- Sys.time()
  print(end-start)
}















