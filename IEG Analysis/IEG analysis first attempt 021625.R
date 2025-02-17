### Attempting IEG analysis

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
}

obj <-  readRDS("~/Desktop/snRNA-seq R Files 122524/RNA Object.rds")

#ok the first thing I want to do is in each neuronal cluster identify genes coexpressed with
#(egr1, npas4a, and fosab)

#create subset object only including neuronal
neuronal_only <- subset(obj, harmony.wnn_res0.4_clusters %notin%
                          c(
                            2, #radial glia
                            4, #olig
                            14, #mg
                            18, #OPC
                            22, #ependymal
                            26, #leuko,
                            29, #fibro
                            15, #exclude
                            30 #OB
                          )
)
DimPlot(neuronal_only, reduction = 'harmony_wnn.umap', label = T)

#OK now for each cluster do correlations with each gene

#### Testing
  mean_expression_data_transposed <- as.data.frame(t(mean_expression_data$RNA))
  
  expression_data <- as.data.frame(t(neuronal_only@assays$RNA$data[,neuronal_only@meta.data$harmony.wnn_res0.4_clusters==19]))
  expression_data_filtered <- expression_data[,which(colSums(expression_data)!=0)]
  
  indices_to_remove <- c()
  for(i in 1:ncol(expression_data_filtered)){
    print(i)
    data_to_test <- expression_data_filtered[,i]
    if(length(unique(data_to_test))==1){
    indices_to_remove <- c(indices_to_remove, i)
      }
    next
    }
  #ok none
  
  #test correlations for each gene
  #LOC111583367 = c-fos like
  #LOC111583368 = c-fos like, also called fosab
  c_fos_index <- which(colnames(expression_data_filtered)=='LOC111583367')
  c_fos_data <- expression_data_filtered[,c_fos_index]
  
  #dont want to test for correlations with itself
  expression_data_filtered <- expression_data_filtered[,-c_fos_index]
  
  ### npas4 
  npas4a_index <-  which(colnames(expression_data_filtered)=='npas4a')
  npas4a_data <- expression_data_filtered[,npas4a_index]
  
  expression_data_filtered <- expression_data_filtered[,-npas4a_index]

    egr1_index <-  which(colnames(expression_data_filtered)=='egr1')
  egr1_data <- expression_data_filtered[,egr1_index]

    expression_data_filtered <- expression_data_filtered[,-egr1_index]
    
    #initialize list
    c_fos_corrs <-data.frame()
    bulk_corr_function <- function(gene ){
      message(
        paste0(
          which(colnames(expression_data_filtered)==gene),
          ' of ',
          length(colnames(expression_data_filtered))))
      
      query_data <- expression_data_filtered[,which(colnames(expression_data_filtered)==gene)]
      
      corr_test <- cor.test(c_fos_data, query_data)
      
      result_data <- data.frame(gene = gene,
                                ieg = 'cfos',
                                estimate = corr_test$estimate,
                                p.value = corr_test$p.value)
      c_fos_corrs <- rbind(result_data, c_fos_corrs)

    }
    test_cfos_1 <- lapply(X =  colnames(expression_data_filtered), FUN = bulk_corr_function)
    test_cfos_1_bound <- do.call(rbind, test_cfos_1)
    test_cfos_1_bound$q.value <- p.adjust(test_cfos_1_bound$p.value, 'fdr',nrow(test_cfos_1_bound))
    test_cfos_1_bound$bonf.p.value <- p.adjust(test_cfos_1_bound$p.value, 'bonferroni',nrow(test_cfos_1_bound))

    ieg_like_genes_cfos <-test_cfos_1_bound$gene[test_cfos_1_bound$q.value<0.05]
    


