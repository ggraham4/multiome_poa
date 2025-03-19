### working on pseudobulk function
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


pseudobulk_function <- function(object = obj,
                                clustering = 'harmony.wnn_res0.4_clusters',
                                cluster = '19'){
  start_time <- Sys.time()  # Start timing
  
  #extract data matrix for cluster of interest
  data <- object@assays[["RNA"]]@layers[["data"]][,object@meta.data[[clustering]] == cluster]# here, we will be using the mean of nromalized data
  
  data_transposed <- t(data)%>%as.data.frame() #make columns genes

  data_transposed$individual <- object@meta.data$individual[object@meta.data[[clustering]] == cluster]
  ## calculate colmeans for each individual
  
  data_transposed_split <- split(data_transposed, f = data_transposed$individual) #split by subject
  
  gene_names <- rownames(object@assays$RNA$data)
  
  individuals <- unique(object$individual)
  
  meta.data <- data.frame(individual = object$individual,
                        Status = object$Status)%>%
    distinct()
  
  Statuses <- as.data.frame(individuals)%>%
    right_join(meta.data, by = join_by('individuals'=='individual'))%>%
    dplyr::select(Status)%>%
    as.vector()%>%
    unname()%>%
    unlist()
  
  message('calculating individual means')
  individual_means <- list() 
  for(i in individuals){
    #print(i)
  data <-   data_transposed_split[[i]]
  if(is.null(data)){next} #here, add which individuals were removed to a list, and remove those individuals from statuses and individuals
  # alternatively, make their values all NA
  individual_expression_means <- colMeans(data[,-ncol(data)])
  names(individual_expression_means)  <- gene_names 
   
  individual_means[[i]]$individual_expression_means = individual_expression_means
  }
  
  message('coalescing data')
  data_means <- data.frame() #bring them together into one df
  for(i in individuals){
    data <- individual_means[[i]]$individual_expression_means
    data_means <- rbind(data_means, data)
  }
  colnames(data_means) <- gene_names
  data_means$individual <- individuals
  data_means$Status <- Statuses
  
  #remove 0 columns
  columns_no_0 <- which(colSums(data_means[,1:length(gene_names)])!=0)
  individual_index <- which(colnames(data_means)=='individual')
  status_index <- which(colnames(data_means)=='Status')
  
  data_means_no_0 <- data_means[,c(columns_no_0, individual_index,status_index)]
  
  good_gene_names <- colnames(data_means[,c(columns_no_0)])
  
  index <- length(good_gene_names)
  
  message('fitting models')
 model_function <- function(i){
   message("Calculating gene ", i, ' of ', index)
   
   model_data <- data_means[,c(i, individual_index, status_index)]%>%
     subset(Status %in% c('M','D','F'))
   colnames(model_data) <- c('Expression', 'individual', 'Status')
   
   model <- lm(Expression ~ Status, data = model_data)
  
    av = anova(model, test ='Chisq')%>%as.data.frame()
    
    anova_p.value = av$`Pr(>F)`[1]
    new_data = data.frame(gene = good_gene_names[i],
                          cluster = cluster,
                          anova_p.value = anova_p.value)
    
    return(new_data)
   
 }

 results <- lapply(1:index, FUN = model_function)

 results_bound <- do.call(rbind, results)
 
 results_bound_nan_omit <- results_bound[!is.nan(results_bound$anova_p.value),]
 
 results_bound_nan_omit$anova_q.value <- p.adjust(results_bound_nan_omit$anova_p.value, 'fdr', nrow(results_bound_nan_omit))
 
 ### do pairwise comparisons for significant values
 
 results_bound_av_0.05 = results_bound_nan_omit[results_bound_nan_omit$anova_p.value<0.05,]
 
 significant_genes = results_bound_av_0.05$gene
 
 pairwise_function <- function(i){
   significant_gene_index = which(colnames(data_means)==i)
   
   message('Gene ', which(significant_genes==i), ' of ', length(significant_genes))
   model_data <- data_means[,c(significant_gene_index, individual_index, status_index)]%>%
     subset(Status %in% c('M','D','F'))
   colnames(model_data) <- c('Expression', 'individual', 'Status')
   
   model <- lm(Expression ~ Status, data = model_data)
   
  pairs = pairs(emmeans(model, "Status")) %>%as.data.frame()
  
  newd= data.frame(gene = i,
                   d_f_p.value = pairs$p.value[pairs$contrast=='D - F'],
                   d_m_p.value =  pairs$p.value[pairs$contrast=='D - M'],
                   f_m_p.value =  pairs$p.value[pairs$contrast=='F - M'],
                   d_f_estimate = pairs$estimate[pairs$contrast=='D - F'],
                   d_m_estimate = pairs$estimate[pairs$contrast=='D - M'],
                   f_m_estimate = pairs$estimate[pairs$contrast=='F - M']
                   
                   )
  return(newd)
  
   }
 pairwise_comparisons = lapply(c(significant_genes), FUN = pairwise_function)
 
 pairwise_comparisons_bound = do.call(rbind,pairwise_comparisons )
 
 complete_data <- results_bound_av_0.05%>%
   right_join(pairwise_comparisons_bound, by = 'gene')
 
 ### classify ####
 complete_data$class = NA
 
 complete_data$class[complete_data$d_m_p.value < 0.05 & 
              complete_data$d_m_estimate > 0 & 
              complete_data$d_f_p.value > 0.05] <- 'Early Upregulated'
 
 
 complete_data$class[complete_data$d_m_p.value > 0.05 & 
              complete_data$d_f_p.value < 0.05 & 
              complete_data$d_f_estimate < 0] <- 'Late Upregulated'
 
 complete_data$class[complete_data$d_m_p.value < 0.05 & 
              complete_data$d_m_estimate < 0 & 
              complete_data$d_f_p.value > 0.05] <- 'Early Downregulated'
 
 
 complete_data$class[complete_data$d_m_p.value > 0.05 & 
              complete_data$d_f_p.value < 0.05 & 
              complete_data$d_f_estimate > 0] <- 'Late Downregulated'
 
 complete_data$class[complete_data$d_m_p.value < 0.05 & 
              complete_data$d_f_p.value < 0.05 & 
              complete_data$d_f_estimate > 0 & 
              complete_data$d_m_estimate > 0] <- 'Transiently Upregulated'
 
 complete_data$class[complete_data$d_m_p.value < 0.05 & 
              complete_data$d_f_p.value < 0.05 & 
              complete_data$d_f_estimate < 0 & 
              complete_data$d_m_estimate < 0] <- 'Transiently Downregulated'
 
 complete_data$class[complete_data$d_m_p.value < 0.05 & 
              complete_data$d_f_p.value < 0.05 & 
              complete_data$d_m_estimate < 0 & 
              complete_data$d_f_estimate > 0] <- 'Progressively Downregulated'
 
 complete_data$class[complete_data$d_m_p.value < 0.05 & 
              complete_data$d_f_p.value < 0.05 & 
              complete_data$d_m_estimate > 0 & 
              complete_data$d_f_estimate < 0] <- 'Progressively Upregulated'
 
 complete_data$class[complete_data$f_m_p.value < 0.05 & 
              complete_data$d_f_p.value > 0.05 & 
              complete_data$d_m_p.value > 0.05 & 
              complete_data$f_m_estimate > 0 ] <- 'Terminally Upregulated'
 
 complete_data$class[complete_data$f_m_p.value < 0.05 & 
              complete_data$d_f_p.value > 0.05 & 
              complete_data$d_m_p.value > 0.05 & 
              complete_data$f_m_estimate < 0  ] <- 'Terminally Downregulated'
 
 complete_data$class <- ifelse(is.na(complete_data$class), 'Tukey > 0.05', complete_data$class)
 
 complete_data_final <- na.omit(complete_data)
 
 end_time <- Sys.time()  # End timing
 message(end_time - start_time)  # Print the time difference
 
 return(complete_data_final)
 
  
}

clust_0 = pseudobulk_function(object = obj,
                              cluster = 0)

for(i in 0:31){
  data = pseudobulk_function(object = obj,
                                cluster = i)
  assign(paste0('out_',i), data)
  
  
}

#BiocManager::install('qvalue')
library(qvalue)

out_2_no_na <- na.omit(out_2)
out_2_no_na$anova_q_diff <- qvalue(p=out_2_no_na$anova_p.value, lambda = 0)$qvalues # i think this only works on the whole list?

