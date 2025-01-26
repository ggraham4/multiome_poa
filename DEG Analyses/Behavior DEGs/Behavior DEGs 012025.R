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
  counts <- t(as.matrix(object@assays$RNA$data[, object@meta.data[[clustering]] == cluster & (object@meta.data$Status == "M" | object@meta.data$Status == "F" | object@meta.data$Status == "D")]))


  
  filtered_cols_matrix <- counts[,which(colSums(counts) != 0)]

  meta_data <- data.frame(
    cells = rownames(object@meta.data[object@meta.data[[clustering]] == cluster & (object@meta.data$Status == "M" | object@meta.data$Status == "F" | object@meta.data$Status == "D"),]),
    Status = object@meta.data$Status[object@meta.data[[clustering]] == cluster & (object@meta.data$Status == "M" | object@meta.data$Status == "F" | object@meta.data$Status == "D")],
    individual =object@meta.data$individual[object@meta.data[[clustering]] == cluster & (object@meta.data$Status == "M" | object@meta.data$Status == "F" | object@meta.data$Status == "D")],
    Behaviors_Day_2 = object@meta.data$Behaviors_Day_2[object@meta.data[[clustering]] == cluster & (object@meta.data$Status == "M" | object@meta.data$Status == "F" | object@meta.data$Status == "D")],
    Time_Day_2 = object@meta.data$Time_Day_2[object@meta.data[[clustering]] == cluster & (object@meta.data$Status == "M" | object@meta.data$Status == "F" | object@meta.data$Status == "D")]
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
                Behaviors_Day_2 = mean(Behaviors_Day_2),
                Time_Day_2 = mean(Time_Day_2)
      )
    if (length(unique(joined_data_for_model$expression)) == 1) return(NULL)
    
    ###Behavior model
  behavior_model <- lm(expression~Behaviors_Day_2*Status, data=joined_data_for_model)
  
  #summary
  behavior_summary <- summary(behavior_model)
  behavior_coefs <- as.data.frame(behavior_summary$coefficients)
  
  #aov
  behavior_aov <- as.data.frame(anova(behavior_model, test= 'Chisq'))
  behavior_aov$term <- rownames(behavior_aov)
  
  #pwise comparisons
  behavior_pairs <- suppressMessages(as.data.frame(pairs(emmeans(behavior_model, 'Status'), adjust = 'none')))
  
  #values
  behavior_p.value <- behavior_aov$`Pr(>F)`[behavior_aov$term == 'Behaviors_Day_2']
  behavior_estimate <- behavior_coefs[2,1]
  status_p.value <- behavior_aov$`Pr(>F)`[behavior_aov$term == 'Status']
  behavior_x_status_p.value <-behavior_aov$`Pr(>F)`[behavior_aov$term == 'Behaviors_Day_2:Status']
  
  d_f_p.value <- behavior_pairs$p.value[behavior_pairs$contrast == 'D - F']
  d_m_p.value <- behavior_pairs$p.value[behavior_pairs$contrast == 'D - M']
  f_m_p.value <- behavior_pairs$p.value[behavior_pairs$contrast == 'F - M']
  
  d_f_estimate <- behavior_pairs$estimate[behavior_pairs$contrast == 'D - F']
  d_m_estimate <- behavior_pairs$estimate[behavior_pairs$contrast == 'D - M']
  f_m_estimate <- behavior_pairs$estimate[behavior_pairs$contrast == 'F - M']
  
  data_for_output <- data.frame(cluster = cluster,
                            gene = gene,
                            behavior_p.value = behavior_p.value,
                            behavior_estimate = behavior_estimate,
                            behavior_x_status_p.value = behavior_x_status_p.value,
                            status_p.value = status_p.value,
                            d_f_p.value = d_f_p.value,
                            d_f_estimate=d_f_estimate,
                            d_m_p.value =d_m_p.value,
                            d_m_estimate = d_m_estimate,
                            f_m_p.value = f_m_p.value,
                            f_m_estimate=f_m_estimate)
  return(data_for_output)
  }
  
  deg_output <- parallelsugar::mclapply(X=genes, FUN=deg_function,mc.cores= detectCores()-1)

  deg_output2 <- do.call(rbind, deg_output)
  
  deg_output2$behavior_q.value <- p.adjust(deg_output2$behavior_p.value,'fdr', nrow(deg_output2)) 
    deg_output2$status_q.value <- p.adjust(deg_output2$status_p.value,'fdr', nrow(deg_output2)) 
  deg_output2$behavior_x_status_q.value <-p.adjust(deg_output2$behavior_x_status_p.value,'fdr', nrow(deg_output2)) 
  
  deg_output2$d_f_q.value <- p.adjust(deg_output2$d_f_p.value,'fdr', nrow(deg_output2)) 
  deg_output2$d_m_q.value <- p.adjust(deg_output2$d_m_p.value,'fdr', nrow(deg_output2)) 
  deg_output2$f_m_q.value <- p.adjust(deg_output2$f_m_p.value,'fdr', nrow(deg_output2))
  
  return(deg_output2)
                            
                          
    }

  
for(i in 0:31){
  start <- Sys.time()
  print(i)
  behavior_deg_data <- behavior_deg_function(object = obj, 
                                             cluster = i,
                                             clustering = 'harmony.wnn_res0.4_clusters')
  assign(paste0('results_cluster_',i), behavior_deg_data, envir = .GlobalEnv)
  
  write.csv(behavior_deg_data,paste0('X:/Fish Lab/Experiments/sex change single nuc POA/Seurat Outputs/012025 Behavior DEGs/behavior_results_cluster_',i,'.csv'))
  end <- Sys.time()
  print(end-start)
  }















