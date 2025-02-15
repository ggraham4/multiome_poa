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
}

obj <- readRDS("C:/Users/Gabe/Desktop/RNA Object.rds")
mean_expression_cluster_data<- readRDS('Functions/mean_expression_cluster_data.rds')


### Read in DEG data ###
degs_list<- list()
for(i in 0:31){
  print(i)
  data <- read.csv(
    paste0('DEG Outputs/012425 Neg Bin w Doms Lower Stringency/cluster_',i,'.csv')
    )
  degs<- data$gene[(data$f_m_q.value<0.1|
                   data$d_m_q.value<0.1|
                   data$d_f_q.value<0.1)]
  degs <- degs[!is.na(degs)]
  degs <- c(degs)
  
  degs_list[[paste0(i)]] = c(degs)
  
}

### Summarize Gene expression of DEGs ### 
expression_data <- list()
for(i in 0:31){
  print(i)
  cluster = paste0(i)
  cluster_data <- data.frame()
  for(g in degs_list[[cluster]]){
    gene_expression_data <- mean_expression_cluster_data(object = obj,
                                                         gene = g,
                                                         cluster = i)
    gene_expression_data$gene = g
    cluster_data <- rbind(cluster_data, gene_expression_data)
  
    }
  expression_data[[paste0(i)]] = cluster_data
  }




