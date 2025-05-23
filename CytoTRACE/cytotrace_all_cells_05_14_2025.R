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
  library(CytoTRACE)  
}

obj <- readRDS('C:/seurat_objects/optimal_clustering_05_06_2025.rds')

##cytotrace for all clusters first
counts_mat <- as.matrix(obj@assays$RNA$counts)

cyto = CytoTRACE(counts_mat,enableFast = T)

obj$cyto= cyto$CytoTRACE

gc()

ggplot(obj@meta.data, aes(x  = res0.8_50nn_40PC_45LSI, y = cyto))+
  geom_boxplot()

## 26 is mature huh? maybe needs to be glia only

ggplot(subset(obj@meta.data,res0.8_50nn_40PC_45LSI %in% c(1,11,26)) , aes(x  = res0.8_50nn_40PC_45LSI, y = cyto))+
  geom_boxplot()
### so maybe it is immature

obj$Status = factor(obj$Status, levels = c('M','D','E','NF','F'))
ggplot(subset(obj@meta.data,res0.8_50nn_40PC_45LSI %in% c(1)& !is.na(Status)) , aes(x  = Status, y = cyto))+
  geom_boxplot()

#saveRDS(obj, 'C:/seurat_objects/optimal_clustering_05_06_2025.rds')


