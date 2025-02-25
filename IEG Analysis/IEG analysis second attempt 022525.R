#IEG analysis second attempt
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

obj$ieg <- ifelse(obj@assays$RNA$data['LOC111583367',] >0, 'ieg', NA)
obj$ieg <- ifelse(obj@assays$RNA$data['egr1',] >0, 'ieg', obj$ieg)
obj$ieg <- ifelse(obj@assays$RNA$data['npas4a',] >0, 'ieg', obj$ieg)
obj$ieg <- ifelse(is.na(obj$ieg), 'no_ieg', obj$ieg)

DimPlot(obj, split.by = 'ieg')

### filter out non-neuronal
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
DimPlot(neuronal_only, split.by = 'ieg')
#### ok a few questions I can now ask

#1) Which clusters are enriched for proportion of cells expressing ieg
#>2) How do the sexes differ in these measures
#>3) Do the findmarkers things zack discusses in his methods
#>4) 

ieg_data_cluster_level <- table(obj@meta.data$ieg, obj@meta.data$harmony.wnn_res0.4_clusters)%>%
  as.data.frame()%>%
  pivot_wider(names_from = 'Var1', values_from = 'Freq')%>%
  mutate(prop_pos = ieg/no_ieg)

###CLUSTER 27 IS THE HIGHEST WAOW
### ITS NOT EVEN CLOSE



