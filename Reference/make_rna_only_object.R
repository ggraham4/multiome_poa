
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

}

multiome_object <- readRDS("~/Desktop/nemo.orig_harmony.integration_all_testd_clusters.rds")

DefaultAssay(multiome_object) <- 'RNA'


multiome_object <- DietSeurat(
  multiome_object,
  assays = "RNA",
  dimreducs = names(multiome_object@reductions)[6], 
  graphs = 'harmony.wsnn', # Optionally include or exclude graphs if present
  misc = TRUE # Retain misc if applicable
  )

Idents(multiome_object) <- "res0.8_50nn_40PC_45LSI"
DimPlot(multiome_object, label = T)

saveRDS(multiome_object, '~/Desktop/optimal_clustering_rna_only.rds')
multiome_object <- readRDS('~/Desktop/optimal_clustering_rna_only.rds')

FeaturePlot(multiome_object, 'crocc2')
