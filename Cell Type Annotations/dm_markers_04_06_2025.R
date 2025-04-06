#Cell type annotations 121924
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
  #library(CytoTRACE)
  # SCTRransform_mean_plot <- readRDS("R/Gabe/SCTRransform_mean_plot.rds")
  #mac.neg.bin <- readRDS(file = 'R/Gabe/mac.neg.bin.rds')
  library('glmGamPoi')
  library(scran)
  library(parallel)
  library(factoextra)
  library(readxl)
  library(factoextra)
  library(forcats)
  library(ggrepel)
  library(biomaRt)
  #mean_cell <- readRDS('R/Gabe/mean_cell.rds')
  library(openxlsx)
  #clown_go <- readRDS('R/Gabe/clown_go.rds')
  library(ComplexHeatmap)
  
}

previous_object <- readRDS('C:/Users/Gabe/Desktop/Old object.rds')
Idents(previous_object) <- 'clusters49'
DimPlot(previous_object, label = T)

marks_dm <- FindMarkers(previous_object, '2')

marks_13d <- FindMarkers(previous_object, '25')

p1 <- FeaturePlot(previous_object, 'tshz1')
p2 <-
  FeaturePlot(previous_object, 'trarg1a')
p1+p2
