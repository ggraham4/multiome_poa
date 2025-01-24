#Testing reading in seurat object

library(Seurat)

multiome_obj <- readRDS("~/Desktop/multiome_poa/multiome_poa/RNA Object.rds")

DimPlot(multiome_obj)
