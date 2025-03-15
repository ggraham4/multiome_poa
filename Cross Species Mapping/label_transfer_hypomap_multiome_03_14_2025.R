#label transfer hypomap multiome

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
multiome_object_mouse_names <- readRDS("C:/Users/Gabe/Desktop/RNA object mouse names.rds")
multiome_object_mouse_names <- FindVariableFeatures(multiome_object_mouse_names)
Idents(multiome_object_mouse_names) <- "harmony.wnn_res0.4_clusters"

hypomap_data <- readRDS("A:/external_data_03_14_2025/hypoMap.rds")
hypomap_data <- FindVariableFeatures(hypomap_data)

### query multiome

#set up reductions
Idents(hypomap_data) <- 'Region_predicted'
DimPlot(hypomap_data, label = T)+ theme(legend.position = 'none')
hypomap_data <- ScaleData(hypomap_data)
hypomap_data <- RunPCA(hypomap_data, 
                                   npcs = 50, 
                                   features = rownames(hypomap_data[['RNA']]))

anchors <- FindTransferAnchors(reference = hypomap_data, query = multiome_object_mouse_names,
                               reference.reduction = "pca")

predictions <- TransferData(anchorset = anchors, refdata = hypomap_data$Region_predicted)
multiome_object_mouse_names <- AddMetaData(multiome_object_mouse_names, metadata = predictions)

Idents(multiome_object_mouse_names) <- 'predicted.id'
DimPlot(multiome_object_mouse_names, label = T)+ theme(legend.position = 'none')
#aw man lmfaooooo



###### query hypomap
Idents(multiome_object_mouse_names) <- 'harmony.wnn_res0.4_clusters'
anchors2 <- FindTransferAnchors(reference = multiome_object_mouse_names, query = hypomap_data,
                                reference.reduction = "pca")

predictions2 <- TransferData(anchorset = anchors2,
                             refdata = multiome_object_mouse_names$harmony.wnn_res0.4_clusters)
hypomap_data <- AddMetaData(hypomap_data, metadata = predictions2)
###it's just not gonna work is it

Idents(hypomap_data) <- 'predicted.id'
DimPlot(hypomap_data, label = T)

### maybe try integration?? ###

merged <- readRDS('A:/mouse_clownfish_merged_harmony.rds')
library(harmony)
merged <- IntegrateLayers(object = merged, method = HarmonyIntegration, orig.reduction = "pca",
                          new.reduction = "integrated.harmony",
                        verbose = FALSE)

merged[["RNA"]] <- JoinLayers(merged[["RNA"]])

merged <- FindNeighbors(merged, reduction = "integrated.harmony", dims = 1:30)
merged <- FindClusters(merged, resolution = 1)

merged <- RunUMAP(merged, dims = 1:30, reduction = "integrated.harmony")

DimPlot(merged, reduction = "umap", group.by = c("source"))
DimPlot(merged, reduction = "umap", group.by = c("harmony.wnn_res0.4_clusters"))
DimPlot(merged, reduction = "umap", group.by = c("clusters49"))






