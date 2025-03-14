#### I want to integrate all of the data we have using mouse gene names and harmony to publicly avlailable data
# that same way that laurent guy and maria tosches do

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

multiome_object_mouse_names <- readRDS("C:/Users/Gabe/Desktop/RNA object mouse names.rds")
DimPlot(multiome_object_mouse_names)
rownames(multiome_object_mouse_names)%>%head()

hypomap_data <- readRDS("A:/external_data_03_14_2025/hypoMap.rds")
DimPlot(hypomap_data)
rownames(hypomap_data)%>%head()


parker_et_al_mouse_names <- readRDS("A:/Anemonefish POA Legacy R Objects/GT_ortho_v07.rds")
rownames(parker_et_al_mouse_names)%>%head()
#there is no dimensional reduction in this object so I will run that quickly

parker_et_al_mouse_names <- NormalizeData(parker_et_al_mouse_names, assay = 'RNA')
parker_et_al_mouse_names <- FindVariableFeatures(parker_et_al_mouse_names, nfeatures = 4000)
parker_et_al_mouse_names <- ScaleData(parker_et_al_mouse_names)

parker_et_al_mouse_names <- RunPCA(parker_et_al_mouse_names, 
                   npcs = 50, 
                   features = rownames(parker_et_al_mouse_names[['RNA']]))

ElbowPlot(parker_et_al_mouse_names)

parker_et_al_mouse_names <- FindNeighbors(parker_et_al_mouse_names, dims = 1:7)

parker_et_al_mouse_names <- FindClusters(parker_et_al_mouse_names, resolution = 0.5)

parker_et_al_mouse_names <- RunUMAP(parker_et_al_mouse_names, dims = 1:7)
  
DimPlot(parker_et_al_mouse_names)
FeaturePlot(parker_et_al_mouse_names, 'Cyp19a1')



#### integrate them with harmony
parker_et_al_mouse_names$source <- 'parker_et_al_mouse_names'
multiome_object_mouse_names$source <- 'multiome_object_mouse_names'
hypomap_data$source <- 'hypomap'

parker_et_al_mouse_names@reductions = list()
multiome_object_mouse_names@reductions <- list()
hypomap_data@reductions <- list()

merged <- merge(hypomap_data, 
                y = c(multiome_object_mouse_names, parker_et_al_mouse_names), 
                add.cell.ids = c('hypo', 'multi', 'parker'))

merged <- NormalizeData(merged)
merged <- FindVariableFeatures(merged, selection.method = "vst", nfeatures = 2000)
merged <- ScaleData(merged)
merged <- RunPCA(merged, npcs = 50)

library(harmony)
merged <- RunHarmony(merged, group.by.vars = "source", dims.use = 1:30)

# Run UMAP on harmony embeddings
merged <- RunUMAP(merged, reduction = "harmony", dims = 1:30)

# Clustering
merged <- FindNeighbors(merged, reduction = "harmony", dims = 1:30)
merged <- FindClusters(merged, resolution = 0.5)

# Visualization
merged$source <- ifelse(!is.na(merged$Author_Class), 'hypomap_data',merged$source)
DimPlot(merged, reduction = "umap", group.by = "source")

#saveRDS(merged, 'A:/mouse_clownfish_merged_harmony.rds')

DimPlot(merged, group.by = 'clusters20', label = T)
DimPlot(merged, group.by = 'harmony.wnn_res0.4_clusters', label = T)
#why is the multiome object still completely excluded

merged <- SCTransform(merged)
