library(Seurat)
library(SeuratData)
library(SeuratDisk)


obj = readRDS("C:/Users/Gabe/Desktop/nemo.orig_harmony.integration_all_testd_clusters.rds")
Idents(obj) = 'res0.8_50nn_40PC_45LSI'
obj = FindSubCluster(obj, 1, graph.name='harmony.wsnn')
Idents(obj) = 'sub.cluster'
DimPlot(obj, reduction = 'harmony_wnn.umap')

obj@meta.data <- data.frame(obj@meta.data)  # Ensure it's a dataframe
DefaultAssay(obj) <- "RNA"

rna.v3 <- CreateAssayObject( 
  data = obj@assays$RNA$data)
  
obj_v3 = CreateSeuratObject(rna.v3)

obj_v3@meta.data = obj@meta.data

SaveH5Seurat(obj_v3, filename = "A:/2025_10_19_obj_rgc_subclusters", overwrite = TRUE)
h5seurat_obj <- LoadH5Seurat('A:/2025_10_19_obj_rgc_subclusters.h5seurat')
Convert('A:/2025_10_19_obj_rgc_subclusters.h5seurat', dest = "h5ad", overwrite = TRUE, verbose = TRUE)

