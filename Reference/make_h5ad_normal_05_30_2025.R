library(Seurat)
library(SeuratData)
library(SeuratDisk)

rna_object  = readRDS("C:/seurat_objects/optimal_clustering_05_06_2025.rds")

DefaultAssay(rna_object) <- "RNA"

rna_object$final_clusters = rna_object$res0.8_50nn_40PC_45LSI

Idents(rna_object) = 'final_clusters'

DimPlot(rna_object)

saveRDS(rna_object, "C:/seurat_objects/optimal_clustering_05_06_2025.rds")

rna.v3 <- CreateAssayObject( 
  data = rna_object@assays$RNA$data)
rna_object_v3 = CreateSeuratObject(rna.v3)
rna_object_v3@meta.data = rna_object@meta.data

SaveH5Seurat(rna_object_v3, filename = 'C:/seurat_objects/optimal_clustering_05_06_2025.h5seurat', overwrite = TRUE)

#it has to be a v3 object oh my christ satija lab should be blown up for doing this

h5seurat_obj <- LoadH5Seurat('C:/seurat_objects/optimal_clustering_05_06_2025.h5seurat')

Convert('C:/seurat_objects/optimal_clustering_05_06_2025.h5seurat', dest = "h5ad", overwrite = TRUE, verbose = TRUE)
