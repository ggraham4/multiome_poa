library(Seurat)
library(SeuratData)
library(SeuratDisk)

rna_object = readRDS('C:/Users/Gabe/Desktop/RNA_object_human_names.rds')

DefaultAssay(rna_object) <- "RNA"

rna.v3 <- CreateAssayObject( 
  data = rna_object@assays$RNA$data)
rna_object_v3 = CreateSeuratObject(rna.v3)
rna_object_v3@meta.data = rna_object@meta.data

rna_object_v3$harmony.wnn_res0.4_clusters = as.character(rna_object_v3$harmony.wnn_res0.4_clusters)

SaveH5Seurat(rna_object_v3, filename = "C:/Users/Gabe/Desktop/RNA_object_human_names_anndata", overwrite = TRUE)

#it has to be a v3 object oh my christ satija lab should be blown up for doing this

h5seurat_obj <- LoadH5Seurat('C:/Users/Gabe/Desktop/RNA_object_human_names_anndata.h5seurat')

Convert('C:/Users/Gabe/Desktop/RNA_object_human_names_anndata.h5seurat', dest = "h5ad", overwrite = TRUE, verbose = TRUE)
