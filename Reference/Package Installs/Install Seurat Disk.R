
remotes::install_github("mojaveazure/seurat-disk", force =T)
remotes::install_github("satijalab/seurat-data", force =T)

library(Seurat)
library(SeuratData)
library(SeuratDisk)

rna_object = readRDS('/Users/ggraham/Desktop/snRNA-seq R Files 122524/RNA Object.rds')
rna_object@meta.data <- data.frame(rna_object@meta.data)  # Ensure it's a dataframe

DefaultAssay(rna_object) <- "RNA"
print(Assays(rna_object))  # Should contain "RNA"
print(colnames(rna_object@meta.data))

print(slotNames(rna_object@assays$RNA))  # Should include "counts" or "data"
rna_object <- JoinLayers(rna_object)
print(slotNames(rna_object@assays$RNA))  # Should include "counts" or "data"

rna_object@assays$RNA@counts = rna_object@assays$RNA$counts

#why the fuck does this not work and why are bioinformaticians such assholes and idiots god I am so tired of this shit

rna.v3 <- CreateAssayObject(counts = rna_object@assays$RNA$counts)
rna_object_v3 = CreateSeuratObject(rna.v3)
rna_object_v3@meta.data = rna_object@meta.data

SaveH5Seurat(rna_object_v3, filename = "/Users/ggraham/Desktop/snRNA-seq R Files 122524/RNA_object_anndata", overwrite = TRUE)

#it has to be a v3 object oh my christ satija lab should be blown up for doing this

h5seurat_obj <- LoadH5Seurat('/Users/ggraham/Desktop/snRNA-seq R Files 122524/RNA_object_anndata.h5seurat')


Convert('/Users/ggraham/Desktop/snRNA-seq R Files 122524/RNA_object_anndata.h5seurat', dest = "h5ad", overwrite = TRUE, verbose = TRUE)
