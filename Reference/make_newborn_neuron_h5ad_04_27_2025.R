library(Seurat)
library(SeuratData)
library(SeuratDisk)

rna_object = readRDS('C:/Users/Gabe/Desktop/RNA Object.rds')
rna_object@meta.data <- data.frame(rna_object@meta.data)  # Ensure it's a dataframe

radial_glia <- FindSubCluster(rna_object, cluster = 2, subcluster.name = 'sub', graph.name = 'harmony.wsnn')

radial_glia_subset <- subset(radial_glia, sub %in% c('2_9', '2_1','2_4','2_7') | harmony.wnn_res0.4_clusters%in% c(0,1,3:31))

DimPlot(radial_glia_subset, reduction = 'harmony_wnn.umap')
DimPlot(rna_object,  reduction = 'harmony_wnn.umap')

labels <- as.character(radial_glia_subset$harmony.wnn_res0.4_clusters)
locs_2_7 <- which(radial_glia_subset$sub == '2_7')
labels[locs_2_7] <- 'new_neurons'

radial_glia_subset$clusters = labels

`%notin%` <- Negate(`%in%`)
radial_glia_subset <- subset(radial_glia_subset, clusters %notin% c(4, 18, 14, 22, 29, 26))

DimPlot(radial_glia_subset,  reduction = 'harmony_wnn.umap', label = T, group.by = 'clusters')

DefaultAssay(radial_glia_subset) <- "RNA"
print(Assays(radial_glia_subset))  # Should contain "RNA"
print(colnames(radial_glia_subset@meta.data))

print(slotNames(radial_glia_subset@assays$RNA))  # Should include "counts" or "data"
print(slotNames(radial_glia_subset@assays$RNA))  # Should include "counts" or "data"

DefaultAssay(radial_glia_subset) <- "RNA"

rna.v3 <- CreateAssayObject(data = radial_glia_subset@assays$RNA$data)
rna_object_v3 = CreateSeuratObject(rna.v3)

# Copy over metadata
rna_object_v3@meta.data = radial_glia_subset@meta.data

# Copy over reductions - this is where the fix is needed
rna_object_v3@reductions = radial_glia_subset@reductions
rna_object_v3@assays$RNA$counts = radial_glia_subset@assays$RNA$counts

# Plot
DimPlot(rna_object_v3, reduction = 'harmony_wnn.umap', label = TRUE, group.by = 'clusters')

SaveH5Seurat(rna_object_v3, filename = "A:/support_vector_classifier_zeppilli_et_al_04_27_2025/data/RNA_object_new_neurons", overwrite = TRUE)

h5seurat_obj <- LoadH5Seurat('A:/support_vector_classifier_zeppilli_et_al_04_27_2025/data/RNA_object_new_neurons.h5seurat')

Convert('A:/support_vector_classifier_zeppilli_et_al_04_27_2025/data/RNA_object_new_neurons.h5seurat', dest = "h5ad", overwrite = TRUE, verbose = TRUE)

