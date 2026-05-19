library(Seurat)
library(SeuratObject)
library(Matrix)
library(ggplot2)

obj        <- readRDS("~/Desktop/optimal_clustering_rna_only.rds")
Idents(obj) <- "final_clusters"

vis_neuron <- readRDS('/Users/ggraham/Desktop/Visium/vis_neuron.rds')


DefaultAssay(obj)        <- "RNA"
DefaultAssay(vis_neuron) <- "Spatial.Polygons"

anchors <- FindTransferAnchors(
  reference            = vis_neuron,
  query                = obj,
  normalization.method = "LogNormalize",
  dims                 = 1:10
)

obj <- TransferData(
  anchorset        = anchors,
  query            = obj,
  reference        = vis_neuron,
  refdata          = list(anatomical = "anatomical_renamed"),
  dims             = 1:10
)

# obj$predicted.anatomical       — winning region per cell
# obj$predicted.anatomical.score — confidence

# plot on multiome UMAP
DimPlot(obj, 
        group.by  = "predicted.anatomical",
        label     = TRUE,
        reduction = "harmony_wnn.umap") +
  labs(title = "Multiome cells labeled by predicted Visium anatomical region")

# optionally mask low confidence calls
obj$predicted.anatomical.confident <- ifelse(
  obj$predicted.anatomical.score > 0.2,
  obj$predicted.anatomical,
  "low confidence"
)

DimPlot(obj,
        group.by  = "predicted.anatomical.confident",
        label     = TRUE,
        reduction = "harmony_wnn.umap")

library(tidyverse)
clusts = table(obj$final_clusters, 
      obj$predicted.anatomical.confident)%>%as.data.frame.matrix()
clusts_norm =clusts/rowSums(clusts)

library(pheatmap)
pheatmap(clusts_norm)
