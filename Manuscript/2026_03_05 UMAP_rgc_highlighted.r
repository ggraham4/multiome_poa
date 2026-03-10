library(Seurat)
library(ggplot2)

obj  = readRDS("~/Desktop/optimal_clustering_rna_only.rds")

obj$RGC = ifelse(obj@meta.data$final_clusters==1, 1, 0)

umap = FeaturePlot(obj, 'RGC', reduction = 'harmony_wnn.umap')+
  theme_void()+
  theme(legend.position = 'none', title = element_blank())
umap

ggsave(plot = umap,
       file = "UMAP_rgc_highlighted.tiff",
       device = "tiff",
       units = "in",
       width = 4,
       height = 4,
       path = "Manuscript/Plots/")

