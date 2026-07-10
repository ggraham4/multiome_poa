library(Seurat)
library(zellkonverter)   # BiocManager::install("zellkonverter")
library(SingleCellExperiment)

# ── Multiome → h5ad ────────────────────────────────────────────────────────
multiome <- readRDS("C:/Users/Gabe/Desktop/nemo.orig_harmony.integration_all_testd_clusters.rds")
DefaultAssay(multiome) <- "RNA"

### subset to rgc
multiome = FindSubCluster(multiome,
                          6, 'harmony.wsnn')

sub_6 = subset(multiome, res0.8_50nn_40PC_45LSI%in%c(6))

Idents(sub_6) = 'sub.cluster'
poa = DimPlot(sub_6, label = F, reduction = 'harmony_wnn.umap')+
  theme_void()+
  theme(legend.position = 'none', title = element_blank())

poa

ggsave(plot = poa,
       file = "UMAP_poa_subclustered.tiff",
       device = "tiff",
       units = "in",
       width = 4,
       height = 4,
       path = "Manuscript/Plots/")
