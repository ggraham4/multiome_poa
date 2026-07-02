library(Seurat)
library(zellkonverter)   # BiocManager::install("zellkonverter")
library(SingleCellExperiment)

# ── Multiome → h5ad ────────────────────────────────────────────────────────
multiome <- readRDS("C:/Users/Gabe/Desktop/nemo.orig_harmony.integration_all_testd_clusters.rds")
DefaultAssay(multiome) <- "RNA"

### subset to rgc
multiome = FindSubCluster(multiome,
                     1, 'harmony.wsnn', resolution = 0.2, subcluster.name = 'sub_res0.2')

sub_1 = subset(multiome, res0.8_50nn_40PC_45LSI%in%c(1,15,26))

Idents(sub_1) = 'sub_res0.2'
DimPlot(sub_1, label = T, reduction = 'harmony_wnn.umap')

# label 
labels  = c(
  '1_1' = '1_NSC',
  '1_0' = '1_NP',
  '1_2' = '1_shha',
  '1_3' = '1_gfap',
  '15' ='EC',
  '26' = 'DG'
)

sub_1@meta.data$rgc_label = unname(labels[sub_1@meta.data$sub_res0.2])


# Seurat v5: join layers first (RNA may be split into counts.1, counts.2, etc.)
sub_1[["RNA"]] <- JoinLayers(sub_1[["RNA"]])

# Convert to SCE, then write h5ad
sce_ref <- as.SingleCellExperiment(sub_1, assay = "RNA")
writeH5AD(sce_ref, file = "C:/Users/Gabe/Desktop/RGC_ref.h5ad")


