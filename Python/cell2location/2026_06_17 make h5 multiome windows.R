library(Seurat)
library(zellkonverter)   # BiocManager::install("zellkonverter")
library(SingleCellExperiment)

# ── Multiome → h5ad ────────────────────────────────────────────────────────
multiome <- readRDS("C:/Users/Gabe/Desktop/nemo.orig_harmony.integration_all_testd_clusters.rds")
DefaultAssay(multiome) <- "RNA"

# Seurat v5: join layers first (RNA may be split into counts.1, counts.2, etc.)
multiome[["RNA"]] <- JoinLayers(multiome[["RNA"]])

# Convert to SCE, then write h5ad
sce_ref <- as.SingleCellExperiment(multiome, assay = "RNA")
writeH5AD(sce_ref, file = "C:/Users/Gabe/Desktop/multiome_ref.h5ad")


