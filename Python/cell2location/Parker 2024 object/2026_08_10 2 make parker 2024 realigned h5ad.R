library(Seurat)
library(zellkonverter)   
library(SingleCellExperiment)

# ── parker → h5ad ────────────────────────────────────────────────────────
parker <- readRDS("A:/Anemonefish POA Legacy R Objects/parker_2024_realigned_ocellaris.rds")
DefaultAssay(parker) <- "RNA"

sce_ref <- as.SingleCellExperiment(parker, assay = "RNA")
writeH5AD(sce_ref, file = "C:/Users/Gabe/Desktop/parker_realigned_ref.h5ad")


