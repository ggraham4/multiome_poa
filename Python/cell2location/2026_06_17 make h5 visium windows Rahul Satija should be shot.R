library(Matrix)
library(tidyverse)

# ── Load raw counts from Seurat object without using Seurat methods ─────────
vis <- readRDS("C:/Users/Gabe/Desktop/Visium/allcells_norm_noInteg_2026-05-14.rds")

# Extract counts matrix directly from the Spatial.Polygons assay
# Seurat v5 stores layers, not slots
counts_mat <- vis[["Spatial.Polygons"]]$counts

cat("Counts matrix dims:", dim(counts_mat), "\n")  # genes x cells

# Extract metadata
meta <- vis@meta.data
cat("Metadata dims:", dim(meta), "\n")
cat("Metadata columns:", colnames(meta), "\n")

# ── Load dissection CSVs ────────────────────────────────────────────────────
sections <- c('6P17', '3P5', '4P5', '4P10')
tog <- data.frame()
for (i in sections) {
  data <- read.csv(paste0('C:/Users/Gabe/Desktop/multiome_poa/Visium/Groupings/', i, ' Coltan Dissection.csv'))
  colnames(data) <- gsub("Coltan\\.D.*ssection", "Coltan.DIssection", colnames(data), ignore.case = TRUE)
  tog <- rbind(tog, data)
}

cat("Dissection CSV rows:", nrow(tog), "\n")

# ── Match barcodes ──────────────────────────────────────────────────────────
# meta rownames have suffix _1, tog$Barcode does not
meta_barcodes <- rownames(meta)
stripped      <- gsub("_[0-9]+$", "", meta_barcodes)

meta$dissection <- unname(tog$Coltan.DIssection[match(stripped, tog$Barcode)])

cat("Matched dissection:\n")
print(table(meta$dissection, useNA = "always"))

# ── Subset to included cells ────────────────────────────────────────────────
keep      <- meta_barcodes[!is.na(meta$dissection)]
meta_sub  <- meta[keep, ]
counts_sub <- counts_mat[, keep]

cat("Subsetted counts:", dim(counts_sub), "\n")
cat("Subsetted metadata:", dim(meta_sub), "\n")

# ── Pure R fallback: write pieces to disk ───────────────────────────────────
library(Matrix)

writeMM(counts_sub, "C:/Users/Gabe/Desktop/Visium/counts.mtx")
write.csv(meta_sub, "C:/Users/Gabe/Desktop/Visium/meta.csv")
write.csv(data.frame(gene = rownames(counts_sub)), "C:/Users/Gabe/Desktop/Visium/genes.csv")
write.csv(data.frame(cell = colnames(counts_sub)), "C:/Users/Gabe/Desktop/Visium/cells.csv")
cat("Done. Pieces written.\n")


writeMM(counts_sub, "C:/Users/Gabe/Desktop/Visium/sketch_counts.mtx")
write.csv(meta_sub, "C:/Users/Gabe/Desktop/Visium/sketch_meta.csv")
write.csv(data.frame(gene = rownames(counts_sub)), "C:/Users/Gabe/Desktop/Visium/sketch_genes.csv")
write.csv(data.frame(cell = colnames(counts_sub)), "C:/Users/Gabe/Desktop/Visium/sketch_cells.csv")
cat("Done. Pieces written.\n")