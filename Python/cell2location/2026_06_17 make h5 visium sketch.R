# 01_export_sketch_for_cell2location.R

library(Matrix)
library(tidyverse)

cat("Loading Seurat object...\n")

vis <- readRDS(
  "C:/Users/Gabe/Desktop/Visium/allcells_norm_noInteg_2026-05-14.rds"
)

# --------------------------------------------------------------------------

# Use SKETCH assay instead of full Spatial.Polygons assay

# --------------------------------------------------------------------------

counts_mat <- vis[["sketch"]]$counts

cat("\nSketch assay loaded\n")
cat("Genes:", nrow(counts_mat), "\n")
cat("Cells:", ncol(counts_mat), "\n")

# --------------------------------------------------------------------------

# Metadata

# --------------------------------------------------------------------------

meta <- vis@meta.data
cat("\nMetadata loaded\n")
cat("Rows:", nrow(meta), "\n")
cat("Columns:", ncol(meta), "\n")

# --------------------------------------------------------------------------

# Verify sketch cells exist in metadata

# --------------------------------------------------------------------------

shared_cells <- intersect(
  colnames(counts_mat),
  rownames(meta)
)

cat(
  "\nSketch cells found in metadata:",
  length(shared_cells),
  "of",
  ncol(counts_mat),
  "\n"
)

if (length(shared_cells) != ncol(counts_mat)) {
  stop("Not all sketch cells were found in metadata.")
}

# --------------------------------------------------------------------------

# Load dissection annotations

# --------------------------------------------------------------------------

cat("\nLoading dissection annotations...\n")

sections <- c(
  "6P17",
  "3P5",
  "4P5",
  "4P10"
)

tog <- data.frame()

for (section in sections) {
  
  f <- paste0(
    "C:/Users/Gabe/Desktop/multiome_poa/Visium/Groupings/",
    section,
    " Coltan Dissection.csv"
  )
  
  d <- read.csv(f)
  
  colnames(d) <- gsub(
    "Coltan\\.D.*ssection",
    "Coltan.DIssection",
    colnames(d),
    ignore.case = TRUE
  )
  
  tog <- rbind(tog, d)
}

cat("Dissection rows loaded:", nrow(tog), "\n")

# --------------------------------------------------------------------------

# Match dissection labels

# --------------------------------------------------------------------------

meta_barcodes <- rownames(meta)

stripped <- gsub(
  "_[0-9]+$",
  "",
  meta_barcodes
)

meta$dissection <- unname(
  tog$Coltan.DIssection[
    match(stripped, tog$Barcode)
  ]
)

cat("\nDissection summary:\n")
print(table(meta$dissection, useNA = "always"))

# --------------------------------------------------------------------------

# Restrict to sketch cells

# --------------------------------------------------------------------------

meta <- meta[shared_cells, , drop = FALSE]

# --------------------------------------------------------------------------

# Keep only cells with dissection labels

# --------------------------------------------------------------------------

keep <- rownames(meta)[
  !is.na(meta$dissection)
]

meta_sub <- meta[keep, , drop = FALSE]

counts_sub <- counts_mat[, keep]

cat("\nAfter filtering:\n")
cat("Genes:", nrow(counts_sub), "\n")
cat("Cells:", ncol(counts_sub), "\n")

# --------------------------------------------------------------------------

# Output directory

# --------------------------------------------------------------------------

outdir <- "C:/Users/Gabe/Desktop/Visium"

dir.create(
  outdir,
  recursive = TRUE,
  showWarnings = FALSE
)

# --------------------------------------------------------------------------

# Export files for AnnData reconstruction

# --------------------------------------------------------------------------

cat("\nWriting matrix...\n")

writeMM(
  counts_sub,
  file.path(outdir, "sketch_counts.mtx")
)

cat("Writing metadata...\n")

write.csv(
  meta_sub,
  file.path(outdir, "sketch_meta.csv")
)

write.csv(
  data.frame(
    gene = rownames(counts_sub)
  ),
  file.path(outdir, "sketch_genes.csv"),
  row.names = FALSE
)

write.csv(
  data.frame(
    cell = colnames(counts_sub)
  ),
  file.path(outdir, "sketch_cells.csv"),
  row.names = FALSE
)

# --------------------------------------------------------------------------

# Summary

# --------------------------------------------------------------------------

cat("\nFinished.\n")

cat(
  "\nMatrix:",
  file.path(outdir, "sketch_counts.mtx"),
  "\n"
)

cat(
  "Metadata:",
  file.path(outdir, "sketch_meta.csv"),
  "\n"
)

cat(
  "Genes:",
  file.path(outdir, "sketch_genes.csv"),
  "\n"
)

cat(
  "Cells:",
  file.path(outdir, "sketch_cells.csv"),
  "\n"
)

cat(
  "\nFinal matrix dimensions:",
  nrow(counts_sub),
  "genes x",
  ncol(counts_sub),
  "cells\n"
)
