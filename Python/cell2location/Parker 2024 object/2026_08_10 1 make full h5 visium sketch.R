# 01_export_sketch_for_cell2location.R

library(Seurat)
library(Matrix)
library(tidyverse)

cat("Loading Seurat object...\n")

vis <- readRDS(
"C:/Users/Gabe/Desktop/Visium/vis_better_barcodes.rds"
)

SpatialDimPlot(vis, images ='s_6P17.polygons')
# perfect its all the cells

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

outdir <- "C:/Users/Gabe/Desktop/Visium/Parker Object"

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
