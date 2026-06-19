library(Seurat)
library(tidyverse)
library(zellkonverter) # Required for writeH5AD
library(SingleCellExperiment)

# === FIX FROM PART 1 ===
vis2 = readRDS("C:/Users/Gabe/Desktop/Visium/allcells_norm_noInteg_2026-05-14.rds")
vis = UpdateSeuratObject(vis2) # Prevents the 'coords_x_orientation' FOV error

images = unique(vis$Slice)
barcodes_included = data.frame()

for(i in images){
  print(i)
  data = read.csv(paste0("Visium/Groupings/", i, " Coltan DIssection.csv")) %>%
    subset(!is.na(Coltan.DIssection))
  data$Barcode_updated = paste0(data$Barcode, i)
  
  # FIX: Correctly append data so you don't overwrite previous iterations
  barcodes_included = rbind(barcodes_included, data) 
}

# ... [Your barcode renaming steps here] ...
barcodes = colnames(vis@assays$Spatial.Polygons)
barcodes_no_ <- sub("_[0-9]+$", "", barcodes)
barcodes_image = paste0(barcodes_no_, vis$Slice)

vis <- RenameCells(vis, new.names = barcodes_image)

vis$dissection = ifelse(colnames(vis) %in% barcodes_included$Barcode_updated, 'In', 'Out')
table(vis$dissection)

# === SAVING AND SUBSETTING ===
saveRDS(vis, "C:/Users/Gabe/Desktop/Visium/vis_better_barcodes.rds")
write.csv(barcodes_included, "Visium/Groupings/Coltan DIssections together and barcodes.csv")
write.csv(data.frame(colnames(vis)), "C:/Users/Gabe/Desktop/Visium/updated_barcodes.csv")

# Set assay and join layers safely
DefaultAssay(vis) <- "sketch"
vis[["sketch"]] <- JoinLayers(vis[["sketch"]])

# Subset the Seurat object
vis_sub = subset(vis, subset = dissection == 'In')

# Convert to SingleCellExperiment safely
# Note: If spatial image metadata causes conversion errors, 
# you can set assay="sketch" explicitly to just port the expression data.

library(Matrix)

# 1. Pull the subsetted metadata (17,589 rows)
metadata_sub <- vis_sub@meta.data

# 2. Bypass Seurat wrappers and pull the RAW counts slot directly
# In Seurat v5, multi-layer matrices are tucked away inside the 'layers' slot.
# Let's find the exact name of the raw counts layer:
layer_names <- names(vis[["sketch"]]@layers)
print(layer_names) # Likely "counts" or "counts.1", etc.

# Grab the first counts layer directly out of the S4 slot
raw_counts_matrix <- vis[["sketch"]]@layers[[layer_names[1]]]

# 3. Pull the TRUE cell names assigned specifically to this layer matrix
# This sidesteps colnames(vis) completely
true_matrix_cell_names <- colnames(vis[["sketch"]])

# 4. Force synchronization safely
colnames(raw_counts_matrix) <- true_matrix_cell_names

# 5. Extract ONLY your subsetted cells
# Intersecting ensures we only grab cells that actually exist in this matrix layer
common_barcodes <- intersect(rownames(metadata_sub), true_matrix_cell_names)
counts_sub <- raw_counts_matrix[, common_barcodes]
metadata_final <- metadata_sub[common_barcodes, ]

# ------------------------------------------------------------------
# WRITE INPUT FILES FOR PYTHON
# ------------------------------------------------------------------
output_dir <- "C:/Users/Gabe/Desktop/Visium"

# A. Save Metadata
write.csv(metadata_final, file = file.path(output_dir, "sketch_meta.csv"), row.names = TRUE)

# B. Save Gene Names
# Rows of the assay are always consistent with the matrix rows
genes_df <- data.frame(gene = rownames(vis[["sketch"]]))
write.csv(genes_df, file = file.path(output_dir, "sketch_genes.csv"), row.names = FALSE)

# C. Save Sparse Matrix
writeMM(counts_sub, file = file.path(output_dir, "sketch_counts.mtx"))

message(paste("Successfully exported", length(common_barcodes), "cells to MTX/CSV!"))