library(Seurat)
library(tidyverse)
library(Matrix)

# ── Paths ──────────────────────────────────────────────────────────────────
visium_dir    <- "C:/Users/Gabe/Desktop/Visium/"
results_dir   <- "C:/Users/Gabe/Desktop/Visium/parker_clusters_49/cell2location_map"
groupings_dir <- "C:/Users/Gabe/Desktop/multiome_poa/Visium/Groupings"
plot_dir      <- "C:/Users/Gabe/Desktop/Visium/plots"
dir.create(plot_dir, showWarnings = FALSE)

# ── Load data ──────────────────────────────────────────────────────────────
cat("Loading Seurat object...\n")
vis_full <- readRDS(file.path(visium_dir, "vis_better_barcodes.rds"))

#vis <- subset(
#  vis_full,
#  cells = Cells(vis)[vis$dissection == "In"] # comment this out because we are no longer doing dissection only
#)
vis = vis_full

cat("Loading c2l high-confidence (q05) abundance...\n")
abundance <- read.csv(
  file.path(results_dir, "abundance_q05.csv"),
  row.names = 1,
  check.names = FALSE
)
cat("Abundance dims:", dim(abundance), "\n")
cat("Abundance barcode example:", head(rownames(abundance), 3), "\n")

# ── Determine Dominant Cell Types ──────────────────────────────────────────
cat("\nIdentifying dominant cell type per sketch spot...\n")

# Find the column index with the highest q05 value for each cell/row
max_idx <- max.col(abundance, ties.method = "first")
dominant_labels <- colnames(abundance)[max_idx]
names(dominant_labels) <- rownames(abundance)

# Add this categorical label to the Seurat object's metadata
vis$dominant_c2l_cluster <- NA
vis$dominant_c2l_cluster[names(dominant_labels)] <- dominant_labels

cat("Unique cell types assigned to sketch:", length(unique(na.omit(vis$dominant_c2l_cluster))), "\n")

# ── Project Labels from Sketch to Full Object ──────────────────────────────
cat("\nProjecting dominant classifications to full object...\n")

vis <- ProjectData(
  object             = vis,
  assay              = "Spatial.Polygons",
  full.reduction     = "full.pca.sketch",
  sketched.assay     = "sketch",          # Targets the original sketch assay
  sketched.reduction = "pca.sketch",
  dims               = 1:15,
  # Pass the metadata column string so Seurat maps the categorical predictions
  refdata            = list(c2l_predicted_id = "dominant_c2l_cluster")
)

cat("Projected metadata column 'c2l_predicted_id' added successfully.\n")

# Verify predictions populated
cat("Assigned spots in full object:", sum(!is.na(vis$c2l_predicted_id)), "\n")

# ── Replot with Projected Categories ───────────────────────────────────────
# Make sure available_fovs is defined or grab them from the object
if (!exists("available_fovs")) {
  available_fovs <- Images(vis)
}

for (fov in available_fovs) {
  fov_clean <- gsub("s_|\\.polygons", "", fov)
  fov_dir   <- file.path(plot_dir, paste0(fov_clean, "_projected_clusters"))
  dir.create(fov_dir, showWarnings = FALSE)
  DefaultFOV(vis) <- fov
  
  cat("Plotting projected clusters for FOV:", fov, "\n")
  
  # Generate a discrete spatial plot of our projected classifications
  tryCatch({
    p <- ImageDimPlot(
      object = vis, 
      group.by = "c2l_predicted_id", 
      size = 1.5,
    ) + ggtitle(paste0(fov_clean, " - Projected C2L Clusters"))
    
    outfile <- file.path(fov_dir, "projected_c2l_map.png")
    ggsave(outfile, p, width = 10, height = 8, dpi = 300)
    cat("  Saved cluster map to:", outfile, "\n")
  }, error = function(e) cat("  Failed plotting FOV map:", fov, "-", conditionMessage(e), "\n"))
}

cat("\nPipeline complete!\n")

saveRDS(vis, file.path(visium_dir, "vis_better_barcodes_c2l05_projection_parker_2024_clusters49.rds"))


