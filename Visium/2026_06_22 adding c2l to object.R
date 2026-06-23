library(Seurat)
library(tidyverse)
library(Matrix)

# ── Paths ──────────────────────────────────────────────────────────────────
visium_dir    <- "C:/Users/Gabe/Desktop/Visium/"
results_dir   <- "C:/Users/Gabe/Desktop/multiome_poa/Python/cell2location/results_2026_06_17/cell2location_map"
groupings_dir <- "C:/Users/Gabe/Desktop/multiome_poa/Visium/Groupings"
plot_dir      <- "C:/Users/Gabe/Desktop/Visium/plots"
dir.create(plot_dir, showWarnings = FALSE)

# ── Load data ──────────────────────────────────────────────────────────────
cat("Loading Seurat object...\n")
vis <- readRDS(file.path(visium_dir, "vis_better_barcodes.rds"))

cat("Loading c2l high-confidence (q05) abundance...\n")
abundance <- read.csv(
  file.path(results_dir, "abundance_q05.csv"),
  row.names = 1,
  check.names = FALSE
)
cat("Abundance dims:", dim(abundance), "\n")

# ── Determine Dominant Cell Types ──────────────────────────────────────────
cat("\nIdentifying dominant cell type per sketch spot...\n")

max_idx <- max.col(abundance, ties.method = "first")
dominant_labels <- colnames(abundance)[max_idx]
names(dominant_labels) <- rownames(abundance)

# Add this categorical label to the Seurat object's metadata
vis$dominant_c2l_cluster <- NA
vis$dominant_c2l_cluster[names(dominant_labels)] <- dominant_labels

cat("Unique cell types assigned to sketch:", length(unique(na.omit(vis$dominant_c2l_cluster))), "\n")

# ── Project Labels Manually via Sketch Weights ─────────────────────────────
cat("\nComputing sketch weight matrix manually...\n")

# 1. Build or locate the exact KNN relation matrix between full object and sketch
vis <- FindSketchWeights(
  object             = vis,
  full.reduction     = "full.pca.sketch",
  sketched.reduction = "pca.sketch",
  dims               = 1:15,
  reduction.model    = NULL  # Forces recalculating clean dimensions
)

# 2. Extract the weight matrix (Rows = Full Cells, Columns = Sketch Cells)
weights <- GetCellRelation(vis, reference = "sketch")

# 3. Pull our sketch labels and match their order to the weight matrix columns
sketch_cells <- colnames(weights)
sketch_labels <- vis$dominant_c2l_cluster
names(sketch_labels) <- rownames(vis@meta.data) # reference global names
sketch_labels <- sketch_labels[sketch_cells]

# 4. Perform categorical label transfer based on max weight
cat("Transferring categorical labels using calculated weights...\n")

# Convert categorical labels to a binary design matrix (Cells x CellTypes)
# This lets us cleanly project categorical data mathematically
label_matrix <- model.matrix(~ 0 + sketch_labels)
colnames(label_matrix) <- gsub("sketch_labels", "", colnames(label_matrix))
rownames(label_matrix) <- sketch_cells

# Multiply weights (Full x Sketch) by label matrix (Sketch x CellTypes)
# Resulting matrix is (Full x CellTypes) representing probability/weight scores
projected_scores <- weights %% label_matrix

# Assign each full cell to its highest scoring projected cell type
full_projected_labels <- colnames(projected_scores)[max.col(projected_scores, ties.method = "first")]
names(full_projected_labels) <- rownames(weights)

# 5. Inject back into full object's metadata
vis$c2l_predicted_id <- NA
vis$c2l_predicted_id[names(full_projected_labels)] <- full_projected_labels

cat("Projected metadata column 'c2l_predicted_id' added successfully.\n")
cat("Assigned spots in full object:", sum(!is.na(vis$c2l_predicted_id)), "\n")

# ── Replot with Projected Categories ───────────────────────────────────────
if (!exists("available_fovs")) {
  available_fovs <- Images(vis)
}

for (fov in available_fovs) {
  fov_clean <- gsub("s_|\\.polygons", "", fov)
  fov_dir   <- file.path(plot_dir, paste0(fov_clean, "_projected_clusters"))
  dir.create(fov_dir, showWarnings = FALSE)
  DefaultFOV(vis) <- fov
  
  cat("Plotting projected clusters for FOV:", fov, "\n")
  
  tryCatch({
    p <- ImageDimPlot(
      object = vis, 
      group.by = "c2l_predicted_id", 
      size = 1.5,
      cols = "alphabet"
    ) + ggtitle(paste0(fov_clean, " - Projected C2L Clusters"))
    
    outfile <- file.path(fov_dir, "projected_c2l_map.png")
    ggsave(outfile, p, width = 10, height = 8, dpi = 300)
    cat("  Saved cluster map to:", outfile, "\n")
  }, error = function(e) cat("  Failed plotting FOV map:", fov, "-", conditionMessage(e), "\n"))
}

cat("\nPipeline complete!\n")