library(Seurat)
library(tidyverse)
library(Matrix)

# ── Paths ──────────────────────────────────────────────────────────────────
visium_dir    <- "~/Desktop/Visium/"
results_dir   <- "/Desktop/Desktop/multiome_poa/Python/cell2location/results_2026_06_17/cell2location_map"
groupings_dir <- "/Desktop/Desktop/multiome_poa/Visium/Groupings"
plot_dir      <- file.path(visium_dir, "plots")
dir.create(plot_dir, showWarnings = FALSE)

# ── Load data ──────────────────────────────────────────────────────────────
cat("Loading Seurat object...\n")
vis <- readRDS(file.path(visium_dir, "vis_dissection.rds"))

cat("Loading c2l abundance...\n")
abundance <- read.csv(
  "Python/cell2location/results_2026_06_17/cell2location_map/abundance_means.csv",
  row.names = 1,
  check.names = FALSE
)
cat("Abundance dims:", dim(abundance), "\n")
cat("Abundance barcode example:", head(rownames(abundance), 3), "\n")

# ── Build full-sized abundance matrix ──────────────────────────────────────
# rows = cell types, cols = all cells in full object
cat("\nBuilding full abundance matrix...\n")

full_barcodes <- rownames(vis@meta.data)
cell_types    <- colnames(abundance)
abund_t       <- t(as.matrix(abundance))   # cell types x sketch cells

glimpse(abund_t)
# Initialise with zeros for all full-object cells
full_abund <- matrix(
  0,
  nrow     = length(cell_types),
  ncol     = length(full_barcodes),
  dimnames = list(cell_types, full_barcodes)
)

# Fill in sketch values where barcodes match
shared <- intersect(colnames(abund_t), full_barcodes)
cat("Sketch cells matching full object:", length(shared), "\n")
full_abund[, shared] <- abund_t[, shared]

cat("Non-zero cells in first cell type:", sum(full_abund[1, ] > 0), "\n")

# ── Add as Seurat assay ────────────────────────────────────────────────────
cat("\nAdding c2l assay to Seurat object...\n")
vis[["c2l"]] <- CreateAssayObject(counts = full_abund)
DefaultAssay(vis) <- "c2l"

# Feature names get underscores replaced with dashes by Seurat — use actual rownames
ct_features <- rownames(vis[["c2l"]])
cat("Cell type features:", head(ct_features), "\n")

# ── Plot each cell type per FOV ────────────────────────────────────────────
available_fovs <- names(vis@images)
cat("\nAvailable FOVs:", available_fovs, "\n")

for (fov in available_fovs) {
  fov_clean <- gsub("s_|\\.polygons", "", fov)
  fov_dir   <- file.path(plot_dir, fov_clean)
  dir.create(fov_dir, showWarnings = FALSE)

  DefaultFOV(vis) <- fov
  cat("Plotting FOV:", fov, "\n")

  for (ct in ct_features) {
    ct_clean <- gsub("meanscell-abundance-w-sf-", "cluster_", ct)
    outfile  <- file.path(fov_dir, paste0(ct_clean, ".png"))

    tryCatch({
      p <- ImageFeaturePlot(
        vis,
        features   = ct,
        max.cutoff = "q95"
      ) + ggtitle(paste0(fov_clean, " - ", ct_clean))

      ggsave(outfile, p, width = 8, height = 8, dpi = 300)

    }, error = function(e) {
      cat("  Failed:", ct, "-", conditionMessage(e), "\n")
    })
  }
  cat("  Done.", length(ct_features), "plots saved to", fov_dir, "\n")
}

SpatialDimPlot(vis, 'dissection', images = available_fovs[1])

# ── Summary heatmap by dissection region ──────────────────────────────────
cat("\nGenerating dissection region summary...\n")

region_abund <- read.csv(
  file.path(results_dir, "abundance_by_dissection_region.csv"),
  row.names = 1
)

if (nrow(region_abund) > 0) {
  library(pheatmap)

  mat <- t(scale(t(as.matrix(region_abund))))
  mat[is.nan(mat)] <- 0

  pheatmap(
    mat,
    cluster_rows = F,
    cluster_cols = F,
    main         = "Cell type abundance by dissection region (z-scored)",
    filename     = file.path(plot_dir, "region_abundance_heatmap.png"),
    width        = 12,
    height       = 6
  )
  cat("Heatmap saved.\n")
}

cat("\nAll done. Plots saved to:", plot_dir, "\n")