# =============================================================================
# RAM-SAFE PSEUDOBULK / AGGREGATE EXPRESSION PIPELINE
# Uses ONLY sketch assay for Visium SCT normalization
# No saving
# No full Spatial.Polygons SCT
# =============================================================================

library(Seurat)
library(dplyr)
library(pheatmap)
library(ggplot2)

options(future.globals.maxSize = Inf)

# =============================================================================
# 1. LOAD OBJECTS
# =============================================================================

parker_2024 <- readRDS("~/Desktop/parker_object.rds")
Idents(parker_2024) <- "clusters49"

vis <- readRDS('/Users/ggraham/Desktop/Visium/vis.rds')
vis_6p17 <- subset(vis, Slice == '6P17')
rm(vis)

# =============================================================================
# 2. LOAD ANATOMICAL ANNOTATIONS
# =============================================================================

maruska_6p17 <- read.csv(
  "Visium/Groupings/Maruska Regions Data Informed/2026_05_16 6p17 Anatomical.csv"
)

ana <- setNames(
  maruska_6p17$Anatomical,
  maruska_6p17$Barcode
)

seurat_barcodes <- gsub(
  "_4$",
  "",
  Cells(vis_6p17)
)

vis_6p17$anatomical <- unname(
  ana[seurat_barcodes]
)

vis_6p17$anatomical[is.na(vis_6p17$anatomical)] <- "Unannotated"

# =============================================================================
# 3. SCT ONLY ON sketch ASSAY
# =============================================================================

DefaultAssay(vis_6p17) <- "sketch"

vis_6p17 <- SCTransform(
  vis_6p17,
  assay = "sketch",
  verbose = FALSE
)

# =============================================================================
# 4. SCT ON REFERENCE
# =============================================================================

DefaultAssay(parker_2024) <- "RNA"

parker_2024 <- SCTransform(
  parker_2024,
  verbose = FALSE
)

# =============================================================================
# 5. OPTIONAL FILTERING
# =============================================================================

vis_6p17 <- subset(
  vis_6p17,
  subset = anatomical != "Unannotated"
)

cluster_sizes <- table(parker_2024$clusters49)

keep_clusters <- names(
  cluster_sizes[cluster_sizes >= 100]
)

parker_2024 <- subset(
  parker_2024,
  subset = clusters49 %in% keep_clusters
)

# =============================================================================
# 6. PSEUDOBULK VISIUM REGIONS
# =============================================================================

# IMPORTANT:
# We aggregate ONLY sketch cells
# because that is the SCT-normalized assay

sketch.cells <- colnames(
  vis_6p17[["sketch"]]
)

vis_sketch <- subset(
  vis_6p17,
  cells = sketch.cells
)

region_avg <- AggregateExpression(
  vis_sketch,
  group.by = "anatomical",
  assays = "SCT",
  slot = "data",
  return.seurat = FALSE
)$SCT

# =============================================================================
# 7. PSEUDOBULK snRNA CLUSTERS
# =============================================================================

cluster_avg <- AggregateExpression(
  parker_2024,
  group.by = "clusters49",
  assays = "SCT",
  slot = "data",
  return.seurat = FALSE
)$SCT

# =============================================================================
# 8. MATCH GENES
# =============================================================================

common_genes <- intersect(
  rownames(region_avg),
  rownames(cluster_avg)
)

region_avg <- region_avg[
  common_genes,
]

cluster_avg <- cluster_avg[
  common_genes,
]

# =============================================================================
# 9. USE VARIABLE GENES ONLY
# =============================================================================

var_genes <- VariableFeatures(parker_2024)

common_var <- intersect(
  common_genes,
  var_genes
)

region_avg <- region_avg[
  common_var,
]

cluster_avg <- cluster_avg[
  common_var,
]

# =============================================================================
# 10. SPEARMAN CORRELATION
# =============================================================================

cor_mat <- cor(
  region_avg%>%as.matrix(),
  cluster_avg%>%as.matrix(),
  method = "spearman"
)

# =============================================================================
# 11. ROW SCALING
# =============================================================================

cor_scaled <- t(
  scale(
    t(cor_mat)
  )
)

# =============================================================================
# 12. HEATMAP
# =============================================================================


# For each anatomical region:
# find best matching snRNA cluster

best_cluster <- apply(
  cor_mat,
  1,
  function(x) names(which.max(x))
)

# order rows by their best cluster
row_order <- order(best_cluster)

# order columns by first appearance in ordered rows
ordered_clusters <- unique(best_cluster[row_order])

# subset matrix
diag_mat <- cor_mat[
  rownames(cor_mat)[row_order],
  ordered_clusters,
  drop = FALSE
]

# -----------------------------------------------------------------------------
# Row-scale for visualization
# -----------------------------------------------------------------------------

diag_scaled <- t(
  scale(
    t(diag_mat)
  )
)

# -----------------------------------------------------------------------------
# Plot
# -----------------------------------------------------------------------------

pheatmap(
  diag_scaled,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  border_color = NA,
  fontsize_row = 10,
  fontsize_col = 10,
  main = "Visium regions → snRNA clusters"
)

# =============================================================================
# 2. snRNA -> VISIUM
# =============================================================================

# transpose matrix
cor_rev <- t(cor_mat)

# for each snRNA cluster:
# find best matching anatomical region

best_region <- apply(
  cor_rev,
  1,
  function(x) names(which.max(x))
)

# order rows
row_order2 <- order(best_region)

# ordered regions
ordered_regions <- unique(best_region[row_order2])

# subset matrix
diag_mat2 <- cor_rev[
  rownames(cor_rev)[row_order2],
  ordered_regions,
  drop = FALSE
]

# -----------------------------------------------------------------------------
# Row scaling
# -----------------------------------------------------------------------------

diag_scaled2 <- t(
  scale(
    t(diag_mat2)
  )
)

# -----------------------------------------------------------------------------
# Plot
# -----------------------------------------------------------------------------

pheatmap(
  diag_scaled2,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  border_color = NA,
  fontsize_row = 10,
  fontsize_col = 10,
  main = "snRNA clusters → Visium regions"
)
# =============================================================================
# 13. MUTUAL BEST MATCHES
# =============================================================================

best_cluster <- apply(
  cor_mat,
  1,
  function(x) names(which.max(x))
)

best_region <- apply(
  cor_mat,
  2,
  function(x) names(which.max(x))
)

mutual_df <- data.frame(
  anatomical = names(best_cluster),
  best_cluster = best_cluster
)

mutual_df$reverse_match <- best_region[
  mutual_df$best_cluster
]

mutual_df$mutual_best <- (
  mutual_df$anatomical ==
    mutual_df$reverse_match
)

mutual_df

# =============================================================================
# 14. VISUALIZE ONLY MUTUAL BEST
# =============================================================================
library(patchwork)

keep_regions <- mutual_df$anatomical[
  mutual_df$mutual_best
]

keep_clusters <- mutual_df$best_cluster[
  mutual_df$mutual_best
]

mutual_cor <- cor_mat[
  keep_regions,
  keep_clusters,
  drop = FALSE
]

a= pheatmap(
  mutual_cor,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  border_color = NA,
  fontsize_row = 10,
  fontsize_col = 10,
  main = "Mutual best anatomical ↔ snRNA matches"
)
b = DimPlot(parker_2024, label = T)

x11(width=8)
b

# Best parker cluster for each anatomical region, from the scaled matrix
best_cluster_scaled <- apply(diag_scaled, 1, function(x) names(which.max(x)))

# Strip g prefix
best_cluster_scaled_clean <- gsub("^g", "", best_cluster_scaled)

# Reverse: for each parker cluster, which region best matches it?
# Multiple clusters can share the same region
best_region_per_cluster_scaled <- apply(diag_scaled, 2, function(x) names(which.max(x)))
names(best_region_per_cluster_scaled) <- gsub("^g", "", names(best_region_per_cluster_scaled))

# Assign to parker
parker_2024$visium_region <- best_region_per_cluster_scaled[as.character(parker_2024$clusters49)]%>%unname()
parker_2024$visium_region[is.na(parker_2024$visium_region)] <- "Unassigned"

# Verify POA clusters are now correct
table(parker_2024$visium_region)

# Plot
x11(width = 14)
DimPlot(parker_2024,
        group.by   = "visium_region",
        label      = TRUE,
        repel      = TRUE,
        label.size = 3) +
  labs(title = "Parker clusters labeled by best Visium region (row-scaled)") +
  theme(legend.position = "right")
