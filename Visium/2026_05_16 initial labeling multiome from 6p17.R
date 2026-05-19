# =============================================================================
# PSEUDOBULK / AGGREGATE EXPRESSION PIPELINE
# Multiome obj (final_clusters) vs Visium 6P17 anatomical regions
# =============================================================================

library(Seurat)
library(dplyr)
library(pheatmap)
library(ggplot2)
library(patchwork)

options(future.globals.maxSize = Inf)

# =============================================================================
# 1. LOAD OBJECTS
# =============================================================================

obj <- readRDS("~/Desktop/optimal_clustering_rna_only.rds")
Idents(obj) <- "final_clusters"

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

seurat_barcodes <- gsub("_4$", "", Cells(vis_6p17))

vis_6p17$anatomical <- unname(ana[seurat_barcodes])
vis_6p17$anatomical[is.na(vis_6p17$anatomical)] <- "Unannotated"

# =============================================================================
# 3. SCT ON sketch ASSAY (Visium)
# =============================================================================

DefaultAssay(vis_6p17) <- "sketch"

vis_6p17 <- SCTransform(
  vis_6p17,
  assay   = "sketch",
  verbose = FALSE
)

# =============================================================================
# 4. SCT ON MULTIOME REFERENCE
# =============================================================================

DefaultAssay(obj) <- "RNA"

obj <- SCTransform(
  obj,
  verbose = FALSE
)

# =============================================================================
# 5. FILTERING
# =============================================================================

vis_6p17 <- subset(
  vis_6p17,
  subset = anatomical != "Unannotated"
)

cluster_sizes <- table(obj$final_clusters)
keep_clusters <- names(cluster_sizes[cluster_sizes >= 100])
obj           <- subset(obj, subset = final_clusters %in% keep_clusters)

# =============================================================================
# 6. PSEUDOBULK VISIUM REGIONS (sketch cells only)
# =============================================================================

sketch_cells <- colnames(vis_6p17[["sketch"]])
vis_sketch   <- subset(vis_6p17, cells = sketch_cells)

region_avg <- AggregateExpression(
  vis_sketch,
  group.by      = "anatomical",
  assays        = "SCT",
  slot          = "data",
  return.seurat = FALSE
)$SCT

# =============================================================================
# 7. PSEUDOBULK MULTIOME CLUSTERS
# =============================================================================

cluster_avg <- AggregateExpression(
  obj,
  group.by      = "final_clusters",
  assays        = "SCT",
  slot          = "data",
  return.seurat = FALSE
)$SCT

# =============================================================================
# 8. MATCH GENES
# =============================================================================

common_genes <- intersect(rownames(region_avg), rownames(cluster_avg))
region_avg   <- region_avg[common_genes, ]
cluster_avg  <- cluster_avg[common_genes, ]

# =============================================================================
# 9. VARIABLE GENES ONLY
# =============================================================================

var_genes  <- VariableFeatures(obj)
common_var <- intersect(common_genes, var_genes)
region_avg <- region_avg[common_var, ]
cluster_avg <- cluster_avg[common_var, ]

cat("Common variable genes used:", length(common_var), "\n")

# =============================================================================
# 10. SPEARMAN CORRELATION
# =============================================================================

cor_mat <- cor(
  as.matrix(region_avg),
  as.matrix(cluster_avg),
  method = "spearman"
)

# =============================================================================
# 11. HEATMAP: Visium regions → Multiome clusters
# =============================================================================

best_cluster <- apply(cor_mat, 1, function(x) names(which.max(x)))

row_order        <- order(best_cluster)
ordered_clusters <- unique(best_cluster[row_order])

diag_mat <- cor_mat[rownames(cor_mat)[row_order], ordered_clusters, drop = FALSE]

diag_scaled <- t(scale(t(diag_mat)))

pheatmap(
  diag_scaled,
  main          = "Visium regions → Multiome final_clusters"
)

# =============================================================================
# 12. HEATMAP: Multiome clusters → Visium regions
# =============================================================================

cor_rev    <- t(cor_mat)
best_region <- apply(cor_rev, 1, function(x) names(which.max(x)))

row_order2      <- order(best_region)
ordered_regions <- unique(best_region[row_order2])

diag_mat2 <- cor_rev[rownames(cor_rev)[row_order2], ordered_regions, drop = FALSE]

diag_scaled2 <- t(scale(t(diag_mat2)))

pheatmap(
  diag_scaled2,
  cluster_rows  = FALSE,
  cluster_cols  = FALSE,
  border_color  = NA,
  fontsize_row  = 10,
  fontsize_col  = 10,
  main          = "Multiome final_clusters → Visium regions"
)

# =============================================================================
# 13. MUTUAL BEST MATCHES
# =============================================================================

best_cluster <- apply(cor_mat, 1, function(x) names(which.max(x)))
best_region  <- apply(cor_mat, 2, function(x) names(which.max(x)))

mutual_df <- data.frame(
  anatomical   = names(best_cluster),
  best_cluster = best_cluster
)

mutual_df$reverse_match <- best_region[mutual_df$best_cluster]
mutual_df$mutual_best   <- mutual_df$anatomical == mutual_df$reverse_match

print(mutual_df)

# =============================================================================
# 14. MUTUAL BEST HEATMAP
# =============================================================================

keep_regions  <- mutual_df$anatomical[mutual_df$mutual_best]
keep_clusters <- mutual_df$best_cluster[mutual_df$mutual_best]

mutual_cor <- cor_mat[keep_regions, keep_clusters, drop = FALSE]

pheatmap(
  mutual_cor,
  cluster_rows  = FALSE,
  cluster_cols  = FALSE,
  border_color  = NA,
  fontsize_row  = 10,
  fontsize_col  = 10,
  main          = "Mutual best: Visium anatomical ↔ Multiome final_clusters"
)


pheatmap(
  cor_mat%>%scale(),
  cluster_rows  = FALSE,
  cluster_cols  = FALSE,
  border_color  = NA,
  fontsize_row  = 10,
  fontsize_col  = 10,
  main          = "Cor Mat"
)

# =============================================================================
# 15. LABEL MULTIOME CLUSTERS BY BEST VISIUM REGION (from diag_scaled)
# =============================================================================

# Use diag_scaled — matches what eye reads off heatmap
# Each multiome cluster gets the region it correlates most with
# Many clusters can share the same region (no 1:1 constraint)

best_region_per_cluster <- apply(diag_scaled, 2, function(x) names(which.max(x)))
names(best_region_per_cluster) <- gsub("^g", "", names(best_region_per_cluster))

obj$visium_region <- unname(
  best_region_per_cluster[as.character(obj$final_clusters)]
)
obj$visium_region[is.na(obj$visium_region)] <- "Unassigned"

cat("\n--- Multiome cluster → Visium region assignments ---\n")
print(table(obj$visium_region))

# =============================================================================
# 16. DIMPLOT
# =============================================================================

x11(width = 16, height = 7)

DimPlot(obj,
        group.by   = "visium_region",
        label      = TRUE,
        repel      = TRUE,
        label.size = 3,
        reduction  = "harmony_wnn.umap") +
  labs(title = "Multiome clusters labeled by best Visium anatomical region") +
  theme(legend.position = "right") +

DimPlot(obj,
        group.by   = "final_clusters",
        label      = TRUE,
        repel      = TRUE,
        label.size = 3,
        reduction  = "harmony_wnn.umap") +
  labs(title = "Multiome final_clusters (numbered)") +
  theme(legend.position = "none")


# ── Full heatmap — all clusters retained ─────────────────────────────────────

# Order rows by best cluster (same as before)
best_cluster <- apply(cor_mat, 1, function(x) names(which.max(x)))
row_order    <- order(best_cluster)

# Order columns numerically — keep ALL clusters
col_order <- order(as.numeric(gsub("^g", "", colnames(cor_mat))))

diag_mat_full <- cor_mat[rownames(cor_mat)[row_order], col_order, drop = FALSE]

# Row scale
diag_scaled_full <- scale(t(scale(t(diag_mat_full))))

pheatmap(
  diag_scaled_full,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  border_color = NA,
  fontsize_row = 10,
  fontsize_col = 10,
  main         = "Visium regions → Multiome final_clusters (all clusters)"
)

# ── Label assignment now uses full matrix ─────────────────────────────────────
best_region_per_cluster <- apply(diag_scaled_full, 2, function(x) names(which.max(x)))
names(best_region_per_cluster) <- gsub("^g", "", names(best_region_per_cluster))

obj$visium_region <- unname(
  best_region_per_cluster[as.character(obj$final_clusters)]
)
obj$visium_region[is.na(obj$visium_region)] <- "Unassigned"

table(obj$visium_region)


vis_6p17$anatomical = ifelse(vis_6p17@assays$Spatial.Polygons$data['s100b',]>0, 'Ependymoglia', vis_6p17$anatomical )


vis_6p17_radial= subset(vis_6p17, anatomical %in% c("Ependymoglia"))

DimPlot(vis_6p17, group.by = 'anatomical')

DimPlot(vis_6p17,
        group.by = 'anatomical',
        cells.highlight = rownames(vis_6p17@meta.data[vis_6p17$anatomical=='Ependymoglia',]),
        pt.size = 0.01,
          sizes.highlight = 0.01
        )
vis_6p17_radial$anatomical = ifelse(vis_6p17_radial@assays$Spatial.Polygons$data['s100b',]>0, 'Ependymoglia', vis_6p17_radial$anatomical )
vis_6p17_radial$anatomical = ifelse(vis_6p17_radial@assays$Spatial.Polygons$data['LOC111577263',]>0, 'Arob+', vis_6p17_radial$anatomical )
vis_6p17_radial$anatomical = ifelse(vis_6p17_radial@assays$Spatial.Polygons$data['gfap',]>0, 'gfap+', vis_6p17_radial$anatomical )
vis_6p17_radial$anatomical = ifelse(vis_6p17_radial@assays$Spatial.Polygons$data['gja1b',]>0, 'gja1b+', vis_6p17_radial$anatomical )
vis_6p17_radial$anatomical = ifelse(vis_6p17_radial@assays$Spatial.Polygons$data['sox2',]>0, 'sox2+', vis_6p17_radial$anatomical )
vis_6p17_radial$anatomical = ifelse(vis_6p17_radial@assays$Spatial.Polygons$data['crocc2',]>0, 'EC', vis_6p17_radial$anatomical )
vis_6p17_radial$anatomical = ifelse(vis_6p17_radial@assays$Spatial.Polygons$data['sv2a',]>0, 'Neuron', vis_6p17_radial$anatomical )
vis_6p17_radial$anatomical = ifelse(vis_6p17_radial@assays$Spatial.Polygons$data['nkx2.1',]>0, 'nkx2.1', vis_6p17_radial$anatomical )
vis_6p17_radial$anatomical = ifelse(vis_6p17_radial@assays$Spatial.Polygons$data['nkx2.1',]>0, 'nkx2.1', vis_6p17_radial$anatomical )

DimPlot(vis_6p17_radial, group.by = 'anatomical')
DimPlot(vis_6p17_radial)


