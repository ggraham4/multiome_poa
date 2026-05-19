library(Seurat)
library(dplyr)
library(pheatmap)
library(ggplot2)
library(patchwork)
`%notin%` = Negate(`%in%`)

# ── 1. Load ──────────────────────────────────────────────────────────────────
vis <- readRDS('/Users/ggraham/Desktop/Visium/vis.rds')

# ── 2. Anatomical labels ──────────────────────────────────────────────────────
barcodes <- list.files(
  '/Users/ggraham/Desktop/multiome_poa/Visium/Groupings/Maruska Regions Data Informed/',
  full.names = TRUE
) %>%
  lapply(read.csv) %>%
  bind_rows()

ana <- setNames(barcodes$Anatomical, barcodes$Barcode)
seurat_barcodes <- gsub("_[0-9]+$", "", Cells(vis))
vis$anatomical <- unname(ana[seurat_barcodes])

# ── 3. QC: how many spots got labels? ────────────────────────────────────────
cat("Labeled spots:", sum(!is.na(vis$anatomical)), "\n")
cat("Unlabeled spots:", sum(is.na(vis$anatomical)), "\n")

# ── 4. Non-neuronal clusters (from marker genes) ─────────────────────────────
Idents(vis) <- 'seurat_clusters.projected0.4'
DotPlot(vis, features = c('elavl3','rbfox3a','gfap','LOC111577263',
                           'mpz','ptprc','p2ry12','cspg4','hbae5')) +
  coord_flip()

non_neuronal <- c(1, 4, 30, 15, 14)

# ── 5. Subset: neurons with anatomical labels only ───────────────────────────
vis_neuron <- subset(vis,
  (seurat_clusters.projected0.4 %notin% non_neuronal) &
  (!is.na(anatomical))
)
cat("Neurons with labels:", ncol(vis_neuron), "\n")
cat("Regions present:", table(vis_neuron$anatomical), "\n")

# ── 6. Recluster on sketch assay ─────────────────────────────────────────────
DefaultAssay(vis_neuron) <- "sketch"
vis_neuron <- FindVariableFeatures(vis_neuron)
vis_neuron <- ScaleData(vis_neuron, verbose = FALSE)
vis_neuron <- RunPCA(vis_neuron, reduction.name = "pca.neuron", verbose = FALSE)
ElbowPlot(vis_neuron, reduction = "pca.neuron")
# set dims after inspecting elbow plot — placeholder 15
ndims <- 15

vis_neuron <- FindNeighbors(vis_neuron, reduction = "pca.neuron", 
                             dims = 1:ndims, graph.name = "neuron_snn")
vis_neuron <- FindClusters(vis_neuron, graph.name = "neuron_snn",
                            resolution = 0.4, cluster.name = "neuron_clusters",
                            verbose = FALSE)
vis_neuron <- RunUMAP(vis_neuron, reduction = "pca.neuron", 
                       dims = 1:ndims, reduction.name = "umap.neuron")

# ── 7. Project to full data ───────────────────────────────────────────────────
vis_neuron <- ProjectData(
  object        = vis_neuron,
  assay         = "Spatial.Polygons",
  full.reduction   = "full.pca.neuron",
  sketched.assay   = "sketch",
  sketched.reduction = "pca.neuron",
  umap.model    = "umap.neuron",
  dims          = 1:ndims,
  refdata       = list(neuron_clusters.projected = "neuron_clusters")
)
Idents(vis_neuron) <- "neuron_clusters.projected"

# ── 8. Visualize ─────────────────────────────────────────────────────────────
DimPlot(vis_neuron, reduction = "umap.neuron", label = TRUE) + NoLegend()
DimPlot(vis_neuron, reduction = "umap.neuron", group.by = "anatomical", label = TRUE)

# ── 9. Anatomy-first heatmap ──────────────────────────────────────────────────
# rows = anatomy, cols = clusters, values = proportion of cluster in that region
regions <- table(vis_neuron$neuron_clusters.projected, vis_neuron$anatomical) %>%
  as.data.frame.matrix()

# normalize: what fraction of each REGION falls in each cluster
region_prop <- t(regions) / rowSums(t(regions))

# diagonal ordering
top_cluster  <- apply(region_prop, 1, which.max)   # best cluster per region
top_anatomy  <- apply(region_prop, 2, which.max)   # best region per cluster
col_order    <- order(top_anatomy)
cluster_pos  <- match(top_cluster, col_order)
row_order    <- order(cluster_pos)

pheatmap(region_prop[row_order, col_order],
         scale        = "none",
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         color        = colorRampPalette(c("grey95", "orange", "darkred"))(100),
         fontsize_row = 8,
         fontsize_col = 7,
         main         = "Fraction of region in each neuronal cluster")

# ── 10. Anatomy-driven markers ───────────────────────────────────────────────
# this is the key question: what genes distinguish regions?
Idents(vis_neuron) <- "anatomical"
region_markers <- FindAllMarkers(vis_neuron, 
                                  only.pos   = TRUE,
                                  min.pct    = 0.25,
                                  logfc.threshold = 0.25)
region_markers %>%
  group_by(cluster) %>%
  slice_max(avg_log2FC, n = 5) %>%
  print(n = Inf)



###
# try a range of resolutions
resolutions <- c(0.2, 0.4, 0.8, 1.2, 1.6, 2.0)

for(r in resolutions){
  vis_neuron <- FindClusters(vis_neuron,
                              graph.name   = "neuron_snn",
                              resolution   = r,
                              cluster.name = paste0("res_", r),
                              verbose      = FALSE)
  cat("Done:", r, "- clusters:", 
      nlevels(vis_neuron@meta.data[[paste0("res_", r)]]), "\n")
}

# quick visual sweep
res_cols <- paste0("res_", resolutions)
entropy <- function(p) {
  p <- p[p > 0]
  -sum(p * log(p))
}

res_scores <- sapply(res_cols, function(res) {
  vals <- vis_neuron@meta.data[[res]]
  if(sum(!is.na(vals)) == 0) return(NA)
  
  cells_use <- rownames(vis_neuron@meta.data)[!is.na(vals)]
  
  regions <- table(vis_neuron@meta.data[cells_use, res],
                   vis_neuron@meta.data[cells_use, "anatomical"]) %>%
    as.data.frame.matrix()
  
  # normalize each cluster to proportions
  region_prop <- regions / rowSums(regions)
  
  # mean entropy across clusters (weighted by cluster size)
  cluster_sizes <- rowSums(regions)
  cluster_entropies <- apply(region_prop, 1, entropy)
  weighted.mean(cluster_entropies, cluster_sizes)
})

data.frame(resolution = res_cols, weighted_entropy = res_scores) %>%
  ggplot(aes(x = as.numeric(gsub("res_", "", res_cols)), y = weighted_entropy)) +
  geom_line() + geom_point() +
  geom_vline(xintercept = as.numeric(gsub("res_", "", names(which.min(res_scores)))), 
             linetype = "dashed", color = "red") +
  labs(x = "Resolution", y = "Weighted mean entropy",
       title = "Lower = clusters map more cleanly to anatomical regions") +
  theme_classic()

# ok so res 2 seems best I guess

DimPlot(vis_neuron, group.by = 'res_2')
  
res = 'res_1.6'

vals <- vis_neuron@meta.data[[res]]

  cells_use <- rownames(vis_neuron@meta.data)[!is.na(vals)]
regions <- table(vis_neuron@meta.data[cells_use, res],
                   vis_neuron@meta.data[cells_use, "anatomical"]) %>%
    as.data.frame.matrix()
  
  region_prop <- t(regions) / rowSums(t(regions))
  region_prop[is.nan(region_prop)] <- 0
  
  top_cluster <- apply(region_prop, 1, which.max)
  top_anatomy <- apply(region_prop, 2, which.max)
  col_order   <- order(top_anatomy)
  row_order   <- order(match(top_cluster, col_order))
  
  pheatmap(region_prop[row_order, col_order],
           cluster_rows = FALSE,
           cluster_cols = FALSE,
           color  = colorRampPalette(c("grey95", "orange", "darkred"))(100),
           main   = paste("Resolution:", gsub("res_", "", res)),
           show_colnames = TRUE,
           fontsize_row = 5,
           fontsize_col = 6,
           silent = F)

  
  
# lets do a lil renaming

renamed = c(
  'could be nPPa'= 'What',
  'could be nMMp' = 'nMMp',
  'could be nPMP' = 'nPMp',
  'dorsolateral nucleus' = 'Dln-r',
  'dorsolateral nucleus caudal part' = 'Dln-c',
  'Dl-v2 dorsal part' = 'Dl-v2d',
  'TGn' = 'TGN'
)

renaming_col = ifelse(vis_neuron$anatomical %in% c(names(renamed)), renamed[vis_neuron$anatomical], vis_neuron$anatomical)  
vis_neuron$anatomical_renamed = renaming_col

res = 'res_2'

vals <- vis_neuron@meta.data[[res]]

  cells_use <- rownames(vis_neuron@meta.data)[!is.na(vals)]
regions <- table(vis_neuron@meta.data[cells_use, res],
                   vis_neuron@meta.data[cells_use, "anatomical_renamed"]) %>%
    as.data.frame.matrix()
  
  region_prop <- t(regions) / rowSums(t(regions))
  region_prop[is.nan(region_prop)] <- 0
  
  top_cluster <- apply(region_prop, 1, which.max)
  top_anatomy <- apply(region_prop, 2, which.max)
  col_order   <- order(top_anatomy)
  row_order   <- order(match(top_cluster, col_order))
  
  pheatmap(region_prop[row_order, col_order],
           cluster_rows = FALSE,
           cluster_cols = FALSE,
           color  = colorRampPalette(c("grey95", "orange", "darkred"))(100),
           main   = paste("Resolution:", gsub("res_", "", res)),
           show_colnames = TRUE,
           fontsize_row = 5,
           fontsize_col = 6,
           silent = F)
  
  Idents(vis_neuron)= 'res_2'
  
  
  #RUN SCT
  options(future.globals.maxSize = Inf)
library(glmGamPoi)
  rm(vis)
  vis_neuron= SCTransform(vis_neuron, assay = 'Spatial.Polygons',
                          conserve.memory = T,
                          vst.flavor = "v2",
                          ncells = 5000)
  
  saveRDS(vis_neuron,'/Users/ggraham/Desktop/Visium/vis_neuron.rds' )
  


