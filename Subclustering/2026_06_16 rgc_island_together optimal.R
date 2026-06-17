library(Seurat)
library(Signac)
library(ggplot2)
library(tidyverse)
library(lme4)
library(cluster)
library(bluster)

obj <- readRDS("~/Desktop/nemo.orig_harmony.integration_all_testd_clusters.rds")

sub_obj <- subset(obj, res0.8_50nn_40PC_45LSI %in% c(0,1,6,9,3,5,24,18,21)& 
                    Status != 'NRM')

phase_to_status = c('M'= 'M',
                    'D'='I',
                    'E' = 'LI',
                    'NF'= 'NF',
                    'F'='F')
sub_obj$Phase = phase_to_status[sub_obj$Status]%>%unname()
sub_obj$Phase = factor(sub_obj$Phase, levels = c('M',
                                                 'I',
                                                 'LI',
                                                 'NF',
                                                 'F'))
###############################################################################
# RNA preprocessing
###############################################################################

DefaultAssay(sub_obj) <- "RNA"

sub_obj <- NormalizeData(sub_obj)
sub_obj <- FindVariableFeatures(sub_obj, selection.method = "vst", nfeatures = 2000)
sub_obj <- ScaleData(sub_obj)

sub_obj <- RunPCA(
  sub_obj,
  features = VariableFeatures(sub_obj),
  reduction.name = "pca"
)

###############################################################################
# ATAC preprocessing
###############################################################################

DefaultAssay(sub_obj) <- "ATAC"

sub_obj <- RunTFIDF(sub_obj)
sub_obj <- FindTopFeatures(sub_obj, min.cutoff = "q0")
sub_obj <- RunSVD(sub_obj, reduction.name = "lsi")

###############################################################################
# Data-driven dim selection
###############################################################################

# RNA: pick PCs explaining >=90% cumulative variance
pct_var <- sub_obj[["pca"]]@stdev^2 / sum(sub_obj[["pca"]]@stdev^2)
n_pcs   <- min(which(cumsum(pct_var) >= 0.90))

# ATAC: extract correlation values from DepthCor plot object
dc_plot <- DepthCor(sub_obj, reduction = "lsi", n = 30)
dc_df   <- dc_plot$data                          # pull the underlying data frame
bad_lsi <- which(abs(dc_df$counts) > 0.75)       # 'counts' col holds correlations
lsi_use <- setdiff(2:30, bad_lsi)                # always skip dim 1

cat("RNA PCs:", n_pcs, "| LSI dims:", paste(range(lsi_use), collapse = "-"), "\n")

###############################################################################
# Dim-combo sweep
###############################################################################

pc_grid  <- unique(c(10, 20, n_pcs))
lsi_grid <- list(2:15, 2:20, lsi_use)

dim_results <- lapply(pc_grid, function(npc) {
  lapply(lsi_grid, function(ld) {
    tag <- paste0("pc", npc, "_lsi", max(ld))
    tmp <- FindMultiModalNeighbors(sub_obj,
             reduction.list       = list("pca", "lsi"),
             dims.list            = list(1:npc, ld),
             modality.weight.name = "RNA.weight")
    tmp <- FindClusters(tmp, graph.name = "wsnn",
                        algorithm = 3, resolution = 0.2, verbose = FALSE)
    tmp <- RunUMAP(tmp, nn.name = "weighted.nn",
                   reduction.name = "wnn.umap", reduction.key = "wnnUMAP_",
                   verbose = FALSE)
    emb  <- Embeddings(tmp, "wnn.umap")
    labs <- as.integer(tmp$seurat_clusters)
    sil  <- if (length(unique(labs)) > 1)
              mean(silhouette(labs, dist(emb))[, 3]) else NA
    data.frame(tag = tag, npc = npc, lsi_max = max(ld),
               sil = sil, k = length(unique(labs)))
  })
})

dim_df   <- do.call(rbind, do.call(c, dim_results))
best_dim <- dim_df[which.max(dim_df$sil), ]
cat("Best dim combo:", best_dim$tag, "| k =", best_dim$k,
    "| sil =", round(best_dim$sil, 3), "\n")

###############################################################################
# Resolution sweep on best dim combo
###############################################################################

sub_obj <- FindMultiModalNeighbors(sub_obj,
             reduction.list       = list("pca", "lsi"),
             dims.list            = list(1:best_dim$npc, 2:best_dim$lsi_max),
             modality.weight.name = "RNA.weight")
sub_obj <- RunUMAP(sub_obj, nn.name = "weighted.nn",
                   reduction.name = "wnn.umap", reduction.key = "wnnUMAP_",
                   verbose = FALSE)
emb <- Embeddings(sub_obj, "wnn.umap")

res_grid   <- seq(0.05, 0.8, by = 0.05)
res_scores <- lapply(res_grid, function(r) {
  tmp  <- FindClusters(sub_obj, graph.name = "wsnn",
                       algorithm = 3, resolution = r, verbose = FALSE)
  labs <- as.integer(tmp$seurat_clusters)
  k    <- length(unique(labs))
  sil  <- if (k > 1 && k < nrow(emb) - 1)
            mean(silhouette(labs, dist(emb))[, 3]) else NA
  data.frame(resolution = r, k = k, mean_sil = sil)
})

res_df   <- do.call(rbind, res_scores)
best_res <- res_df$resolution[which.max(res_df$mean_sil)]

print(res_df)
cat("Best resolution:", best_res, "\n")

###############################################################################
# Final clustering
###############################################################################

sub_obj <- FindClusters(sub_obj, graph.name = "wsnn",
                        algorithm = 3, resolution = best_res)
DimPlot(sub_obj, reduction = "wnn.umap", label = TRUE)

###############################################################################
# Evaluate
###############################################################################

VlnPlot(sub_obj, features = "RNA.weight", group.by = "seurat_clusters")

purity <- neighborPurity(emb, sub_obj$seurat_clusters)
cat("Median neighbor purity:", median(purity$purity), "\n")

ggplot(res_df, aes(x = resolution, y = mean_sil)) +
  geom_line(color = "#7F77DD", linewidth = 1) +
  geom_point(aes(size = k), color = "#7F77DD") +
  geom_vline(xintercept = best_res, linetype = "dashed",
             color = "#D85A30", linewidth = 0.8) +
  scale_size_continuous(range = c(2, 6), name = "# clusters") +
  labs(x = "Resolution", y = "Mean silhouette",
       title = "Resolution sweep — dashed = selected") +
  theme_bw()


DefaultAssay(sub_obj) = 'RNA' 
Idents(sub_obj) = 'wsnn_res.0.5'
marks = FindAllMarkers(sub_obj, only.pos = T)

supertype =c(
  '0'= 'cckb+',
  '1'='gal+',
  '2'='sim1+',
  '3'= 'gnrh+'
)

cyto = CytoTRACE::CytoTRACE(sub_obj@assays$RNA$data%>%as.matrix())$CytoTRACE
sub_obj$cyto = cyto 

FeaturePlot(sub_obj, 'cyto', reduction ='wnn.umap')
FeaturePlot(sub_obj, 's100b', reduction ='wnn.umap')

FeaturePlot(sub_obj, 'gad2', reduction ='wnn.umap')
FeaturePlot(sub_obj, 'slc17a6b', reduction ='wnn.umap')
FeaturePlot(sub_obj, 'LOC111571064', reduction ='wnn.umap')
FeaturePlot(sub_obj, 'rbfox3a', reduction ='wnn.umap')

DimPlot(sub_obj, reduction   = 'wnn.umap', label = T)
DimPlot(sub_obj, reduction   = 'wnn.umap', label = T, group.by = 'res0.8_50nn_40PC_45LSI')

marx = FindAllMarkers(sub_obj, only.pos = T)

saveRDS(sub_obj, '~/Desktop/sub_obj_rgc_island_optim_clustering.rds')

