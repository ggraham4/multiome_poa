# =============================================================================
# LABEL TRANSFER — two approaches
# =============================================================================
# =============================================================================
options(future.globals.maxSize = Inf)

library(Seurat)
library(dplyr)
library(pheatmap)
library(ggplot2)
library(patchwork)
# =============================================================================
# APPROACH 1: Formal FindTransferAnchors
# multiome obj (reference) → vis_neuron (query)
# =============================================================================
#obj        <- readRDS("~/Desktop/optimal_clustering_rna_only.rds")
#Idents(obj) <- "final_clusters"

vis_neuron <- readRDS('/Users/ggraham/Desktop/Visium/vis_neuron.rds')

# ==================

#DefaultAssay(vis_neuron) <- "SCT"
#DefaultAssay(obj) <- "RNA"
#obj <- SCTransform(obj, verbose = FALSE)

#DefaultAssay(obj)        <- "SCT"

#obj= RunPCA(obj)
#DefaultAssay(obj)        <- "RNA"
#obj$Status = factor(obj$Status , levels = c('NRM', 'M', "D",'E','NF','F'))

#saveRDS(obj,"~/Desktop/optimal_clustering_rna_only.rds")


obj        <- readRDS("~/Desktop/optimal_clustering_rna_only.rds")
DefaultAssay(obj)        <- "SCT"
Idents(obj) <- "final_clusters"

anchors <- FindTransferAnchors(
  reference          = obj,
  query              = vis_neuron,
  normalization.method = "SCT",
  reference.reduction  = "pca",
  reduction = 'pcaproject',
  dims               = 1:30
)

vis_neuron <- TransferData(
  anchorset       = anchors,
  query           = vis_neuron,
  reference       = obj,
  refdata         = list(
    transferred_cluster = "final_clusters"
  ),
  weight.reduction = "pcaproject",
  dims             = 1:30
)

# now vis_neuron has:
# vis_neuron$predicted.transferred_cluster       — winning label
# vis_neuron$predicted.transferred_cluster.score — confidence 0-1

# check score distribution — low scores mean ambiguous assignments
hist(vis_neuron$predicted.transferred_cluster.score, breaks = 50,
     main = "Label transfer confidence", xlab = "Score")

# heatmap: transferred cluster vs anatomical region
regions_lt <- table(vis_neuron$predicted.transferred_cluster,
                    vis_neuron$anatomical_renamed) %>%
  as.data.frame.matrix()

region_prop_lt <- t(regions_lt) / rowSums(t(regions_lt))

diag_order <- function(mat) {
  best_col   <- apply(mat, 1, which.max)
  col_order  <- order(apply(mat, 2, which.max))
  row_order  <- order(match(best_col, col_order))
  list(mat   = mat[row_order, col_order],
       rows   = rownames(mat)[row_order],
       cols   = colnames(mat)[col_order])
}

plot_cor <- function(mat, title, ...) {
  ord <- diag_order(mat)
  pheatmap(
    scale(t(scale(t(ord$mat)))),
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    border_color = NA,
    color        = colorRampPalette(c("blue", "white", "darkred"))(100),
    fontsize_row = 7,
    fontsize_col = 7,
    main         = title,
    ...
  )
}

ord <- diag_order(region_prop_lt)

pheatmap(
  ord$mat,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  border_color = NA,
  color        = colorRampPalette(c("grey95", "orange", "darkred"))(100),
  fontsize_row = 7,
  fontsize_col = 7,
  main         = "Approach 1: Transferred cluster vs anatomical region"
)

# spatial plot — where does each cluster land on the tissue?
#SpatialDimPlot(vis_neuron, 
 #              group.by = "predicted.transferred_cluster",
  #             label    = FALSE)

# =============================================================================
# APPROACH 2: seurat_clusters.projected0.4 already in vis_neuron
# this is the sketch-based projection Seurat did internally
# =============================================================================

# check it exists
table(vis_neuron$seurat_clusters.projected0.4)

regions_proj <- table(vis_neuron$seurat_clusters.projected0.4,
                      vis_neuron$anatomical_renamed) %>%
  as.data.frame.matrix()

region_prop_proj <- t(regions_proj) / rowSums(t(regions_proj))

ord2 <- diag_order(region_prop_proj)

pheatmap(
  ord2$mat,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  border_color = NA,
  color        = colorRampPalette(c("grey95", "orange", "darkred"))(100),
  fontsize_row = 7,
  fontsize_col = 7,
  main         = "Approach 2: Projected 0.4 clusters vs anatomical region"
)

#SpatialDimPlot(vis_neuron,
 #              group.by = "seurat_clusters.projected0.4",
  #             label    = FALSE)

# =============================================================================
# COMPARE THE TWO APPROACHES
# do they agree on cell assignments?
# =============================================================================

agreement <- table(
  vis_neuron$predicted.transferred_cluster,
  vis_neuron$seurat_clusters.projected0.4
)

# what fraction of cells get concordant top assignments?
# map each projected cluster to its best transferred cluster
best_transfer <- apply(agreement, 2, which.max)
cat("Cluster agreement between approaches:\n")
print(best_transfer)

# score distribution split by whether assignment was high confidence
vis_neuron$transfer_confident <- 
  vis_neuron$predicted.transferred_cluster.score > 0.5

cat("\nHigh confidence assignments:", 
    sum(vis_neuron$transfer_confident), "\n")
cat("Low confidence assignments:", 
    sum(!vis_neuron$transfer_confident), "\n")


# projected 0.4 clusters vs transferred multiome clusters
agreement_prop <- t(agreement) / rowSums(t(agreement))

ord3 <- diag_order(agreement_prop)

pheatmap(
  ord3$mat,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  border_color = NA,
  color        = colorRampPalette(c("grey95", "orange", "darkred"))(100),
  fontsize_row = 7,
  fontsize_col = 7,
  main         = "Projected 0.4 clusters vs Transferred multiome clusters"
)


### this sucks, downsample ###
# downsample obj to equal cells per cluster
cells_per_cluster <- 200  # adjust based on your smallest cluster size

set.seed(42)
cells_keep <- obj@meta.data %>%
  tibble::rownames_to_column("cell") %>%
  group_by(final_clusters) %>%
  group_modify(~ slice_sample(.x, n = min(cells_per_cluster, nrow(.x)))) %>%
  pull(cell)

`%notin%`= Negate(`%in%`)

obj_balanced <- subset(obj, cells = cells_keep)
non_neuronal <- c(1,  #rgc
                  2, #ol
                  11, # mg 
                  15, #ec 
                  20, # leuko 
                  13,  # opc
                  26, #dg
                  # remove 9 cause its a dick
                  9)

obj_balanced_neuron = subset(obj_balanced, final_clusters %notin% non_neuronal)

# check balance
table(obj_balanced$final_clusters)

# now run transfer on balanced reference
anchors <- FindTransferAnchors(
  reference            = obj_balanced_neuron,
  query                = vis_neuron,
  normalization.method = "SCT",
  reference.reduction  = "pca",
  reduction            = "pcaproject",
  dims                 = 1:30
)

vis_neuron <- TransferData(
  anchorset        = anchors,
  query            = vis_neuron,
  reference        = obj_balanced_neuron,
  refdata          = list(transferred_cluster = "final_clusters"),
  weight.reduction = "pcaproject",
  dims             = 1:30
)

agreement <- table(
  vis_neuron$predicted.transferred_cluster,
  vis_neuron$res_2
)

# what fraction of cells get concordant top assignments?
# map each projected cluster to its best transferred cluster
best_transfer <- apply(agreement, 2, which.max)
cat("Cluster agreement between approaches:\n")
print(best_transfer)

# score distribution split by whether assignment was high confidence
vis_neuron$transfer_confident <- 
  vis_neuron$predicted.transferred_cluster.score > 0.3

cat("\nHigh confidence assignments:", 
    sum(vis_neuron$transfer_confident), "\n")
cat("Low confidence assignments:", 
    sum(!vis_neuron$transfer_confident), "\n")


  vis_neuron$predicted.multiome =   vis_neuron$predicted.transferred_cluster

pheatmap(
  t(diag_order((agreement)/rowSums(agreement))$mat),
  cluster_rows = F,
  cluster_cols = F,
  border_color = NA,
  color        = colorRampPalette(c("grey95", "orange", "darkred"))(100)
)

agreement_vis_clusters  =table(vis_neuron$anatomical_renamed,
                               vis_neuron$res_2)

pheatmap(
  diag_order(((agreement_vis_clusters)/rowSums(agreement_vis_clusters)))$mat%>%na.omit(),
  cluster_rows = T,
  cluster_cols = T,
  border_color = NA,
  color        = colorRampPalette(c("grey95", "orange",'red', "darkred"))(100)
)



#saveRDS(vis_neuron,'/Users/ggraham/Desktop/Visium/vis_neuron.rds' )

library(cowplot)
x11(width = 16)
SpatialDimPlot(vis_neuron%>%
                 subset(predicted.multiome==6), group.by = 'predicted.multiome',
               pt.size.factor = 0.5,
               images = 's_6P17.polygons')

x11(width = 16)
SpatialDimPlot(vis_neuron%>%
                 subset(predicted.multiome==0), group.by = 'predicted.multiome',
               pt.size.factor = 0.5,
               images = 's_6P17.polygons')

x11(width = 16)
SpatialDimPlot(vis_neuron%>%
                 subset(predicted.multiome==3), group.by = 'predicted.multiome',
               pt.size.factor = 0.5,
               images = 's_6P17.polygons')

x11(width = 16)
SpatialDimPlot(vis_neuron%>%
                 subset(predicted.multiome==24), group.by = 'predicted.multiome',
               pt.size.factor = 0.5,
               images = 's_6P17.polygons')

x11(width = 16)
SpatialDimPlot(vis_neuron%>%
                 subset(predicted.multiome==7), group.by = 'predicted.multiome',
               pt.size.factor = 0.5,
               images = 's_6P17.polygons')

x11(width = 16)
SpatialDimPlot(vis_neuron%>%
                 subset(predicted.multiome==8), group.by = 'predicted.multiome',
               pt.size.factor = 0.5,
               images = 's_6P17.polygons')

library(stringr)

Slices = unique(vis_neuron$Slice)
for(i in Slices){
  meta_slice = subset(vis_neuron@meta.data, Slice== i)
  
multiome_predictions = meta_slice[,'predicted.multiome']%>%as.data.frame()
multiome_predictions$Barcode = rownames(meta_slice)
#fix barcodes
multiome_predictions$Barcode <- str_sub(multiome_predictions$Barcode, end = -3)
multiome_predictions = multiome_predictions[,c(2,1)]

colnames(multiome_predictions)= c('Barcode', 'multiome_predicted')
write.csv(multiome_predictions, paste0('Visium/Groupings/',i,' multiome predicted.csv'), row.names = FALSE)
}

Slices = unique(vis_neuron$Slice)
for(i in Slices){
  meta_slice = subset(vis_neuron@meta.data, Slice== i)
  
res_2 = meta_slice[,'res_2']%>%as.data.frame()
res_2$Barcode = rownames(meta_slice)
#fix barcodes
res_2$Barcode <- str_sub(res_2$Barcode, end = -3)
res_2 = res_2[,c(2,1)]

colnames(res_2)= c('Barcode', 'res_2')
write.csv(res_2, paste0('Visium/Groupings/',i,'res_2.csv'), row.names = FALSE)
}








