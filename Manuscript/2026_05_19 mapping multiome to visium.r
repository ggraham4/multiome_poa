#### label transfer multiome -> visium res 2 heatmap production by gjg

options(future.globals.maxSize = Inf)
# libs
library(Seurat)
library(dplyr)
library(pheatmap)
library(ggplot2)
library(patchwork)

`%notin%` = Negate(`%in%`)

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

diag_order <- function(mat) {
  best_col   <- apply(mat, 1, which.max)
  col_order  <- order(apply(mat, 2, which.max))
  row_order  <- order(match(best_col, col_order))
  list(mat   = mat[row_order, col_order],
       rows   = rownames(mat)[row_order],
       cols   = colnames(mat)[col_order])
}


# read in data
vis_neuron <- readRDS('/Users/ggraham/Desktop/Visium/vis_neuron.rds')
DefaultAssay(vis_neuron)        <- "SCT"

# remove clusters that confuse me, otherwise aggregate ones that look similar
vis_neuron = subset(vis_neuron, anatomical_renamed %notin% 
                      c('ventral diffuse',
                        'Diffuse',
                        'not sure 2',
                        'not sure 4',
                        'not sure 3',
                        'not sure 1',
                        'What'
                        ))


# read in multiome object, downsample, and remove glia
obj        <- readRDS("~/Desktop/optimal_clustering_rna_only.rds")
DefaultAssay(obj)        <- "SCT"
Idents(obj) <- "final_clusters"

cells_per_cluster <- 200

set.seed(42)
cells_keep <- obj@meta.data %>%
  tibble::rownames_to_column("cell") %>%
  group_by(final_clusters) %>%
  group_modify(~ slice_sample(.x, n = min(cells_per_cluster, nrow(.x)))) %>%
  pull(cell)



obj_balanced <- subset(obj, cells = cells_keep)
non_neuronal <- c(1,  #rgc
                  2, #ol
                  11, # mg 
                  15, #ec 
                  20, # leuko 
                  13,  # opc
                  26, #dg
                 9 #poorly defined
                   )

obj_balanced_neuron = subset(obj_balanced, final_clusters %notin% non_neuronal)

# Find Transfer anchors using SCT assay and transfer data
DefaultAssay(obj)        <- "SCT"
Idents(obj) <- "final_clusters"

rm(obj)
rm(obj_balanced)

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


### vis_cluster to region heatmap ####

agreement_vis_clusters  =table(vis_neuron$anatomical_renamed,
                               vis_neuron$res_2)


mat <- diag_order(((agreement_vis_clusters) / rowSums(agreement_vis_clusters)))$mat %>%
  na.omit() %>%
  as.data.frame.matrix()

# keep clusters (cols) where at least one region maps to them strongly
col_max <- apply(mat, 2, max)
mat <- mat[, col_max > 00]  # at least one region is >15% in this cluster

# keep regions (rows) where they map strongly to at least one cluster
row_max <- apply(mat, 1, max)
mat <- mat[row_max > 0.0, ]  # region puts >15% of its cells in one cluster

vis_to_region = pheatmap(mat,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         border_color = NA,
         color = colorRampPalette(c("grey95", "orange", "red", "darkred"))(100))  

#ggsave(plot = vis_to_region,
 #      file = "vis_to_region.svg",
  #     device = "svg",
   #    units = "in",
    #   width = 3.5,
     #  height = 3.5,
      # path = "Manuscript/Plots/Manuscript v1.2.1/visium/")

#### multiome to region ####

agreement_mul  =table(vis_neuron$anatomical_renamed,
                               vis_neuron$predicted.multiome)


mul <- diag_order(((agreement_mul) / rowSums(agreement_mul)))$mat %>%
  as.data.frame.matrix()%>%
    na.omit()
  
mul = mul[-which(rownames(mul)==('not sure 21')),]

# keep clusters (cols) where at least one region maps to them strongly
col_max <- apply(mul, 2, max)
mul <- mul[, col_max > 0.0]  # at least one region is >15% in this cluster

# keep regions (rows) where they map strongly to at least one cluster
row_max <- apply(mul, 1, max)
mul <- mul[row_max > 0.0, ]  # region puts >15% of its cells in one cluster

mul = mul[rownames(mat),]

mul_to_region = pheatmap(mul[,order(colnames(mul)%>%as.numeric())],
         cluster_rows = F,
         cluster_cols = FALSE,
         border_color = NA,
         color = colorRampPalette(c("grey95", "orange", "red", "darkred"))(100),
         treeheight_row = 0)  


##ggsave(plot = mul_to_region,
  #     file = "mul_to_region.svg",
   #    device = "svg",
    #   units = "in",
     #  width = 4.8,
      # height = 4,
       #path = "Manuscript/Plots/Manuscript v1.2.1/visium/")

#mul to vis

agreement_mul_vis  =table(vis_neuron$res_2,
                               vis_neuron$predicted.multiome)


mul_vis <- diag_order(((agreement_mul_vis) / rowSums(agreement_mul_vis)))$mat %>%
  as.data.frame.matrix()%>%
    na.omit()
  

# keep clusters (cols) where at least one region maps to them strongly
col_max <- apply(mul_vis, 2, max)
mul_vis <- mul_vis[, col_max > 0.0]  # at least one region is >15% in this cluster

# keep regions (rows) where they map strongly to at least one cluster
row_max <- apply(mul_vis, 1, max)
mul_vis <- mul_vis[row_max > 0.0, ]  # region puts >15% of its cells in one cluster


mul_to_vis = pheatmap(mul_vis[order(rownames(mul_vis)%>%as.numeric()),order(colnames(mul_vis)%>%as.numeric())],
         cluster_rows = F,
         cluster_cols = FALSE,
         border_color = NA,
         color = colorRampPalette(c("grey95", "orange", "red", "darkred"))(100),
         treeheight_row = 0)  


#ggsave(plot = mul_to_vis,
#       file = "mul_to_vis.svg",
#       device = "svg",
#       units = "in",
#       width = 3.5,
#       height = 3.5,
#       path = "Manuscript/Plots/Manuscript v1.2.1/visium/")





##### Spatial Plot
clusters = c(5, 6, 7, 10,14,22,23)

#9 10 14 15 17 18 26
colors <- c('#CAE0AB', '#1965B0', '#7BAFDE', '#4EB265','#882E72' ,'#EE8026','#DC050C')
named_colors <- setNames(colors[1:length(clusters)], clusters)

p =SpatialDimPlot(vis_neuron%>%
                 subset(predicted.multiome%in%clusters), group.by = 'predicted.multiome',
               pt.size.factor = 0.5,
               images = 's_6P17.polygons',
               cols = named_colors)

#ggsave(plot = p,
#       file = "6p17_clusters.tiff",
#       device = "tiff",
#       units = "in",
#       width = 5,
#       height = 5,
#       path = "Manuscript/Plots/Manuscript v1.2.1/visium/")






