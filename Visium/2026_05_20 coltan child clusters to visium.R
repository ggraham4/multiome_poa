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


# read in parker object, downsample, and remove glia
obj        <- readRDS("~/Desktop/parker_object.rds")
obj = SCTransform(obj)
DefaultAssay(obj)        <- "SCT"

Idents(obj) <- "clusters49"

cells_per_cluster <- 200

set.seed(42)
cells_keep <- obj@meta.data %>%
  tibble::rownames_to_column("cell") %>%
  group_by(clusters49) %>%
  group_modify(~ slice_sample(.x, n = min(cells_per_cluster, nrow(.x)))) %>%
  pull(cell)



obj_balanced <- subset(obj, cells = cells_keep)
non_neuronal <- c(19, # problematic neurons
                  20,
                  45,
                  22, # oligo
                  40, # leuko
                  44, 15, 3, 32,46, #rgc,
                  28 # opc
                  )

obj_balanced_neuron = subset(obj_balanced, clusters49 %notin% non_neuronal)

# Find Transfer anchors using SCT assay and transfer data
DefaultAssay(obj)        <- "SCT"
Idents(obj) <- "clusters49"

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


vis_neuron_sub <- TransferData(
  anchorset        = anchors,
  query            = vis_neuron,
  reference        = obj_balanced_neuron,
  refdata          = list(transferred_cluster = "clusters49"),
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
vis_neuron_sub$transfer_confident <- 
  vis_neuron_sub$predicted.transferred_cluster.score > 0.3

cat("\nHigh confidence assignments:", 
    sum(vis_neuron_sub$transfer_confident), "\n")
cat("Low confidence assignments:", 
    sum(!vis_neuron_sub$transfer_confident), "\n")


  vis_neuron_sub$predicted.multiome =   vis_neuron_sub$predicted.transferred_cluster

pheatmap(
  t(diag_order((agreement)/rowSums(agreement))$mat),
  cluster_rows = F,
  cluster_cols = F,
  border_color = NA,
  color        = colorRampPalette(c("grey95", "orange", "darkred"))(100)
)


agreement_mul  =table(vis_neuron_sub$anatomical_renamed,
                               vis_neuron_sub$predicted.multiome)


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


#mul_to_region = pheatmap(mul[,order(colnames(mul)%>%as.numeric())],
#         cluster_rows = T,
#         cluster_cols = FALSE,
#         border_color = NA,
 #        color = colorRampPalette(c("grey95", "orange", "red", "darkred"))(100),
#         treeheight_row = 0, 
#         main = 'Percent of Each Region Child Clusters')  


ggsave(plot = mul_to_region,
       file = "pakrer_child_to_region.svg",
       device = "svg",
       units = "in",
       width = 3.5,
       height = 3.5,
       path = "Manuscript/Plots/Manuscript v1.2.1/visium/")


mul_sub <- mul

cluster_top_region <- apply(mul_sub, 2, function(x) rownames(mul_sub)[which.max(x)])

obj_balanced_neuron@meta.data$top_anatomical_region <- 
  cluster_top_region[as.character(obj_balanced_neuron@meta.data$clusters49)]

DimPlot(obj_balanced_neuron, group.by = "top_anatomical_region", label = TRUE) +
  ggtitle("Top anatomical region per cluster")

# second highest hit per cluster
cluster_second_region <- apply(mul_sub, 2, function(x) {
  sorted <- sort(x, decreasing = TRUE)
  rownames(mul_sub)[which(x == sorted[2])[1]]
})

obj_balanced_neuron@meta.data$second_anatomical_region <- 
  cluster_second_region[as.character(obj_balanced_neuron@meta.data$clusters49)]

DimPlot(obj_balanced_neuron, group.by = "second_anatomical_region", label = TRUE) +
  ggtitle("Second anatomical region per  cluster")


clusters = c(25, # poa
             2 #dm
             )

#9 10 14 15 17 18 26
colors <- c('#1965B0','#DC050C')
named_colors <- setNames(colors[1:length(clusters)], clusters)

p =SpatialDimPlot(vis_neuron_sub%>%
                 subset(predicted.multiome%in%clusters), group.by = 'predicted.multiome',
               pt.size.factor = 0.5,
               images = 's_6P17.polygons',
               cols = named_colors)

p

ggsave(plot = p,
       file = "6p17_coltan_poa_dm3_clusters.tiff",
       device = "tiff",
       units = "in",
       width = 5,
       height = 5,
       path = "Manuscript/Plots/Manuscript v1.2.1/visium/")
