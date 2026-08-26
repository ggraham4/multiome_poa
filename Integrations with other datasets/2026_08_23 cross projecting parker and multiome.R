library(Seurat)

parker_2024 <- readRDS("A:/Anemonefish POA Legacy R Objects/parker_2024_realigned_ocellaris.rds")
multiome        <- readRDS("C:/Users/Gabe/Desktop/nemo.orig_harmony.integration_all_testd_clusters.rds")

multiome <- NormalizeData(multiome, assay = "RNA")
parker_2024 <- NormalizeData(parker_2024, assay = 'RNA')

parker_2024 =FindVariableFeatures(parker_2024)
multiome = FindVariableFeatures(multiome)

my_colors <- c("gray", 'white', "orange", "red")
color_func <- colorRampPalette(my_colors)
continuous_palette <- color_func(50)


parker_anchors = FindTransferAnchors(reference = parker_2024,
                                     query = multiome,
                                     dims = 1:30,
                                     reference.reduction = "pca",
                                     features = (rownames(parker_2024)))

predictions <- TransferData(anchorset = parker_anchors,
                            refdata = parker_2024$parent_clusters,
                            dims = 1:30)

multiome <- AddMetaData(multiome, metadata = predictions)
multiome$Predicted_Parent_Cluster = multiome$predicted.id

parent_on_multiome=DimPlot(multiome, 
        group.by ='Predicted_Parent_Cluster',
        reduction = 'harmony_wnn.umap',
        label =T,
        raster =T,
        pt.size = 6,
        raster.dpi = c(2400,2400))
parent_on_multiome

ggsave(plot = parent_on_multiome,
       file = paste0('parent_on_multiome.svg'),
       device = "svg",
       units = "in",
       width =8,
       height = 8,
       path = "Manuscript/Plots")

multiome_native=DimPlot(multiome, 
                           group.by ='res0.8_50nn_40PC_45LSI',
                           reduction = 'harmony_wnn.umap',
                           label =T,
                           raster =T,
                           pt.size = 6,
                           raster.dpi = c(2400,2400))
multiome_native

ggsave(plot = multiome_native,
       file = paste0('multiome_native.svg'),
       device = "svg",
       units = "in",
       width =8,
       height = 8,
       path = "Manuscript/Plots")


predictions2 <- TransferData(anchorset = parker_anchors,
                            refdata = parker_2024$child_clusters,
                            dims = 1:30)

multiome <- AddMetaData(multiome, metadata = predictions2)
multiome$Predicted_Child_Cluster = multiome$predicted.id

library(colorspace)

generate_nonadjacent_palette <- function(n, seed = 1) {
  set.seed(seed)
  
  # 1. Base colors: full hue circle, cycling chroma/luminance so
  #    it's not just hue doing the work of differentiating clusters
  hues    <- seq(15, 375, length.out = n + 1)[1:n]
  lum_cyc <- rep(c(60, 45, 80, 55), length.out = n)
  chr_cyc <- rep(c(100, 75, 55, 90), length.out = n)
  base_colors <- hcl(h = hues, c = chr_cyc, l = lum_cyc)
  
  # 2. Coprime-stride permutation: pick k with gcd(k, n) == 1
  gcd <- function(a, b) if (b == 0) a else gcd(b, a %% b)
  k <- floor(n / 2)
  while (gcd(k, n) != 1) k <- k - 1
  
  perm <- ((seq(0, n - 1) * k) %% n) + 1
  base_colors[perm]
}

my_49_colors <- generate_nonadjacent_palette(49)

child_on_multiome=DimPlot(multiome, 
        group.by ='Predicted_Child_Cluster',
        reduction = 'harmony_wnn.umap',
        label =T,
        raster =T,
        pt.size = 6,
        raster.dpi = c(2400,2400),
        cols = my_49_colors)

child_on_multiome

ggsave(plot = child_on_multiome,
       file = paste0('child_on_multiome.svg'),
       device = "svg",
       units = "in",
       width =8,
       height = 8,
       path = "Manuscript/Plots")


child_tab = table(multiome$res0.8_50nn_40PC_45LSI,
      multiome$Predicted_Child_Cluster)%>%
  as.data.frame.matrix()

child_to_multiome_prop_heatmap =pheatmap((child_tab/rowSums(child_tab)),
         cluster.rows = T,
         cluster.cols = T,
         color = continuous_palette,
         angle_col = 90
)


svg(filename = 'C:/Users/Gabe/Desktop/multiome_poa/Manuscript/Plots/child_to_multiome_prop_heatmap.svg',
    width = 8,
    height = 6.5)
child_to_multiome_prop_heatmap
dev.off()

parent_tab = table(multiome$res0.8_50nn_40PC_45LSI,
                  multiome$Predicted_Parent_Cluster)%>%
  as.data.frame.matrix()

parent_to_multiome_prop_heatmap =pheatmap((parent_tab/rowSums(parent_tab)),
                                          cluster.rows = T,
                                          cluster.cols = T,
                                          color = continuous_palette,
                                          angle_col = 90
)

svg(filename = 'C:/Users/Gabe/Desktop/multiome_poa/Manuscript/Plots/parent_to_multiome_prop_heatmap.svg',
    width = 8,
    height = 6.5)
parent_to_multiome_prop_heatmap
dev.off()


### other way around ###
multiome =RunPCA(multiome)
multiome_anchors = FindTransferAnchors(reference = multiome,
                                     query = parker_2024,
                                     dims = 1:30,
                                     reference.reduction = "pca",
                                     features = (rownames(multiome)))


predictions3 <- TransferData(anchorset = multiome_anchors,
                            refdata = multiome$res0.8_50nn_40PC_45LSI,
                            dims = 1:30)

parker_2024 <- AddMetaData(parker_2024, metadata = predictions3)
parker_2024$Predicted_Cluster = parker_2024$predicted.id

multiome_on_parker = DimPlot(parker_2024,
                             group.by ='Predicted_Cluster',
                             label =T,
                             raster =T,
                             pt.size = 6,
                             raster.dpi = c(2400,2400))
multiome_on_parker

ggsave(plot = multiome_on_parker,
       file = paste0('multiome_on_parker.svg'),
       device = "svg",
       units = "in",
       width =8,
       height = 8,
       path = "Manuscript/Plots")

parker_parent = DimPlot(parker_2024,
                             group.by ='parent_clusters',
                             label =T,
                             raster =T,
                             pt.size = 6,
                             raster.dpi = c(2400,2400))
parker_parent

ggsave(plot = parker_parent,
       file = paste0('parker_parent.svg'),
       device = "svg",
       units = "in",
       width =8,
       height = 8,
       path = "Manuscript/Plots")

parker_child = DimPlot(parker_2024,
                        group.by ='child_clusters',
                        label =T,
                        raster =T,
                        pt.size = 6,
                        raster.dpi = c(2400,2400))
parker_child

ggsave(plot = parker_child,
       file = paste0('parker_child.svg'),
       device = "svg",
       units = "in",
       width =8,
       height = 8,
       path = "Manuscript/Plots")



mult_tab = table(multiome$Predicted_Parent_Cluster,
                  multiome$res0.8_50nn_40PC_45LSI)%>%
  as.data.frame.matrix()

mult_tab_prop_heatmap =pheatmap((mult_tab/rowSums(mult_tab)),
                                cluster.rows = T,
                                cluster.cols = T,
                                color = continuous_palette,
                                angle_col = 90
)
mult_tab_prop_heatmap

svg(filename = 'C:/Users/Gabe/Desktop/multiome_poa/Manuscript/Plots/multiome_to_parent_prop_heatmap.svg',
    width = 8,
    height = 6.5)
mult_tab_prop_heatmap
dev.off()

mult_tab_ch = table(multiome$Predicted_Child_Cluster,
                 multiome$res0.8_50nn_40PC_45LSI)%>%
  as.data.frame.matrix()

mult_tab_child_prop_heatmap =pheatmap((mult_tab_ch/rowSums(mult_tab_ch)),
                                      cluster.rows = T,
                                      cluster.cols = T,
                                      color = continuous_palette,
                                      angle_col = 90
)
mult_tab_child_prop_heatmap

svg(filename = 'C:/Users/Gabe/Desktop/multiome_poa/Manuscript/Plots/multiome_to_child_prop_heatmap.svg',
    width = 8,
    height = 6.5)
mult_tab_child_prop_heatmap
dev.off()


### ============================================================
### Mutual best-match (reciprocal best hit) cluster analysis
### Continues from the multiome <-> parker_2024 label-transfer
### code you already ran.
### ============================================================

library(dplyr)

## ---- 1. Helper: build a row-normalized proportion table ------
## rows   = "query" cluster identity (character/factor)
## cols   = predicted identity from the other dataset
prop_table <- function(cluster_ids, predicted_ids) {
  tab <- table(cluster_ids, predicted_ids) %>% as.data.frame.matrix()
  tab / rowSums(tab)
}

## ---- 2. Helper: best match per row (which column has max proportion) ----
best_match <- function(prop_tab) {
  data.frame(
    cluster      = rownames(prop_tab),
    best_match   = colnames(prop_tab)[apply(prop_tab, 1, which.max)],
    best_match_prop = apply(prop_tab, 1, max),
    row.names = NULL
  )
}

## ---- 3. Helper: intersect two best-match tables to get mutual best hits ----
## dirA: data.frame(cluster = A_id, best_match = B_id, best_match_prop)
## dirB: data.frame(cluster = B_id, best_match = A_id, best_match_prop)
mutual_best_matches <- function(dirA, dirB) {
  colnames(dirA) <- c("A_cluster", "A_best_in_B", "A_prop")
  colnames(dirB) <- c("B_cluster", "B_best_in_A", "B_prop")
  
  merged <- merge(dirA, dirB,
                  by.x = "A_best_in_B", by.y = "B_cluster",
                  all = FALSE)
  
  merged %>%
    filter(A_cluster == B_best_in_A) %>%   # reciprocal condition
    transmute(
      multiome_cluster = A_cluster,
      parker_cluster    = A_best_in_B,
      multiome_to_parker_prop = A_prop,
      parker_to_multiome_prop = B_prop
    )
}

### ============================================================
### CHILD-LEVEL mutual best matches
### ============================================================

## Direction 1: multiome cluster -> best-matching parker child cluster
## (you already built this un-normalized as `child_tab`; rebuild normalized)
child_prop_A <- prop_table(multiome$res0.8_50nn_40PC_45LSI,
                           multiome$Predicted_Child_Cluster)
child_bestA  <- best_match(child_prop_A)

## Direction 2: parker child cluster -> best-matching multiome cluster
## Need the reverse contingency table using predictions3 results,
## but predictions3 was run predicting multiome clusters onto parker
## cells based on the *overall* PCA transfer -- reuse parker_2024$Predicted_Cluster
child_prop_B <- prop_table(parker_2024$child_clusters,
                           parker_2024$Predicted_Cluster)
child_bestB  <- best_match(child_prop_B)

child_mutual <- mutual_best_matches(child_bestA, child_bestB)
print(child_mutual)

### ============================================================
### PARENT-LEVEL mutual best matches
### ============================================================

parent_prop_A <- prop_table(multiome$res0.8_50nn_40PC_45LSI,
                            multiome$Predicted_Parent_Cluster)
parent_bestA  <- best_match(parent_prop_A)

parent_prop_B <- prop_table(parker_2024$parent_clusters,
                            parker_2024$Predicted_Cluster)
parent_bestB  <- best_match(parent_prop_B)

parent_mutual <- mutual_best_matches(parent_bestA, parent_bestB)
print(parent_mutual)


#### Export #### 
write.csv(parent_mutual, "Manuscript/Manuscript v.2/Supplemental Files/Appendix SX. Multiome Parker Parent Clusters Mutual Best Hits.csv")
write.csv(child_mutual, "Manuscript/Manuscript v.2/Supplemental Files/Appendix SX. Multiome Parker Child Clusters Mutual Best Hits.csv")

