library(Seurat)

parker_2024 <- readRDS("A:/Anemonefish POA Legacy R Objects/parker_2024_realigned_ocellaris.rds")
multiome        <- readRDS("C:/Users/Gabe/Desktop/nemo.orig_harmony.integration_all_testd_clusters.rds")

multiome <- NormalizeData(multiome, assay = "RNA")
parker_2024 <- NormalizeData(parker_2024, assay = 'RNA')

parker_2024 =FindVariableFeatures(parker_2024)
multiome = FindVariableFeatures(multiome)

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

DimPlot(multiome, 
        group.by ='Predicted_Parent_Cluster',
        reduction = 'harmony_wnn.umap',
        label =T)

predictions2 <- TransferData(anchorset = parker_anchors,
                            refdata = parker_2024$child_clusters,
                            dims = 1:30)

multiome <- AddMetaData(multiome, metadata = predictions2)
multiome$Predicted_Child_Cluster = multiome$predicted.id

DimPlot(multiome, 
        group.by ='Predicted_Child_Cluster',
        reduction = 'harmony_wnn.umap',
        label =T)


child_tab = table(multiome$res0.8_50nn_40PC_45LSI,
      multiome$Predicted_Child_Cluster)%>%
  as.data.frame.matrix()

pheatmap((child_tab/rowSums(child_tab)),
         cluster.rows = T,
         cluster.cols = T)

parent_tab = table(multiome$res0.8_50nn_40PC_45LSI,
                  multiome$Predicted_Parent_Cluster)%>%
  as.data.frame.matrix()

pheatmap((parent_tab/rowSums(parent_tab)),
         cluster.rows = T,
         cluster.cols = T)


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

DimPlot(parker_2024, group.by ='Predicted_Cluster', label =T)

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




