#labels parker et al 2024
{
  library(parallel)
  library(clusterProfiler)
  library(blme)
  library(Seurat)
  library(tidyverse)
  library(tidyr)
  library(lme4)
  library(dplyr)
  library(MASS)
  library(SeuratObject)
  library(Signac)
  library('glmGamPoi')
  library(scran)
  library(parallel)
  library(factoextra)
  library(readxl)
  library(factoextra)
  library(forcats)
  library(ggrepel)
  library(biomaRt)
  library(openxlsx)
  library(ComplexHeatmap)
}

multiome_object <- readRDS("C:/seurat_objects/optimal_clustering_05_06_2025.rds")
Idents(multiome_object) <- "res0.8_50nn_40PC_45LSI"


#### Find anchorsets and transfer labels #####
previous_object <- readRDS('C:/Users/Gabe/Desktop/Old object.rds')

#### Clusters 20 ####
Idents(previous_object) <- 'clusters20'
anchors <- FindTransferAnchors(reference = previous_object, query = multiome_object,
                               reference.reduction = "pca")

predictions <- TransferData(anchorset = anchors, refdata = previous_object$clusters20)
multiome_object <- AddMetaData(multiome_object, metadata = predictions)

Idents(multiome_object) <- 'predicted.id'
DimPlot(multiome_object, label = T, reduction = 'harmony_wnn.umap')

##other way around
Idents(multiome_object) <- 'res0.8_50nn_40PC_45LSI'
anchors2 <- FindTransferAnchors(reference = multiome_object, query = previous_object,
                                reference.reduction = "rnaPCA")

predictions2 <- TransferData(anchorset = anchors2, refdata = multiome_object$res0.8_50nn_40PC_45LSI)
previous_object <- AddMetaData(previous_object, metadata = predictions2)

Idents(previous_object) <- 'predicted.id'
DimPlot(previous_object, label = T)

predictions2$query_id <- previous_object$clusters20

old_obj_predictions_matrix <- predictions2%>%
  group_by(query_id)%>%
  summarise(across(starts_with("prediction.score"), mean, na.rm = TRUE))

colnames(old_obj_predictions_matrix) <- str_remove(colnames(old_obj_predictions_matrix), 'prediction.score.')

colnames(old_obj_predictions_matrix) <- as.numeric(colnames(old_obj_predictions_matrix))

old_obj_predictions_matrix%>%
  order()

old_obj_predictions_matrix[,2:(ncol(old_obj_predictions_matrix)-1)]%>%
  as.matrix()%>%
  scale()%>%
  Heatmap(
    cluster_rows = FALSE, 
    cluster_columns = FALSE,
    row_order = 1:20,
    row_labels = 0:19,
    column_order = order(as.numeric(colnames(old_obj_predictions_matrix[,2:(ncol(old_obj_predictions_matrix)-1)]))),
    row_title = 'previous',
    column_title = 'multiome')

DimPlot(previous_object, label = T)


### align to clusters49
Idents(previous_object) <- 'clusters49'
anchors_49 <- FindTransferAnchors(reference = previous_object, query = multiome_object,
                               reference.reduction = "pca")

predictions_49 <- TransferData(anchorset = anchors_49, refdata = previous_object$clusters49)
multiome_object <- AddMetaData(multiome_object, metadata = predictions_49)

Idents(multiome_object) <- 'predicted.id'
DimPlot(multiome_object, label = T, reduction = 'harmony_wnn.umap')

##other way around
Idents(multiome_object) <- 'res0.8_50nn_40PC_45LSI'
anchors_49_2 <- FindTransferAnchors(reference = multiome_object, query = previous_object,
                                reference.reduction = "rnaPCA")

predictions_49_2 <- TransferData(anchorset = anchors_49_2, refdata = multiome_object$res0.8_50nn_40PC_45LSI)
previous_object <- AddMetaData(previous_object, metadata = predictions_49_2)

Idents(previous_object) <- 'predicted.id'
DimPlot(previous_object, label = T)

predictions_49_2$query_id <- previous_object$clusters49

old_obj_predictions_matrix_49 <- predictions_49_2%>%
  group_by(query_id)%>%
  summarise(across(starts_with("prediction.score"), mean, na.rm = TRUE))

colnames(old_obj_predictions_matrix_49) <- str_remove(colnames(old_obj_predictions_matrix_49), 'prediction.score.')

colnames(old_obj_predictions_matrix_49) <- as.numeric(colnames(old_obj_predictions_matrix_49))

old_obj_predictions_matrix_49%>%
  order()

old_obj_predictions_matrix_49[,2:(ncol(old_obj_predictions_matrix_49)-1)]%>%
  as.matrix()%>%
  scale()%>%
  Heatmap(
    cluster_rows = FALSE, 
    cluster_columns = FALSE,
    row_order = 1:49,
    row_labels = 0:48,
    column_order = order(as.numeric(colnames(old_obj_predictions_matrix_49[,2:(ncol(old_obj_predictions_matrix_49)-1)]))),
    row_title = 'previous',
    column_title = 'multiome')


DimPlot(previous_object, group.by = 'clusters49', label =T)
DimPlot(multiome_object, group.by = 'res0.8_50nn_40PC_45LSI', label = T, reduction = 'harmony_wnn.umap')
