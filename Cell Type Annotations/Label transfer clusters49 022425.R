#Cell type annotations 121924
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
  #library(CytoTRACE)
  # SCTRransform_mean_plot <- readRDS("R/Gabe/SCTRransform_mean_plot.rds")
  #mac.neg.bin <- readRDS(file = 'R/Gabe/mac.neg.bin.rds')
  library('glmGamPoi')
  library(scran)
  library(parallel)
  library(factoextra)
  library(readxl)
  library(factoextra)
  library(forcats)
  library(ggrepel)
  library(biomaRt)
  #mean_cell <- readRDS('R/Gabe/mean_cell.rds')
  library(openxlsx)
  #clown_go <- readRDS('R/Gabe/clown_go.rds')
  library(ComplexHeatmap)
  
}
multiome_object <- readRDS('C:/Users/Gabe/Desktop/RNA Object.rds')
Idents(multiome_object) <- "harmony.wnn_res0.4_clusters"


#### Find anchorsets and transfer labels #####
previous_object <- readRDS('C:/Users/Gabe/Desktop/Old object.rds')

#### Clusters 49 ####
Idents(previous_object) <- 'clusters49'
anchors <- FindTransferAnchors(reference = previous_object, query = multiome_object,
                               reference.reduction = "pca")

predictions <- TransferData(anchorset = anchors, refdata = previous_object$clusters49)
multiome_object <- AddMetaData(multiome_object, metadata = predictions)

Idents(multiome_object) <- 'predicted.id'
DimPlot(multiome_object, label = T)

##other way around
Idents(multiome_object) <- 'harmony.wnn_res0.4_clusters'
anchors2 <- FindTransferAnchors(reference = multiome_object, query = previous_object,
                               reference.reduction = "pca")

predictions2 <- TransferData(anchorset = anchors2, refdata = multiome_object$harmony.wnn_res0.4_clusters)
previous_object <- AddMetaData(previous_object, metadata = predictions2)

Idents(previous_object) <- 'predicted.id'
DimPlot(previous_object, label = T)


old_obj_predictions_matrix <- predictions2%>%
  group_by(query_id)%>%
  summarise(across(starts_with("prediction.score"), mean, na.rm = TRUE))

colnames(old_obj_predictions_matrix) <- str_remove(colnames(old_obj_predictions_matrix), 'prediction.score.')

old_obj_predictions_matrix[,2:33]%>%
  as.matrix()%>%
  scale()%>%
  Heatmap(row_labels = old_obj_predictions_matrix$query_id,
          cluster_rows = FALSE, 
          cluster_columns = FALSE,
          row_title = 'previous',
          column_title = 'multiome')

DimPlot(previous_object, label = T)
##cluster 27 does not align strongly to any clusters, maybe try with the relabeled object

