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

renamed_clusters <- c("0" = "9",
                      "1" = "1", 
                      "2" = "2",
                      "3" = "13",
                      "4" = "14",
                      "5" = "15",
                      "6" = "3",
                      "7" = "10",
                      "8" = "17",
                      "9" = "4",
                      "10" = "5",
                      "11" = "16",
                      '12' ='6',
                      '13' = NA,
                      '14'='11',
                      '15'='7',
                      '16'='8',
                      '17'= '18',
                      '18'= '12',
                      '19' = '19'
                      )

previous_object@meta.data$new_cluster_ids <- NULL
previous_object@meta.data$new_cluster_ids <- renamed_clusters[as.character(previous_object@meta.data$clusters20)]


Idents(previous_object) = 'new_cluster_ids'
DimPlot(previous_object, label = T)

#### Find anchorsets and transfer labels #####

#### new_cluster_ids ####
Idents(previous_object) <- 'new_cluster_ids'
anchors <- FindTransferAnchors(reference = previous_object, query = multiome_object,
                               reference.reduction = "pca")

predictions <- TransferData(anchorset = anchors, refdata = previous_object$new_cluster_ids)
multiome_object <- AddMetaData(multiome_object, metadata = predictions)

Idents(multiome_object) <- 'predicted.id'

arr <- list(x = -15, y = -15, x_len = 4, y_len = 4)

label_transferred_object = DimPlot(multiome_object, label = T, reduction = 'harmony_wnn.umap')+
  theme_void()+
  theme(legend.position = 'none')+
  annotate('text',
           x = -17, y = -12.5, label = 'UMAP_2', angle = 90, size =3)+
  annotate('text',
           x = -12.5, y = -17, label = 'UMAP_1', size =3)+
  annotate("segment", 
           x = arr$x, xend = arr$x + c(arr$x_len, 0), 
           y = arr$y, yend = arr$y + c(0, arr$y_len), 
           arrow = arrow(type = "closed", length = unit(10, 'pt'))) 
  

label_transferred_object

ggsave(plot = label_transferred_object,
       file = "label_transfer.tiff",
       device = "tiff",
       units = "in",
       width = 2.5,
       height = 2.5,
       path = "Bachelors Thesis/Plots/Figure 2",
       dpi = 300)


DimPlot(previous_object, label =T)


