library(Seurat)
library(patchwork)
library(ggplot2)
library(tidyverse)

vis2 = readRDS("C:/Users/Gabe/Desktop/Visium/vis_better_barcodes_dissection_c2l95_projection.rds")
vis = vis2
### here, I want to remove the final_ from each cell barcode and append with
### the Slice that it comes from, then, I want to do the same to each grouping
images = unique(vis$Slice)

barcodes_included =data.frame()
for(i in images){
  print(i)
  data = read.csv(paste0("C:/Users/Gabe/Desktop/multiome_poa/Visium/Groupings/Maruska Regions Data Informed/2026_05_17 ",tolower(i)," Anatomical.csv"))%>%
    subset(!is.na(Anatomical))
  data$Barcode_updated =paste0(data$Barcode, i)
  barcodes_included = rbind(data, barcodes_included)
}
vis$Anatomical <- "No Region"

idx <- match(colnames(vis), barcodes_included$Barcode_updated)

vis$Anatomical[!is.na(idx)] <- barcodes_included$Anatomical[idx[!is.na(idx)]]

renamed = c(
  'could be nPPa'= 'What',
  'could be nMMp' = 'nMMp',
  'could be nPMP' = 'nPMp',
  'dorsolateral nucleus' = 'Dln-r',
  'dorsolateral nucleus caudal part' = 'Dln-c',
  'Dl-v2 dorsal part' = 'Dl-v2d',
  'TGn' = 'TGN'
)

renaming_col = ifelse(vis$Anatomical %in% c(names(renamed)), renamed[vis$Anatomical], 
                      vis$Anatomical)  
vis$Anatomical = renaming_col

vis$Anatomical = ifelse(vis$Anatomical%in%c(
                                            'not sure 21',
                                            'OT'),
                        NA,
                        vis$Anatomical)

vis$Predicted_final_clusters <- sub(".*_", "", vis$dominant_c2l_cluster)

#SpatialDimPlot(vis, 'Anatomical', images = "s_6P17.polygons")

heatmat = table(vis$Predicted_final_clusters,
                vis$Anatomical)%>%
  as.data.frame.matrix()

library(pheatmap)
pheatmap((heatmat/rowSums(heatmat))%>%scale(),
         cluster_rows = T,
         cluster_cols = T)


saveRDS(vis, "C:/Users/Gabe/Desktop/Visium/vis_better_barcodes_dissection_c2l05_projection_anatomical.rds")

#write.csv(barcodes_included, "Visium/Groupings/Coltan DIssections together and barcodes.csv")

