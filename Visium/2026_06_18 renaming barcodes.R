library(Seurat)

vis2 = readRDS("~/Desktop/Visium/vis.rds")
vis = vis2
### here, I want to remove the final_ from each cell barcode and append with
### the Slice that it comes from, then, I want to do the same to each grouping
images = unique(vis$Slice)

barcodes_included =data.frame()
for(i in images){
  print(i)
data = read.csv(paste0("Visium/Groupings/",i," Coltan DIssection.csv"))%>%
  subset(!is.na( Coltan.DIssection))
data$Barcode_updated =paste0(data$Barcode, i)
barcodes_included = rbind(data, barcodes_included)
}

print(head(data$Barcode))
print(head(data$Barcode_updated))

#Second
barcodes = colnames(vis@assays$Spatial.Polygons)
print(barcodes%>%head())
barcodes_no_ <- sub("_[0-9]+$", "", barcodes)
print(barcodes_no_%>%head())
barcodes_image = paste0(barcodes_no_, vis$Slice)
vis <- RenameCells(vis, new.names = barcodes_image)

vis$dissection = ifelse(colnames(vis)%in%barcodes_included$Barcode_updated,
                        'In',
                        'Out')
table(vis$dissection)

saveRDS(vis, "~/Desktop/Visium/vis_better_barcodes.rds")
write.csv(barcodes_included, "Visium/Groupings/Coltan DIssections together and barcodes.csv")

SpatialDimPlot(vis, 'dissection', images = "s_6P17.polygons")
