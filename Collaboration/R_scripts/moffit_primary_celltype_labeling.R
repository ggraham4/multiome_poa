
library(Seurat)
library(Signac)
library(ggplot2)
library(scCustomize)
library(patchwork)

nemo <- readRDS('/Users/kathrynleatherbury/GaTech Dropbox/CoS/BioSci/BioSci-Streelman/MultiomeProjects/Katie/Clownfish/Seurat_Data/seurat_objects/nemo_gex.sct.atac_integration_10.27.24.rds')

DefaultAssay(nemo) <- 'RNA'
Idents(nemo) <- 'orig.ident'

gc()
gc()



# increase clustering resolution to 0.6 to separate glut and gaba clusters more accurately 


# build a joint neighbor graph using both gex + atac assays - WNN
nemo <- FindMultiModalNeighbors(
  object = nemo,
  reduction.list = list("harmony.rna", 'atacLSI'),
  dims.list = list(1:50, 2:50),
  verbose = TRUE,
  prune.SNN = 0,
  knn.graph.name = "harmony.wknn",
  snn.graph.name = "harmony.wsnn",
  weighted.nn.name = "harmony.wnn"
)

# build a joint UMAP visualization
nemo <- RunUMAP(
  object = nemo,
  nn.name = "harmony.wnn",
  reduction.name = "harmony_wnn.umap",
  reduction.key = "harmony_wnnUMAP_",
  verbose = TRUE,
  min.dist = 0.4,
  metric = "euclidean"
)

nemo <- FindClusters(nemo,
                     graph.name = "harmony.wsnn",
                     algorithm = 3,
                     verbose = T,
                     resolution = 0.6,
                     cluster.name = 'harmony.wnn_res0.6_clusters')



####################################################



##### Primary Cell Type Labeling based on Moffit et al. (2018). Molecular, spatial, and functional single-cell profiling of the hypothalamic preoptic region. #####


nemo$primary_celltypes_Moffit <- nemo$harmony.wnn_res0.6_clusters
nemo$primary_celltypes_Moffit <- as.character(nemo$primary_celltypes_Moffit)
Idents(nemo) <- 'harmony.wnn_res0.6_clusters'


p1 <- DimPlot(object = nemo, reduction = "harmony_wnn.umap", pt.size = .3, group.by = "harmony.wnn_res0.6_clusters", shuffle = T, label = F) + NoLegend() + NoAxes() +
  labs(title = 'nemo Multiome', subtitle = 'POA Clusters \n') + 
  theme(plot.title = element_text(hjust=0.5, face="bold", size = 15, vjust = 0), plot.subtitle = 
          element_text(hjust = 0.5, 
                       size = 12, 
                       face="italic"))
p1 <- LabelClusters(p1, id = "harmony.wnn_res0.6_clusters",  fontface = "bold", color = "black")
p1

colors <- c("gray87", "lightblue", "#41B6C4", 'deepskyblue3', "#225EA8", "#253494")



# Newly Formed Oligodendrocytes (NFOs) = 33
FeaturePlot_scCustom(seurat_object = nemo,
                     features = 'gpr17',
                     order = TRUE,
                     label = T,
                     repel = TRUE,
                     reduction = 'harmony_wnn.umap',
                     colors_use = colors,
                     na_cutoff = 1, pt.size = 0.3)
nemo$primary_celltypes_Moffit[nemo$primary_celltypes_Moffit == 33] = 'NFO'

p2 <- DimPlot(nemo, reduction = "harmony_wnn.umap", label = TRUE, pt.size = 0.5, group.by = 'primary_celltypes_Moffit') + labs(title = 'Moffit Primary Cell Type Predictions') + NoAxes()
p2

p3 <- FeaturePlot_scCustom(seurat_object = nemo,
                           features = 'gpr17',
                           order = TRUE,
                           label = F,
                           repel = TRUE,
                           reduction = 'harmony_wnn.umap',
                           colors_use = colors,
                           na_cutoff = 1, 
                           pt.size = 0.3)+ NoAxes()
p2 + p3

# mature oligo = 4, 34
FeaturePlot_scCustom(seurat_object = nemo,
                           features = 'myrf',
                           order = TRUE,
                           label = T,
                           repel = TRUE,
                           reduction = 'harmony_wnn.umap',
                           colors_use = colors,
                           na_cutoff = 1, pt.size = 0.3)
nemo$primary_celltypes_Moffit[nemo$primary_celltypes_Moffit == 4] = 'MO'
nemo$primary_celltypes_Moffit[nemo$primary_celltypes_Moffit == 34] = 'MO'

p2 <- DimPlot(nemo, reduction = "harmony_wnn.umap", label = TRUE, pt.size = 0.5, group.by = 'primary_celltypes_Moffit') + labs(title = 'Moffit Primary Cell Type Predictions') + NoAxes()
p2

p3 <- FeaturePlot_scCustom(seurat_object = nemo,
                           features = 'myrf',
                           order = TRUE,
                           label = F,
                           repel = TRUE,
                           reduction = 'harmony_wnn.umap',
                           colors_use = colors,
                           na_cutoff = 1, 
                           pt.size = 0.3)+ NoAxes()
p2 + p3



# Radial Glia (astrocytes in moffit) = 2
p1 + FeaturePlot_scCustom(seurat_object = nemo,
                          features = 's1pr1',
                          order = TRUE,
                          label = T,
                          repel = TRUE,
                          reduction = 'harmony_wnn.umap',
                          colors_use = colors,
                          na_cutoff = 1, pt.size = 0.3)
nemo$primary_celltypes_Moffit[nemo$primary_celltypes_Moffit == 2] = 'Astrocytes'

p2 <- DimPlot(nemo, reduction = "harmony_wnn.umap", label = TRUE, pt.size = 0.5, group.by = 'primary_celltypes_Moffit') + labs(title = 'Moffit Primary Cell Type Predictions') + NoAxes()


p3 <- FeaturePlot_scCustom(seurat_object = nemo,
                           features = 's1pr1',
                           order = TRUE,
                           label = F,
                           repel = TRUE,
                           reduction = 'harmony_wnn.umap',
                           colors_use = colors,
                           na_cutoff = 1, 
                           pt.size = 0.3)+ NoAxes()
p2 + p3



# microglia = 15
FeaturePlot_scCustom(seurat_object = nemo,
                     features = 'csf1rb',
                     order = TRUE,
                     label = T,
                     repel = TRUE,
                     reduction = 'harmony_wnn.umap',
                     colors_use = colors,
                     na_cutoff = 1, pt.size = 0.3)
nemo$primary_celltypes_Moffit[nemo$primary_celltypes_Moffit == 15] = 'Microglia'

p2 <- DimPlot(nemo, reduction = "harmony_wnn.umap", label = TRUE, pt.size = 0.5, group.by = 'primary_celltypes_Moffit') + labs(title = 'Moffit Primary Cell Type Predictions') + NoAxes()


p3 <- FeaturePlot_scCustom(seurat_object = nemo,
                           features = 'csf1rb',
                           order = TRUE,
                           label = F,
                           repel = TRUE,
                           reduction = 'harmony_wnn.umap',
                           colors_use = colors,
                           na_cutoff = 1, 
                           pt.size = 0.3)+ NoAxes()
p2 + p3



# OPCs = 19
FeaturePlot_scCustom(seurat_object = nemo,
                     features = 'aplnra',
                     order = TRUE,
                     label = T,
                     repel = TRUE,
                     reduction = 'harmony_wnn.umap',
                     colors_use = colors,
                     na_cutoff = 1, pt.size = 0.3)
nemo$primary_celltypes_Moffit[nemo$primary_celltypes_Moffit == 19] = 'OPC'

p2 <- DimPlot(nemo, reduction = "harmony_wnn.umap", label = TRUE, pt.size = 0.5, group.by = 'primary_celltypes_Moffit') + labs(title = 'Moffit Primary Cell Type Predictions') + NoAxes()


p3 <- FeaturePlot_scCustom(seurat_object = nemo,
                           features = 'aplnra',
                           order = TRUE,
                           label = F,
                           repel = TRUE,
                           reduction = 'harmony_wnn.umap',
                           colors_use = colors,
                           na_cutoff = 1, 
                           pt.size = 0.3)+ NoAxes()
p2 + p3



# Macrophages (Leukocytes) = 27
FeaturePlot_scCustom(seurat_object = nemo,
                     features = 'tbx21',
                     order = TRUE,
                     label = T,
                     repel = TRUE,
                     reduction = 'harmony_wnn.umap',
                     colors_use = colors,
                     na_cutoff = 1, pt.size = 0.3)
nemo$primary_celltypes_Moffit[nemo$primary_celltypes_Moffit == 27] = 'Macrophages'

p2 <- DimPlot(nemo, reduction = "harmony_wnn.umap", label = TRUE, pt.size = 0.5, group.by = 'primary_celltypes_Moffit') + labs(title = 'Moffit Primary Cell Type Predictions') + NoAxes()


p3 <- FeaturePlot_scCustom(seurat_object = nemo,
                           features = 'tbx21',
                           order = TRUE,
                           label = F,
                           repel = TRUE,
                           reduction = 'harmony_wnn.umap',
                           colors_use = colors,
                           na_cutoff = 1, 
                           pt.size = 0.3)+ NoAxes()
p2 + p3


# pericytes = 30
FeaturePlot_scCustom(seurat_object = nemo,
                     features = 'osr1',
                     order = TRUE,
                     label = T,
                     repel = TRUE,
                     reduction = 'harmony_wnn.umap',
                     colors_use = colors,
                     na_cutoff = 1, pt.size = 0.3)
nemo$primary_celltypes_Moffit[nemo$primary_celltypes_Moffit == 30] = 'Pericytes'

p2 <- DimPlot(nemo, reduction = "harmony_wnn.umap", label = TRUE, pt.size = 0.5, group.by = 'primary_celltypes_Moffit') + labs(title = 'Moffit Primary Cell Type Predictions') + NoAxes()


p3 <- FeaturePlot_scCustom(seurat_object = nemo,
                           features = 'osr1',
                           order = TRUE,
                           label = F,
                           repel = TRUE,
                           reduction = 'harmony_wnn.umap',
                           colors_use = colors,
                           na_cutoff = 1, 
                           pt.size = 0.3)+ NoAxes()
p2 + p3

# Fibroblast marker = 25
FeaturePlot_scCustom(seurat_object = nemo,
                     features = 'dcn',
                     order = TRUE,
                     label = T,
                     repel = TRUE,
                     reduction = 'harmony_wnn.umap',
                     colors_use = colors,
                     na_cutoff = 1, pt.size = 0.3)
nemo$primary_celltypes_Moffit[nemo$primary_celltypes_Moffit == 25] = 'Fibroblasts'

p2 <- DimPlot(nemo, reduction = "harmony_wnn.umap", label = TRUE, pt.size = 0.5, group.by = 'primary_celltypes_Moffit') + labs(title = 'Moffit Primary Cell Type Predictions') + NoAxes()


p3 <- FeaturePlot_scCustom(seurat_object = nemo,
                           features = 'dcn',
                           order = TRUE,
                           label = F,
                           repel = TRUE,
                           reduction = 'harmony_wnn.umap',
                           colors_use = colors,
                           na_cutoff = 1, 
                           pt.size = 0.3)+ NoAxes()
p2 + p3

# Ependymal cells = 23
FeaturePlot_scCustom(seurat_object = nemo,
                     features = 'fsip1',
                     order = TRUE,
                     label = T,
                     repel = TRUE,
                     reduction = 'harmony_wnn.umap',
                     colors_use = colors,
                     na_cutoff = 1, pt.size = 0.3)
nemo$primary_celltypes_Moffit[nemo$primary_celltypes_Moffit == 23] = 'Ependymal Cells'

p2 <- DimPlot(nemo, reduction = "harmony_wnn.umap", label = TRUE, pt.size = 0.5, group.by = 'primary_celltypes_Moffit') + labs(title = 'Moffit Primary Cell Type Predictions') + NoAxes()


p3 <- FeaturePlot_scCustom(seurat_object = nemo,
                           features = 'fsip1',
                           order = TRUE,
                           label = F,
                           repel = TRUE,
                           reduction = 'harmony_wnn.umap',
                           colors_use = colors,
                           na_cutoff = 1, 
                           pt.size = 0.3)+ NoAxes()
p2 + p3

# Endothelial cells = 7
FeaturePlot_scCustom(seurat_object = nemo,
                     features = 'ephb4a',
                     order = TRUE,
                     label = T,
                     repel = TRUE,
                     reduction = 'harmony_wnn.umap',
                     colors_use = colors,
                     na_cutoff = 1, pt.size = 0.3)
nemo$primary_celltypes_Moffit[nemo$primary_celltypes_Moffit == 7] = 'Endothelial Cells'

p2 <- DimPlot(nemo, reduction = "harmony_wnn.umap", label = TRUE, pt.size = 0.5, group.by = 'primary_celltypes_Moffit') + labs(title = 'Moffit Primary Cell Type Predictions') + NoAxes()


p3 <- FeaturePlot_scCustom(seurat_object = nemo,
                           features = 'ephb4a',
                           order = TRUE,
                           label = F,
                           repel = TRUE,
                           reduction = 'harmony_wnn.umap',
                           colors_use = colors,
                           na_cutoff = 1, 
                           pt.size = 0.3)+ NoAxes()
p2 + p3

DimPlot(nemo, reduction = "harmony_wnn.umap", label = TRUE, pt.size = 0.5, group.by = 'primary_celltypes_Moffit') + labs(title = 'Moffit Primary Cell Type Predictions') + NoAxes()


# GLUT markers
p2 <- FeaturePlot_scCustom(seurat_object = nemo,
                           features = 'slc17a6b',
                           order = TRUE,
                           label = T,
                           repel = TRUE,
                           reduction = 'harmony_wnn.umap',
                           colors_use = colors,
                           na_cutoff = 1, pt.size = 0.3)+ NoAxes()
p1 + p2
 nemo$primary_celltypes_Moffit[nemo$primary_celltypes_Moffit == 10] = 'Glutamatergic'
nemo$primary_celltypes_Moffit[nemo$primary_celltypes_Moffit == 11] = 'Glutamatergic'
nemo$primary_celltypes_Moffit[nemo$primary_celltypes_Moffit == 7] = 'Glutamatergic'
nemo$primary_celltypes_Moffit[nemo$primary_celltypes_Moffit == 12] = 'Glutamatergic'
nemo$primary_celltypes_Moffit[nemo$primary_celltypes_Moffit == 18] = 'Glutamatergic'
nemo$primary_celltypes_Moffit[nemo$primary_celltypes_Moffit == 1] = 'Glutamatergic'
nemo$primary_celltypes_Moffit[nemo$primary_celltypes_Moffit == 24] = 'Glutamatergic'
nemo$primary_celltypes_Moffit[nemo$primary_celltypes_Moffit == 17] = 'Glutamatergic'
nemo$primary_celltypes_Moffit[nemo$primary_celltypes_Moffit == 20] = 'Glutamatergic'
nemo$primary_celltypes_Moffit[nemo$primary_celltypes_Moffit == 22] = 'Glutamatergic'
nemo$primary_celltypes_Moffit[nemo$primary_celltypes_Moffit == 28] = 'Glutamatergic'
nemo$primary_celltypes_Moffit[nemo$primary_celltypes_Moffit == 13] = 'Glutamatergic'
nemo$primary_celltypes_Moffit[nemo$primary_celltypes_Moffit == 16] = 'Glutamatergic'
nemo$primary_celltypes_Moffit[nemo$primary_celltypes_Moffit == 14] = 'Glutamatergic'
nemo$primary_celltypes_Moffit[nemo$primary_celltypes_Moffit == 6] = 'Glutamatergic'
nemo$primary_celltypes_Moffit[nemo$primary_celltypes_Moffit == 9] = 'Glutamatergic'
nemo$primary_celltypes_Moffit[nemo$primary_celltypes_Moffit == 3] = 'Glutamatergic'

p2 <- DimPlot(nemo, reduction = "harmony_wnn.umap", label = TRUE, pt.size = 0.5, group.by = 'primary_celltypes_Moffit') + labs(title = 'Moffit Primary Cell Type Predictions') + NoAxes()


p3 <- FeaturePlot_scCustom(seurat_object = nemo,
                           features = 'slc17a6b',
                           order = TRUE,
                           label = F,
                           repel = TRUE,
                           reduction = 'harmony_wnn.umap',
                           colors_use = colors,
                           na_cutoff = 1, 
                           pt.size = 0.3)+ NoAxes()
p2 + p3




# GABA markers 
p2 <- FeaturePlot_scCustom(seurat_object = nemo,
                           features = 'gad2',
                           order = TRUE,
                           label = T,
                           repel = TRUE,
                           reduction = 'harmony_wnn.umap',
                           colors_use = colors,
                           na_cutoff = 1, pt.size = 0.3) + NoAxes()
p1 + p2
nemo$primary_celltypes_Moffit[nemo$primary_celltypes_Moffit == 0] = 'GABAergic'
nemo$primary_celltypes_Moffit[nemo$primary_celltypes_Moffit == 5] = 'GABAergic'
nemo$primary_celltypes_Moffit[nemo$primary_celltypes_Moffit == 8] = 'GABAergic'
nemo$primary_celltypes_Moffit[nemo$primary_celltypes_Moffit == 21] = 'GABAergic'
nemo$primary_celltypes_Moffit[nemo$primary_celltypes_Moffit == 26] = 'GABAergic'
nemo$primary_celltypes_Moffit[nemo$primary_celltypes_Moffit == 31] = 'GABAergic'
nemo$primary_celltypes_Moffit[nemo$primary_celltypes_Moffit == 29] = 'GABAergic'
nemo$primary_celltypes_Moffit[nemo$primary_celltypes_Moffit == 32] = 'GABAergic'

p2 <- DimPlot(nemo, reduction = "harmony_wnn.umap", label = TRUE, pt.size = 0.5, group.by = 'primary_celltypes_Moffit') + labs(title = 'Moffit Primary Cell Type Predictions') + NoAxes()


p3 <- FeaturePlot_scCustom(seurat_object = nemo,
                           features = 'gad2',
                           order = TRUE,
                           label = F,
                           repel = TRUE,
                           reduction = 'harmony_wnn.umap',
                           colors_use = colors,
                           na_cutoff = 1, 
                           pt.size = 0.3)+ NoAxes()
p2 + p3

DimPlot(nemo, reduction = "harmony_wnn.umap", label = TRUE, pt.size = 0.5, group.by = 'primary_celltypes_Moffit') + labs(title = 'Moffit Primary Cell Type Predictions') + NoAxes()


Idents(nemo) <- "primary_celltypes_Moffit"
new.cluster.ids <- c("Endothelial","Astrocytes", "Glutamatergic", "GABAergic", "Ependymal", "OPC",  "Microglia", "MO", "Macrophages", "NFO", "Pericytes", "Fibroblasts")
nemo <- RenameIdents(nemo, new.cluster.ids)
nemo$primary_celltypes_Moffit <- Idents(nemo)

a <- DimPlot_scCustom(seurat_object = nemo, 
                      reduction = "harmony_wnn.umap", 
                      label = F, 
                      pt.size = 0.5, repel = T,
                      group.by = 'primary_celltypes_Moffit',
                      colors_use = c('gold3','orange2', 'deepskyblue','lightcoral', 'chartreuse3','lightslateblue', 'turquoise3','green3', 'magenta1', 'dodgerblue2', 'seagreen2', 'hotpink2')) + 
  labs(title = 'Parent Cell Type IDs Predicted by Mofit et. al (2018)') + NoAxes() + NoLegend()
a <- LabelClusters(a, id = "primary_celltypes_Moffit",  fontface = "bold", color = "black")
a