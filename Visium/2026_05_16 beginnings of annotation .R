library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(dittoSeq) 
library(cowplot)

parker_2024 = readRDS("~/Desktop/parker_object.rds")
Idents(parker_2024) = 'clusters49'


obj  = readRDS("~/Desktop/optimal_clustering_rna_only.rds")
Idents(obj) = 'final_clusters'
DimPlot(obj)

vis = readRDS('/Users/ggraham/Desktop/Visium/vis.rds')
Idents(vis) = 'seurat_clusters.projected0.4'
DimPlot(vis)

x11(width = 14)
SpatialFeaturePlot(vis, features = "nCount_Spatial.Polygons", ncol = 4, pt.size.factor = 0.01) 

x11(width = 14)
SpatialDimPlot(vis, images = c('s_6P17.polygons'))
# Okabe-Ito + extended colorblind-safe palette
cb_colors <- c(
  "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2",
  "#D55E00", "#CC79A7", "#000000", "#882255", "#117733",
  "#44AA99", "#88CCEE", "#DDCC77", "#332288", "#AA4499",
  "#999933", "#661100", "#6699CC", "#888888", "#E8601C",
  "#F6C141", "#1965B0", "#7BAFDE", "#4EB265", "#CAE0AB",
  "#F7F056", "#DC050C", "#72190E", "#8DD593", "#C6DEC7"
)

DotPlot(vis, 
        c('gja1b',
          's100b',
          'gad2',
          'slc17a6b',
          'slc17a7a',
          'p2ry12',
          'ptprc',
          'mpz',
          'th',
          'th2',
          'meis2a',
          'hmx2',
          'hmx3a',
          'npy',
          'sst1.1'
          ),
        assay = 'sketch')+
  coord_flip()

markers_6 = FindMarkers(vis, 6)
markers_20 = FindMarkers(vis, 20)

markers_24 = FindMarkers(vis, 24)

annots <- c(
  '0'  = 'Unannotated_0',   # placeholder — fill in what you know
  '1'  = 'Ependymoglia',
  '2'  = 'Diencephalon', # this one is problematic because it is very heavily in OT so needs to be split I believe
  '3'  = 'Vd_r_Vs_m',
  '4'  = 'Diencephalon',
  '5'  = 'Diencephalon',
  '6'  = 'Dm',
  '7'  = 'Dl_blob_v2',
  '8'  = 'aPOA',
  '9'  = 'Blood',
  '10' = 'Ependymoglia',
  '11' = 'Dd_Dp',
  '12' = 'Diencephalon',
  '13' = 'caudal_subpallium',
  '14' = 'Blood',
  '15' = 'Immune',
  '16' = 'Diencephalon',
  '17' = 'Diencephalon',
  '18' = 'Diencephalon',
  '19' = 'Dl_g',
  '20' = 'Diencephalon', 
  '21' = 'Vl',
  '22' = 'pPOA',
  '23' = 'Ependymoglia',
  '24' = 'Interneuron',
  '25' = 'Ependymoglia',
  '26' = 'Ependymoglia',
  '27' = 'Interneuron',
  '28' = 'Diencephalon',
  '29' = 'Diencephalon',
  '30' = 'Diencephalon'
)

# Safe mapping — NAs become 'Unannotated' instead of silently breaking
vis$annots_0.4 <- ifelse(
  vis$seurat_clusters.projected0.4 %in% names(annots),
  annots[as.character(vis$seurat_clusters.projected0.4)],
  paste0("Unannotated_", vis$seurat_clusters.projected0.4)
)


# Get number of unique annotations and assign colors
annot_levels <- sort(unique(vis$annots_0.4))
n_annots     <- length(annot_levels)
color_map    <- setNames(cb_colors[seq_len(n_annots)], annot_levels)

x11(width = 16)
SpatialDimPlot(vis,
               images         = 's_4P10.polygons',
               group.by       = 'annots_0.4',
               cols           = color_map,
               pt.size.factor = 0.5) +
  guides(fill = guide_legend(override.aes = list(size = 6)))

DotPlot(vis, group.by = 'annots_0.4', 
        features = 'mpz')

SpatialDimPlot(vis, images = c('s_6P17.polygons'), group.by = 'annots_0.4',
                              pt.size.factor = 10)
x11(width = 14)
SpatialDimPlot(vis%>%subset(annots_0.4=='Dorsal_Subpallium'),
               images = c('s_6P17.polygons'),
               group.by = 'annots_0.4', 
               pt.size.factor = 1)
x11(width = 14)
SpatialDimPlot(vis%>%subset(seurat_clusters.projected0.4=='17'),
               group.by = 'seurat_clusters.projected0.4', 
               pt.size.factor = 0.01)




