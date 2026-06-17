library(Seurat)
library(zellkonverter)   # BiocManager::install("zellkonverter")
library(SingleCellExperiment)

# ── Multiome → h5ad ────────────────────────────────────────────────────────
multiome <- readRDS("/Users/ggraham/Desktop/optimal_clustering_rna_only.rds")
DefaultAssay(multiome) <- "RNA"

# Seurat v5: join layers first (RNA may be split into counts.1, counts.2, etc.)
multiome[["RNA"]] <- JoinLayers(multiome[["RNA"]])

# Convert to SCE, then write h5ad
sce_ref <- as.SingleCellExperiment(multiome, assay = "RNA")
writeH5AD(sce_ref, file = "/Users/ggraham/Desktop/multiome_ref.h5ad")


# ── Visium HD → h5ad ───────────────────────────────────────────────────────
# here, also add dissection metadata column


vis <- readRDS("/Users/ggraham/Desktop/Visium/vis.rds")
DefaultAssay(vis) <- "sketch"   # try "RNA" if this errors
vis <- FindVariableFeatures(vis)
vis <- ScaleData(vis, verbose = FALSE)
vis <- RunPCA(vis, reduction.name = "pca", verbose = FALSE)

ndims <- 15

vis <- FindNeighbors(vis, reduction = "pca", 
                             dims = 1:ndims, graph.name = "snn")
vis <- FindClusters(vis, graph.name = "snn",
                            resolution = 0.4, cluster.name = "clusters",
                            verbose = FALSE)
vis <- RunUMAP(vis, reduction = "pca", 
                       dims = 1:ndims, reduction.name = "umap")

vis <- ProjectData(
  object        = vis,
  assay         = "Spatial.Polygons",
  full.reduction   = "full.pca",
  sketched.assay   = "sketch",
  sketched.reduction = "pca",
  dims          = 1:ndims,
  refdata       = list(clusters.projected = "clusters")
)

sections =c(
  '6P17',
  '3P5',
  '4P5',
  '4P10'
)

tog = data.frame()
for(i in sections){
  data = read.csv(
    paste0('/Users/ggraham/Desktop/multiome_poa/Visium/Groupings/',i,' Coltan Dissection.csv')
    )
  tog = rbind(tog, data)
  
}

included = tog$Coltan.DIssection
names(included) = tog$Barcode
seurat_barcodes <- gsub("_[0-9]+$", "", Cells(vis))

vis$dissection <- unname(included[seurat_barcodes])

#saveRDS(vis, "/Users/ggraham/Desktop/Visium/vis_meta.rds")
vis_dissection =subset(vis, !is.na(vis$dissection))

# Drop reductions that were computed on the full object (mismatched after subset)
vis_dissection@reductions <- list()

# Also drop graphs for the same reason
vis_dissection@graphs <- list()

vis_dissection[["Spatial.Polygons"]] <- JoinLayers(vis_dissection[["Spatial.Polygons"]])

sce_vis <- as.SingleCellExperiment(vis_dissection, assay = "Spatial.Polygons")
sce_vis$ID = paste0(                              
  sce_vis$Slice,
   sce_vis$Sex,
  sce_vis$individual)

writeH5AD(sce_vis, file = "/Users/ggraham/Desktop/vis.h5ad")

saveRDS(vis_dissection, "/Users/ggraham/Desktop/Visium/vis_dissection.rds")



