library(Seurat)
integrated = readRDS( "~/Desktop/2025_11_24_poa_tel_integrated.rds")

# vim LOC111579533

FeaturePlot(integrated, 'LOC111579533', reduction = 'harmony.clusters')
FeaturePlot(integrated, 'gfap', reduction = 'harmony.clusters')
DimPlot(integrated, reduction = 'harmony.clusters', label = T)

integrated = JoinLayers(integrated)
marks_17 = FindMarkers(integrated, 17) # ependymal cells

FeaturePlot(integrated, 'nFeature_RNA', reduction = 'harmony.clusters')

marks_5 = FindMarkers(integrated, 5, only.pos = T) 

library(hdWGCNA)

FeaturePlot(integrated, 'npy', reduction = 'harmony.clusters')
FeaturePlot(integrated, 'oxt', reduction = 'harmony.clusters')
FeaturePlot(integrated, 'gal', reduction = 'harmony.clusters')
FeaturePlot(integrated, 'oxt', reduction = 'harmony.clusters')

marks_3 = FindMarkers(integrated, 3, only.pos = T) 
marks_16 = FindMarkers(integrated, 16, only.pos = T) 

DotPlot(integrated, c('gal',
                      'npy',
                      'cckb',
                      'tac3a',
                      'kiss1',
                      'esr2b',
                      'ar',
                      'pgr',
                      'sst1.1',
                      'hmx3a',
                      'hmx2',
                      'galr1a'))+
  coord_flip()


# ok its cluster 2

DimPlot(integrated, cells.highlight = WhichCells(integrated, idents = c("17")), reduction = 'harmony.clusters')

integrated = FindSubCluster(integrated, '17', graph.name = 'RNA_snn')

DimPlot(integrated, group.by = 'sub.cluster', reduction = 'harmony.clusters', label = T)

marks_17_0 = FindMarkers(integrated, '17_0', group.by ='sub.cluster', only.pos =T)



### reclustering 
library(Seurat)
library(harmony)

# Make sure mito percent is calculated if not already
integrated[["percent.mt"]] <- PercentageFeatureSet(integrated, pattern = "^mt-")  # use "^MT-" for human

# Re-run standard preprocessing
integrated <- NormalizeData(integrated)
integrated <- FindVariableFeatures(integrated)
integrated <- ScaleData(integrated, vars.to.regress = c("nFeature_RNA", "percent.mt"))
integrated <- RunPCA(integrated)

# Harmony integration correcting for source (+ optionally continuous covariates)
# Split layers by source before integrating
integrated[["RNA"]] <- split(integrated[["RNA"]], f = integrated$source)

# Then run IntegrateLayers
integrated <- IntegrateLayers(
  object = integrated,
  method = HarmonyIntegration,
  orig.reduction = "pca",
  new.reduction = "harmony",
  group.by.vars = "source",
  verbose = FALSE
)

# Rejoin after
integrated <- JoinLayers(integrated)


# Recluster
integrated <- FindNeighbors(integrated, reduction = "harmony", dims = 1:30)
integrated <- FindClusters(integrated, resolution = 0.5, cluster.name = "harmony_reclustered")
integrated <- RunUMAP(integrated, reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony")

DimPlot(integrated, group.by = "harmony_reclustered", reduction = "umap.harmony", label = TRUE)
DimPlot(integrated, group.by = "source", reduction = "umap.harmony", shuffle = T)

POA = readRDS("~/Desktop/nemo.orig_harmony.integration_all_testd_clusters.rds")

Idents(POA) = 'res0.8_50nn_40PC_45LSI'
poa_anchors <- FindTransferAnchors(reference = POA, query = integrated, dims = 1:30,
    reference.reduction = "rnaPCA")
predictions <- TransferData(anchorset = poa_anchors, refdata = POA$res0.8_50nn_40PC_45LSI, dims = 1:30)
integrated <- AddMetaData(integrated, metadata = predictions)
integrated$multiome_prediction = ifelse(
predictions$prediction.score.max<0.9, NA, integrated$predicted.id)


DimPlot(integrated, group.by = "harmony_reclustered", reduction = "umap.harmony", label = TRUE)
DimPlot(integrated, group.by = "multiome_prediction", reduction = "umap.harmony", label = T)
# so 6 is more or less 1, I guess it really is a true cluster, 9 is no longer downstream from rgcs 
# not really sure what is going on with 3
parker = readRDS("~/Desktop/parker_object.rds")

parker_anchors <- FindTransferAnchors(reference = parker, query = integrated, dims = 1:30,
    reference.reduction = "pca")
predictions <- TransferData(anchorset = parker_anchors, refdata = parker$clusters49, dims = 1:30)
integrated <- AddMetaData(integrated, metadata = predictions)
integrated$parker_predictions49 = ifelse(
predictions$prediction.score.max<0.8, NA, integrated$predicted.id)

DimPlot(integrated, group.by = "parker_predictions49", reduction = "umap.harmony", label = T)
# im realizing now that this is equivalent to overtraining lmfao, but it has no idea what 6 is

DimPlot(parker, group.by = 'clusters49', label = T)
