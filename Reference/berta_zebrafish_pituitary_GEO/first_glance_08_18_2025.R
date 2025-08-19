# ---- packages ----
library(GEOquery)
library(Seurat)
library(Matrix)
library(dplyr)
library(stringr)

# -----------------------------
# 1) Download supplementary files from GEO
# -----------------------------
setwd('A:/')
acc <- "GSE272806"
getGEOSuppFiles(acc, makeDirectory = TRUE)   # creates ./GSE272806/
supp_dir <- file.path(getwd(), acc)

# Helper: find files by the names shown on GEO
f_genes   <- list.files(supp_dir, pattern = "Integrated_filtered_genes.tsv.gz$",   full.names = TRUE)
f_barc    <- list.files(supp_dir, pattern = "Integrated_filtered_barcodes.tsv.gz$",full.names = TRUE)
f_mtx     <- list.files(supp_dir, pattern = "Integrated_filtered_matrix.mtx.gz$",  full.names = TRUE)
f_meta    <- "A:/GSE272806/GSE272806_Integrated_filtered_metadata.csv.gz"

stopifnot(length(f_genes)==1, length(f_barc)==1, length(f_mtx)==1, length(f_meta)==1)

# -----------------------------
# 2) Read the MTX bundle into Seurat
# -----------------------------
# Read features/genes (zebrafish). These files can be 1–3 columns; use the first non-empty as ID.
genes <- read.delim(f_genes, header = FALSE, stringsAsFactors = FALSE)
# pick an ID column (prefer Ensembl if present)
id_col <- which.max(colSums(genes != ""))   # crude but robust
gene_ids <- make.unique(genes[[id_col]])    # avoid duplicate rownames

barcodes <- read.delim(f_barc, header = FALSE, stringsAsFactors = FALSE)[[1]]
mtx <- readMM(f_mtx)

# Sanity: align dimnames
stopifnot(nrow(mtx) == length(gene_ids), ncol(mtx) == length(barcodes))
rownames(mtx) <- gene_ids
colnames(mtx) <- barcodes

# Create Seurat object for GEO
geo <- CreateSeuratObject(counts = mtx, project = "GSE272806")

# Add GEO-provided per-cell metadata (align by barcode)
geo_meta <- read.csv(f_meta)
# assume metadata has a column with cell barcodes; guess column by exact or near match
bc_col <- which(colnames(geo_meta) %in% c("barcode","barcodes","cell","cell_id","Cell","CellID"))
if (length(bc_col)==0) {
  # try to find a column that matches most barcodes
  bc_col <- which.max(colSums(sapply(geo_meta, function(x) as.character(x) %in% colnames(geo))))
}
geo_meta$barcode_tmp <- as.character(geo_meta[[bc_col]])
geo_meta <- geo_meta[match(colnames(geo), geo_meta$barcode_tmp), ]

geo <- AddMetaData(geo, metadata = geo_meta)


#make DimPlot  - seurat vignette
# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
geo[["percent.mt"]] <- PercentageFeatureSet(geo, pattern = "^MT-")
# Visualize QC metrics as a violin plot
VlnPlot(geo, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.

plot1 <- FeatureScatter(geo, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(geo, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
geo <- subset(geo, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
geo <- NormalizeData(geo, normalization.method = "LogNormalize", scale.factor = 10000)
geo <- FindVariableFeatures(geo, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(geo), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(geo)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2
all.genes <- rownames(geo)
geo <- ScaleData(geo, features = all.genes)
geo <- RunPCA(geo, features = VariableFeatures(object = geo))
# Examine and visualize PCA results a few different ways
print(geo[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(geo, dims = 1:2, reduction = "pca")
DimPlot(geo, reduction = "pca") + NoLegend()
DimHeatmap(geo, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(geo, dims = 1:15, cells = 500, balanced = TRUE)
ElbowPlot(geo) #20 was used in manuscript

geo <- FindNeighbors(geo, dims = 1:20)
geo <- FindClusters(geo, resolution = 0.4)
# Look at cluster IDs of the first 5 cells
head(Idents(geo), 5)
# Look at cluster IDs of the first 5 cells
head(Idents(geo), 5)
geo <- RunUMAP(geo, dims = 1:20)

DimPlot(geo, label = T, group.by = 'celltaxonomy')
#saveRDS(geo, 'A:/GSE272806/seurat_object.rds')

markers= FindAllMarkers(geo)
markers_signif = subset(markers, p_val_adj<0.05)

### try transfer anchors
obj_6_only <- readRDS( 'A:/WGCNA_6/wgcna_seurat_obj_08_12_2025.rds')

anchors <- FindTransferAnchors(reference = geo, query = obj_6_only,
                               reference.reduction = "pca")

predictions <- TransferData(anchorset = anchors, refdata = geo$RNA_snn_res.0.4)
obj_6_only <- AddMetaData(obj_6_only, metadata = predictions)

DimPlot(obj_6_only, reduction = 'harmony_wnn.umap', group.by = 'predicted.id', label =T)

### aligns most strongly to 0 it seems

table(obj_6_only$sub.cluster, obj_6_only$predicted.id)

table(geo$cluster, geo$RNA_snn_res.0.4)
#cluster 10

FeaturePlot(obj_6_only,'her6', reduction = 'harmony_wnn.umap')

#allegedly melanotropes, though Im suspicious of how biologicaly possible that is, cause POA innervates pit