library(Seurat)
library(tidyverse)

# read in hypomap
hypomap = readRDS('C:/Users/Gabe/Downloads/hypoMap.rds')
#read in clonwfish data
obj = readRDS("C:/Users/Gabe/Desktop/nemo.orig_harmony.integration_all_testd_clusters.rds")

# Read in clownfish - mus lookup table
m_ocellaris = read.csv("Reference/mmusculus_to_aocellaris.csv")

# convert ocellaris names to mus names, subset to common features ####
# get mus genes
hypomap_genes = rownames(hypomap)

# extract ocellaris matrix
ocellaris_counts_matrix = obj@assays$RNA$counts

# get ocellaris genes
ocellaris_genes = rownames(ocellaris_counts_matrix)

# match ocellaris genes to mus genes, filter to those that exist in both the object and with a mus homolog
m_ocellaris_filt <- m_ocellaris %>%
  filter(aocellaris_name %in% ocellaris_genes) %>%   # keep only genes present in your data
  filter(mmusculus_name != "" & !is.na(mmusculus_name)) %>%# remove empty strings and NAs
  filter(mmusculus_name %in%hypomap_genes)  # keep genes that are in hypomap

#subset the counts matrix to those genes
ocellaris_counts_mus_match <- ocellaris_counts_matrix[m_ocellaris_filt$aocellaris_name, ]

#assign those genes as rownames
rownames(ocellaris_counts_mus_match) <- m_ocellaris_filt$mmusculus_name

#make dense matrix
ocellaris_counts_mus_match_dense <- as.matrix(ocellaris_counts_mus_match)

#coalesce across common rownames
coalesced_matrix <- rowsum(ocellaris_counts_mus_match_dense, group = rownames(ocellaris_counts_mus_match_dense))

### make new clownfish object ####
mus_clown = CreateSeuratObject(coalesced_matrix,
                               meta.data = obj@meta.data)

# recluster
mus_clown =  NormalizeData(mus_clown, normalization.method = "LogNormalize", scale.factor = 10000)
mus_clown <- FindVariableFeatures(mus_clown, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(mus_clown)
mus_clown <- ScaleData(mus_clown, features =  VariableFeatures(mus_clown))
mus_clown <- RunPCA(mus_clown, features = VariableFeatures(object = mus_clown))

# add clustering back 
mus_clown@reductions$harmony_wnn.umap = obj@reductions$harmony_wnn.umap

# use the same clusters we have calculated previously, but different UMAP cause new features, though it literally does not matter
Idents(mus_clown) = 'res0.8_50nn_40PC_45LSI'
DimPlot(mus_clown, reduction = 'harmony_wnn.umap')

### Do the same w hypomap ###
DimPlot(hypomap, label =T)
hypomap =  NormalizeData(hypomap, normalization.method = "LogNormalize", scale.factor = 10000)
hypomap <- FindVariableFeatures(hypomap, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(hypomap)
hypomap <- ScaleData(hypomap, features = VariableFeatures(object = hypomap))
hypomap <- RunPCA(hypomap, features = VariableFeatures(object = hypomap))

### Label Transfer ####
transfer_anchors <- FindTransferAnchors(
  reference            = hypomap,
  query                = mus_clown,
  dims                 = 1:30,
  reference.reduction  = "pca",
  normalization.method = "LogNormalize",
  features             = intersect(VariableFeatures(hypomap),
                                   VariableFeatures(mus_clown)))  


predictions <- TransferData(
  anchorset = transfer_anchors,
  refdata   = list(
    annotated           = hypomap$Region_summarized
  ),
  dims = 1:30
)

all_predictions <- do.call(cbind, predictions)
mus_clown <- AddMetaData(mus_clown, metadata = all_predictions)

tab = table(mus_clown$res0.8_50nn_40PC_45LSI, mus_clown$predicted.id)%>%as.data.frame.matrix()%>%t()
tab = tab/rowSums(tab)

library(pheatmap)
pheatmap(scale(tab))
 # still no dice
