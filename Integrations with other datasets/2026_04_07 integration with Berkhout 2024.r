library(Seurat)
library(tidyverse)
library(biomaRt)

#### Read in Stuff##################


# read in clownfish data
obj = readRDS("C:/Users/Gabe/Desktop/nemo.orig_harmony.integration_all_testd_clusters.rds")

# read in berkhout pvn atlas
pvn = readRDS("C:/Users/Gabe/Downloads/Berkhout2023_PVN_Atlas.rds")

Idents(pvn) = 'annotated'
DimPlot(pvn, reduction = 'mde')

# Read in clownfish - mus lookup table
m_ocellaris = read.csv("Reference/mmusculus_to_aocellaris.csv")


#### Cross species conversions##################
# Subset both objects to common genes
pvn_genes_riken = rownames(pvn)
ocellaris_genes_symbol = rownames(obj)

## Rename mmusculus riken genes to symbol so I can use lookup table
# use biomart
mart_mus <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")

riken_to_symbol <- getBM(
  attributes = c("mgi_symbol", "external_gene_name"),
  values = pvn_genes_riken,
  mart = mart_mus
)

#make naming vector
riken_to_sym_vec <- setNames(
  riken_to_symbol$mgi_symbol,
  riken_to_symbol$external_gene_name
)

# Map mus gene names to symbols
pvn_genes_symbol <- riken_to_sym_vec[pvn_genes_riken]

# Now filter using lookup table
common_mus_symbols <- intersect(pvn_genes_symbol, m_ocellaris$mmusculus_name)

#define common genes among both organisms
mus_symbols_keep = common_mus_symbols

# for ocellaris, keep all genes that are in the mus symbols keep
ocellaris_symbols_keep =  m_ocellaris$aocellaris_name[m_ocellaris$mmusculus_name%in% mus_symbols_keep & m_ocellaris$aocellaris_name!='']
ocellaris_symbols_keep

#now we have to go back to mus and give them real names
# Get the riken names back for the common genes
mus_to_riken_vec <- setNames(
  riken_to_symbol$external_gene_name,
  riken_to_symbol$mgi_symbol
)
riken_symbols_keep <- mus_to_riken_vec[mus_symbols_keep]

# Extract counts matrix and subset to common genes
pvn_counts <- GetAssayData(pvn, layer = "counts")%>%as.matrix()
pvn_counts_sub <- pvn_counts[riken_symbols_keep[riken_symbols_keep%in%rownames(pvn_counts)], ]

# Rename rows from riken to symbol
rownames(pvn_counts_sub) <- (riken_symbols_keep[riken_symbols_keep%in%rownames(pvn_counts)])

# Make new Seurat object
pvn_new <- CreateSeuratObject(counts = pvn_counts_sub)

# Transfer metadata
pvn_new <- AddMetaData(pvn_new, metadata = pvn@meta.data)

#### Make ocellaris have mmusculus names for label transfer#######
# Subset ocellaris to common genes
obj_counts <- GetAssayData(obj, layer = "counts") %>% as.matrix()
obj_counts_sub <- obj_counts[ocellaris_symbols_keep[ocellaris_symbols_keep %in% rownames(obj_counts)], ]

# Build renaming vector: ocellaris -> mus
ocellaris_to_mus <- setNames(
  m_ocellaris$mmusculus_name[m_ocellaris$aocellaris_name %in% ocellaris_symbols_keep & m_ocellaris$aocellaris_name != ''],
  m_ocellaris$aocellaris_name[m_ocellaris$aocellaris_name %in% ocellaris_symbols_keep & m_ocellaris$aocellaris_name != '']
)

# Keep only rows that have a valid mapping
valid_genes <- rownames(obj_counts_sub)[rownames(obj_counts_sub) %in% names(ocellaris_to_mus)]

# Build renaming vector: ocellaris -> mus
ocellaris_to_mus <- setNames(
  m_ocellaris$mmusculus_name[m_ocellaris$aocellaris_name %in% ocellaris_symbols_keep & m_ocellaris$aocellaris_name != ''],
  m_ocellaris$aocellaris_name[m_ocellaris$aocellaris_name %in% ocellaris_symbols_keep & m_ocellaris$aocellaris_name != '']
)

# Drop any entries where mus name is empty or whitespace
ocellaris_to_mus <- ocellaris_to_mus[trimws(ocellaris_to_mus) != '']

# Keep only rows that have a valid mapping
valid_genes <- rownames(obj_counts_sub)[rownames(obj_counts_sub) %in% names(ocellaris_to_mus)]
obj_counts_sub <- obj_counts_sub[valid_genes, ]

rownames(obj_counts_sub) <- ocellaris_to_mus[rownames(obj_counts_sub)]

# Coalesce duplicate mus symbols by summing counts
obj_counts_sub <- rowsum(obj_counts_sub, row.names(obj_counts_sub))

# Make new Seurat object
obj_new <- CreateSeuratObject(counts = obj_counts_sub)
obj_new <- AddMetaData(obj_new, metadata = obj@meta.data)

## Normalization ############
# I want to preserve my clustering for the obj
obj_new@reductions = obj@reductions
obj_new@graphs = obj@graphs

# ok lets do some label transfer

# normalize obj new
obj_new =  NormalizeData(obj_new, normalization.method = "LogNormalize", scale.factor = 10000)
obj_new <- FindVariableFeatures(obj_new, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(obj_new)
obj_new <- ScaleData(obj_new, features = all.genes)
obj_new <- RunPCA(obj_new, features = VariableFeatures(object = obj_new))
#ElbowPlot(obj_new)
Idents(obj_new) = 'res0.8_50nn_40PC_45LSI'
DimPlot(obj_new, reduction = 'harmony_wnn.umap')

# normalize pvn
pvn_new =  NormalizeData(pvn_new, normalization.method = "LogNormalize", scale.factor = 10000)
pvn_new <- FindVariableFeatures(pvn_new, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(pvn_new)
pvn_new <- ScaleData(pvn_new, features = all.genes)
pvn_new <- RunPCA(pvn_new, features = VariableFeatures(object = pvn_new))
ElbowPlot(pvn_new) # woah that is shit 
# lets say 8
pvn_new <- FindNeighbors(pvn_new, dims = 1:8)
pvn_new <- FindClusters(pvn_new, resolution = 0.1)
pvn_new <- RunUMAP(pvn_new, dims = 1:8)

Idents(pvn_new) ='annotated'
DimPlot(pvn_new, reduction = "umap", label = T)
# genuinely how did this happen
# ok I am just bringing their old reductions
pvn_new@reductions$mde = pvn@reductions$mde
DimPlot(pvn_new, reduction = "mde", label = T)

### Label Transfer ####
transfer_anchors <- FindTransferAnchors(
  reference            = pvn_new,
  query                = obj_new,
  dims                 = 1:30,
  reference.reduction  = "pca",
  normalization.method = "LogNormalize",
  features             = intersect(VariableFeatures(pvn_new),
                                   VariableFeatures(obj_new))  # shared HVGs only
)

predictions <- TransferData(
  anchorset = transfer_anchors,
  refdata   = list(
    annotated           = pvn_new$annotated
  ),
  dims = 1:30
)

obj_new$predictions = predictions$predicted.id
obj_new$prediction.score.max = predictions$prediction.score.max

DimPlot(obj_new, group.by = "predictions", label = TRUE, reduction = 'harmony_wnn.umap') 

FeaturePlot(obj_new, features = "prediction.score.max", reduction = 'harmony_wnn.umap')
hist(obj_new$prediction.score.max, breaks = 50)

obj_new$predictions_filtered <- ifelse(
  obj_new$prediction.score.max > 0.5,
  obj_new$predictions,
  "low_confidence"
)

DimPlot(obj_new, group.by = "predictions_filtered", label = TRUE, reduction = 'harmony_wnn.umap') 

### heatmap #####
# predictions is a list with a data frame, grab it directly
score_df_raw <- predictions

score_cols <- grep("^prediction.score.", colnames(score_df_raw), value = TRUE)
score_cols <- score_cols[score_cols != "prediction.score.max"]

score_df <- score_df_raw %>%
  mutate(cluster = obj_new$res0.8_50nn_40PC_45LSI) %>%
  dplyr::select(cluster, all_of(score_cols)) %>%
  group_by(cluster) %>%
  summarise(across(everything(), mean)) %>%
  column_to_rownames("cluster")

colnames(score_df) <- gsub("prediction.score.", "", colnames(score_df))

pheatmap::pheatmap(
  as.matrix(score_df),
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  color = colorRampPalette(c("white", "darkred"))(100),
  fontsize = 8,
  main = "Mean prediction score per cluster"
)
# literally no confident clusters lmao













