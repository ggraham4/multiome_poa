library(Seurat)
library(Signac)
library(ggplot2)
library(tidyverse)
library(lme4)

obj <- readRDS("~/Desktop/nemo.orig_harmony.integration_all_testd_clusters.rds")

sub_obj <- subset(obj, res0.8_50nn_40PC_45LSI %in% c(1,6))

###############################################################################
# RNA preprocessing
###############################################################################

DefaultAssay(sub_obj) <- "RNA"

sub_obj <- NormalizeData(sub_obj)
sub_obj <- FindVariableFeatures(sub_obj, selection.method = "vst", nfeatures = 2000)

sub_obj <- ScaleData(sub_obj)

sub_obj <- RunPCA(
  sub_obj,
  features = VariableFeatures(sub_obj),
  reduction.name = "pca"
)

ElbowPlot(sub_obj)

###############################################################################
# ATAC preprocessing
###############################################################################

DefaultAssay(sub_obj) <- "ATAC"

sub_obj <- RunTFIDF(sub_obj)

sub_obj <- FindTopFeatures(
  sub_obj,
  min.cutoff = "q0"
)

sub_obj <- RunSVD(
  sub_obj,
  reduction.name = "lsi"
)

DepthCor(sub_obj)

###############################################################################
# WNN integration
###############################################################################

# Typical choices:
# RNA PCs 1:30
# ATAC LSI 2:30 (skip first component because it often reflects depth)

sub_obj <- FindMultiModalNeighbors(
  sub_obj,
  reduction.list = list("pca", "lsi"),
  dims.list = list(1:30, 2:30),
  modality.weight.name = "RNA.weight"
)

###############################################################################
# Clustering on WNN graph
###############################################################################

sub_obj <- FindClusters(
  sub_obj,
  graph.name = "wsnn",
  algorithm = 3,
  resolution = 0.2
)

###############################################################################
# UMAP on WNN graph
###############################################################################

sub_obj <- RunUMAP(
  sub_obj,
  nn.name = "weighted.nn",
  reduction.name = "wnn.umap",
  reduction.key = "wnnUMAP_"
)

DimPlot(
  sub_obj,
  reduction = "wnn.umap",
  label = TRUE
)

####
DefaultAssay(sub_obj) <- "RNA"
Marks = FindAllMarkers(sub_obj, only.pos = T)

DotPlot(sub_obj, c(
  'slc17a6b',
  'gad2',
  'elavl3',
  'LOC111577263',
  's100b',
  'glul',
  'fgfr3',
  'sox4a',
  'gja1b',
  'pcna',
  'mki67',
  'gfap',
  'hmx2',
  'hmx3a',
  'oxt',
  'avp',
  'cckb',
  'gal',
    'LOC111571064',
  'crhb'
))+
  coord_flip()

marks_0_2 = FindMarkers(sub_obj, ident.1 = '0', ident.2 = '2')
marks_8 = FindMarkers(sub_obj,'8')
marks_3 = FindMarkers(sub_obj,'3')
marks_4 = FindMarkers(sub_obj,'4')

labs = c(
  '0' ='RGC_arob',
  '1' = 'POA_cckb',
  '2' = 'RGC_arob_lrp4',
  '3' ='RGC_quiescent',
  '4'= 'RGC_wnt',
  '5' = 'POA_gal',
  '6' = 'POA_oxt',
  '7' = 'RGC_gfap_ptch1',
  '8' ='RGC_aqp',
  '9'='POA_gnrh1_roraa',
  '10' = 'RGC_gfap_grid2',
  '11' = 'POA_avp',
  '12' = 'POA_gnrh1_th'
)

sub_obj@meta.data$named_clust = labs[sub_obj@meta.data$wsnn_res.0.2]%>%unname()

#marks_gnrh = FindMarkers(sub_obj, ident.1 = 'POA_gnrh1_1', ident.2 = 'POA_gnrh1_2')

#marks_gfap = FindMarkers(sub_obj, ident.1 = 'RGC_gfap_1', ident.2 = 'RGC_gfap_2')

Idents(sub_obj)= 'named_clust'
DimPlot(sub_obj, reduction = 'wnn.umap', label = T)

Cyto = CytoTRACE::CytoTRACE(sub_obj@assays$RNA$data%>%as.matrix())$CytoTRACE
sub_obj$Cyto =Cyto
FeaturePlot(sub_obj, 'Cyto', reduction = 'wnn.umap', label = T)

DotPlot(sub_obj, 'Cyto')
FeaturePlot(sub_obj, c('sox4b',
                       'pcna',
                       'LOC111574271'), 
            reduction = 'wnn.umap',
            label = F, 
            order =T)


DotPlot(sub_obj, c('sox4b',
                       'pcna',
                       'LOC111574271'))

library(AUCell)

# Create gene set
gene_set <- list(sox4b_pcna_module = c('sox4b', 'pcna', 'LOC111574271'))

# Get expression matrix
expr_matrix <- sub_obj@assays$RNA$data %>% as.matrix()

# Build rankings
cell_rankings <- AUCell_buildRankings(expr_matrix, plotStats = FALSE)

# Calculate AUC
cell_AUC <- AUCell_calcAUC(gene_set, cell_rankings)

# Add to metadata
sub_obj$sox4b_pcna_AUC <- as.numeric(getAUC(cell_AUC)["sox4b_pcna_module", ])

# DotPlot
DotPlot(sub_obj, features = "sox4b_pcna_AUC")
sub_obj@meta.data$Status = factor(sub_obj@meta.data$Status,
                                  levels = c('NRM',
                                             'M',
                                             'D',
                                             'E',
                                             'NF',
                                             'F'))
ggplot(sub_obj@meta.data%>%
  group_by(individual, Status,named_clust )%>%
  summarize(mean_mod = mean(sox4b_pcna_AUC)),
  
  aes(x= Status, y = mean_mod))+
  geom_boxplot()+
  geom_point()+
facet_wrap(~named_clust, scales = 'free')

te =sub_obj@meta.data%>%
  group_by(individual, Status,named_clust )%>%
  summarize(mean_mod = mean(sox4b_pcna_AUC))%>%
  subset(named_clust == 'POA_cckb'&
           Status != 'NRM')
model = lm(mean_mod~Status, data = te)
anova(model, test= 'Chisq')

teer = sub_obj@meta.data%>%
    subset(named_clust == 'POA_cckb'&
           Status != 'NRM')

model2 = lmer(sox4b_pcna_AUC~Status+(1|individual), data = teer)
car::Anova(model2, 3)
# also no
