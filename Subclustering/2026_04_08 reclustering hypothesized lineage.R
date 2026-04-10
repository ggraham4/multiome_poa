### im gonna waste some time and do a whole re analysis of cluster6
library(Seurat)
library(tidyverse)

# read in
obj <- readRDS("~/Desktop/nemo.orig_harmony.integration_all_testd_clusters.rds")

# subset 
obj_6 = subset(obj, res0.8_50nn_40PC_45LSI ==6)

# renormalize and recluster
obj_6 =  NormalizeData(obj_6, normalization.method = "LogNormalize", scale.factor = 10000)
obj_6 <- FindVariableFeatures(obj_6, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(obj_6)
obj_6 <- ScaleData(obj_6, features = all.genes)
obj_6 <- RunPCA(obj_6, features = VariableFeatures(object = obj_6))
ElbowPlot(obj_6)
# 4
obj_6 <- FindNeighbors(obj_6, dims = 1:4)
obj_6 <- FindClusters(obj_6, 
                      resolution = 0.2, 
                      graph.name = 'RNA_nn')
obj_6 <- RunUMAP(obj_6, dims = 1:4)
DimPlot(obj_6, reduction = "umap")
# Fuck it this is as good as its gonna get
Idents(obj_6)= 'RNA_nn_res.0.2'
marks_obj_6 = FindAllMarkers(obj_6, only.pos = T)

DimPlot(subset(obj_6, Status %in% c('M','D','F')), group.by = 'Status', reduction = 'umap')

# plot markers
# Get top 4 markers per cluster
top4 <- marks_obj_6 %>%
  group_by(cluster) %>%
  filter(p_val_adj < 0.05) %>%
  slice_max(order_by = avg_log2FC, n = 4) %>%
  ungroup()

genes_to_plot <- unique(top4$gene)

# Dot plot
DotPlot(obj_6,
        features   = genes_to_plot,
        group.by   = "RNA_nn_res.0.2")+
  coord_flip()

DotPlot(obj_6,
        features   = c('cckb',
                       'LOC111571064',
                       'esr2b',
                       'tacr3a'),
        group.by   = "RNA_nn_res.0.2")+
  coord_flip()

obj_6$annot = ifelse(obj_6$RNA_nn_res.0.2 == 3, 'gnrh_ndnf', NA)
obj_6$annot = ifelse(obj_6$RNA_nn_res.0.2 == 2, 'avp_sim1a', obj_6$annot)
obj_6$annot = ifelse(obj_6$RNA_nn_res.0.2 == 1, 'thrhra_gal', obj_6$annot)
obj_6$annot = ifelse(obj_6$RNA_nn_res.0.2 == 0, 'bcl11aa_SR', obj_6$annot) # also lhx1 + so seems to be immature?

Idents(obj_6) = 'annot'

### CytoTRACE####
library(CytoTRACE)
obj_6$CytoTRACE = CytoTRACE(obj_6@assays$RNA$counts%>%as.matrix())$CytoTRACE
FeaturePlot(obj_6, 'CytoTRACE', reduction = 'umap')

cyt_plt = obj_6@meta.data%>%
  group_by(individual, Status, annot)%>%
  summarize(mean_cyto = mean(CytoTRACE))%>%
  subset(Status %in% c('M','D','F'))
cyt_plt$Status = factor(cyt_plt$Status, levels = c('M','D','F'))

ggplot(cyt_plt, aes(x = annot,y = mean_cyto, color = Status))+
  geom_boxplot()+
  geom_point(position = position_dodge(0.75))

DimPlot(obj_6, reduction = "umap", label = T)

### Proportions ####
prop <- obj_6@meta.data %>%
  group_by(Status, individual, annot) %>%
  summarise(n = n()) %>%
  group_by(individual) %>%
  mutate(prop = n / sum(n))%>%
    subset(Status %in% c('M','D','F'))
prop$Status = factor(prop$Status, levels = c('M','D','F'))

ggplot(prop, aes(x = annot, y = prop, color = Status)) +
  geom_boxplot()+
    geom_point(position = position_dodge(0.75))

# cckb in the bcl11 cluster

#### If I recluster with just radial glia, 9, and 6, is there a lineage?
obj_lineage = subset(obj, res0.8_50nn_40PC_45LSI %in% c(6, 1, 9))

# renormalize and recluster
obj_lineage =  NormalizeData(obj_lineage, normalization.method = "LogNormalize", scale.factor = 10000)
obj_lineage <- FindVariableFeatures(obj_lineage, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(obj_lineage)
obj_lineage <- ScaleData(obj_lineage, features = all.genes)
obj_lineage <- RunPCA(obj_lineage, features = VariableFeatures(object = obj_lineage))
ElbowPlot(obj_lineage)
# 5
obj_lineage <- FindNeighbors(obj_lineage, dims = 1:5)
obj_lineage <- FindClusters(obj_lineage, 
                      resolution = 0.2, 
                      graph.name = 'RNA_nn')
obj_lineage <- RunUMAP(obj_lineage, dims = 1:5)
DimPlot(obj_lineage, reduction = "umap")

Idents(obj_lineage) = 'RNA_nn_res.0.2'
marks_lineage = FindAllMarkers(obj_lineage)
FeaturePlot(obj_lineage, 'LOC111577263', reduction = 'umap')
DimPlot(obj_lineage, reduction = 'harmony_wnn.umap', label = T)
marks_lineage = FindAllMarkers(obj_lineage)

obj_lineage$annotation = ifelse(obj_lineage$RNA_nn_res.0.2 == 0, '1_RGC_1', NA)
obj_lineage$annotation = ifelse(obj_lineage$RNA_nn_res.0.2 == 1, '1_RGC_2', obj_lineage$annotation)
obj_lineage$annotation = ifelse(obj_lineage$RNA_nn_res.0.2 == 2, '9_1', obj_lineage$annotation)
obj_lineage$annotation = ifelse(obj_lineage$RNA_nn_res.0.2 == 3, '6_POA_1', obj_lineage$annotation)
obj_lineage$annotation = ifelse(obj_lineage$RNA_nn_res.0.2 == 4, '6_POA_2', obj_lineage$annotation )
obj_lineage$annotation = ifelse(obj_lineage$RNA_nn_res.0.2 == 5, '9_2', obj_lineage$annotation )
obj_lineage$annotation = ifelse(obj_lineage$RNA_nn_res.0.2 == 6, '6_POA_3', obj_lineage$annotation )

Idents(obj_lineage) = 'annotation'
DimPlot(obj_lineage, reduction = 'umap', label = T)
DimPlot(subset(obj_lineage, Status %in% c('M','D','F')), reduction = 'umap', label = F, group.by = 'Status')

DimPlot(obj_lineage, reduction = 'atacUMAP', label = T)
DimPlot(obj_lineage, reduction = 'atacUMAP', label = F)
DimPlot(subset(obj_lineage, Status %in% c('M','D','F')), reduction = 'atacUMAP', label = F, group.by = 'Status')


obj_lineage$CytoTRACE = CytoTRACE(obj_lineage@assays$RNA$counts%>%as.matrix())$CytoTRACE
FeaturePlot(obj_lineage, 'CytoTRACE', reduction = 'umap')
#oh yeah I forgot RGCs dont express much lmao oh no

prop2 <- obj_lineage@meta.data %>%
  group_by(Status, individual, annotation) %>%
  summarise(n = n()) %>%
  group_by(individual) %>%
  mutate(prop = n / sum(n))%>%
    subset(Status %in% c('M','D','F'))
prop2$Status = factor(prop2$Status, levels = c('M','D','F'))



# well I dont understand so what can you do

cyt2 <- obj_lineage@meta.data %>%
  group_by(Status, individual, annotation) %>%
  mutate(cyt = mean(CytoTRACE))%>%
    subset(Status %in% c('M','D','F'))
cyt2$Status = factor(cyt2$Status, levels = c('M','D','F'))

ggplot(cyt2, aes(x = annotation, y = cyt, color = Status)) +
  geom_boxplot()+
    geom_point(position = position_dodge(0.75))

ggplot(prop2, aes(x = annotation, y = prop, color = Status)) +
  geom_boxplot()+
    geom_point(position = position_dodge(0.75))


#SUPPORT OH YEAH LETS FUCKING GO
Idents(obj) = 'res0.8_50nn_40PC_45LSI'
FeaturePlot(obj, 'sox2', reduction = 'harmony_wnn.umap')

FeaturePlot(obj_lineage, 'sox2', reduction = 'umap')
DimPlot(obj_lineage, label = T, reduction = 'umap')

FeaturePlot(obj, 'pcna', reduction = "harmony_wnn.umap")
FeaturePlot(obj, 'LOC111574271', reduction = "harmony_wnn.umap") #dcx

DotPlot(obj, 'LOC111574271') 
DimPlot(obj, label = T, reduction = "harmony_wnn.umap")
#