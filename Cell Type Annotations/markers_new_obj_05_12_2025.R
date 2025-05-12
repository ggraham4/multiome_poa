{
  library(parallel)
  library(clusterProfiler)
  library(blme)
  library(Seurat)
  library(tidyverse)
  library(tidyr)
  library(lme4)
  library(dplyr)
  library(MASS)
  library(SeuratObject)
  library(Signac)
  library('glmGamPoi')
  library(scran)
  library(parallel)
  library(factoextra)
  library(readxl)
  library(factoextra)
  library(forcats)
  library(ggrepel)
  library(biomaRt)
  library(openxlsx)

}

obj <- readRDS('~/Desktop/optimal_clustering_rna_only.rds')

markers <- unique(c(
  'elavl3',#neuron
  'gad2',#GABA 
  'LOC111588076', #gad1
  'LOC111584103', #vglut2.1
  'slc17a6b', #vglut
  'slc17a7a', #vglut1
  'sst1.1', #interneuron marker
    'slc18a3b', #ach marker
  'hmx2',
  'hmx3a',
  'kiss1',
  'kiss1ra',
  'kiss1rb',
  'tac1',
  'tac3a',
  'tacr1a',
  'tacr3a',
  'tacr3l',
  'npy',
  'npy8ar',
  'npy2rl',
  'LOC111581412', #npy4r
  'esr1',
  'esr2a',
  'esr2b',
  'ar',
  'LOC111562384', #ccka
  'cckb',
  'cckar',
  'cckbra',
  'cckbrb',
  'gal',
  'galr1a',
  'oxt',
  'oxtrb',
  'avp',
  'avpr2aa',
  'LOC111577263', #brain aromatase - radial glia
  'gfap', #astrocyte marker
  'crocc2', #ependymal cell marker
  'mbpa', #oligo marker
  'cspg4', #OPC marker
  'p2ry12', #microglia marker
  'ptprc', #leukocyte marker
  'col15a1b',
  'gdpd5a'

) )


marker_gene_plot <- DotPlot(object = obj, 
                 features = markers,
        dot.min = 0.1
) + 
  coord_flip()+
  scale_size(range = c(0,2))+
  theme(axis.text.x = element_text(angle = -90))
marker_gene_plot

#the hell is going on with 24

DimPlot(obj, label = T)

markers_24 = FindMarkers(obj, 24)

table_data = as.data.frame.matrix(table(obj$res0.8_50nn_40PC_45LSI, obj$individual))

table_data[25,]

#this really seems like an immature neuron IDK what the fuck this is

library(CytoTRACE)

`%notin%`  = Negate(`%in%`)
neuron_only = subset(obj,res0.8_50nn_40PC_45LSI %notin% c(1, #rgc
                                                          2, #olig
                                                          11, #micro
                                                          13, #opc
                                                          20, # leuko
                                                          26 #another glia
                                                          
                                                          ))

cyto = CytoTRACE(as.matrix(neuron_only@assays$RNA$counts))
neuron_only$cyto = cyto$CytoTRACE

ggplot(neuron_only@meta.data, aes(x = res0.8_50nn_40PC_45LSI, y = cyto))+
  geom_boxplot()
#24 is NOT an immature neuron?? But 22 is?

DimPlot(neuron_only, label = T)

markers_22 = FindMarkers(obj, 22)
#22 likely 27

DimPlot(obj, label = T, group.by = 'gabe_celltype_ids') #wtf does NA mean
# it looks like 24 is just straight up new

DimPlot(obj, label = T) #wtf does NA mean

### i think its label transfer time

oldish_obj = readRDS("/Volumes/jrhodes/Fish Lab/Experiments/sex change single nuc POA/R/RNA Object.rds")
Idents(oldish_obj) ='harmony.wnn_res0.4_clusters'

anchors <- FindTransferAnchors(reference = oldish_obj, query = obj,
                               reference.reduction = "pca")

predictions <- TransferData(anchorset = anchors, refdata = oldish_obj$harmony.wnn_res0.4_clusters)
obj <- AddMetaData(obj, metadata = predictions)

DimPlot(obj, label = T, group.by = 'predicted.id')
# its 15??


