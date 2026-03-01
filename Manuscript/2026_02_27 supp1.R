
#supp fig 1
#Dotplot 
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
multiome_object <- readRDS("~/Desktop/nemo.orig_harmony.integration_all_testd_clusters.rds") 
Idents(multiome_object) = 'res0.8_50nn_40PC_45LSI'
`%notin%` <- Negate(`%in%`)

#### By marker genes #####


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
  'LOC111577263', #brain aromatase - radial glia
  'gfap', #astrocyte marker
  'crocc2', #ependymal cell marker
  'mbpa', #oligo marker
  'cspg4', #OPC marker
  'p2ry12', #microglia marker
  'ptprc', #leukocyte marker
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
  'avpr2aa'

) )


marker_gene_plot <- DotPlot(object = multiome_object, 
                 group.by = "res0.8_50nn_40PC_45LSI", 
                 features = markers,
        dot.min = 0.1
) + 
  coord_flip()+
  scale_size(range = c(0,2))+
  theme(axis.text.x = element_text(angle = -90))+
  scale_x_discrete(labels= c(
'elavl3',#neuron
'gad2',#GABA 
'gad1', #gad1
'vglut2.1', #vgat2.1
'slc17a6b', #vglut
'slc17a7a', #vglut1
'sst1.1', #interneuron marker
'slc18a3b', #ach marker
  'hmx2',
  'hmx3a',
'cyp19a1b', #brain aromatase - radial glia
'gfap', #astrocyte marker
'crocc2', #ependymal cell marker
'mbpa', #oligo marker
'cspg4', #OPC marker
'p2ry12', #microglia marker
'ptprc', #leukocyte marker
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
'npy4r', #npy4r
'esr1',
'esr2a',
'esr2b',
'ar',
'ccka', #ccka
'cckb',
'cckar',
'cckbra',
'cckbrb',
'gal',
'galr1a',
'oxt',
  'oxtrb',
'avp',
'avpr2aa'
))+
theme(axis.title.y = element_blank())+
  labs(y = 'Cluster')+
  theme(axis.text.x = element_text(size = 10, vjust =0.4),
        axis.text.y = element_text(size = 8),
        legend.text = element_text(size=10),
        legend.title  = element_text(size=12))
marker_gene_plot

ggsave(plot = marker_gene_plot,
       file = "marker_gene_plot.svg",
       device = "svg",
       units = "in",
       width = 6,
       height = 4.5,
       path = "Manuscript/Plots/Supp 1")




