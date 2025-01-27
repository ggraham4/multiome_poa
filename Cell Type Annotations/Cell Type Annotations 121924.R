#Cell type annotations 121924
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
  #library(CytoTRACE)
 # SCTRransform_mean_plot <- readRDS("R/Gabe/SCTRransform_mean_plot.rds")
  #mac.neg.bin <- readRDS(file = 'R/Gabe/mac.neg.bin.rds')
  library('glmGamPoi')
  library(scran)
  library(parallel)
  library(factoextra)
  library(readxl)
  library(factoextra)
  library(forcats)
  library(ggrepel)
  library(biomaRt)
  #mean_cell <- readRDS('R/Gabe/mean_cell.rds')
  library(openxlsx)
  #clown_go <- readRDS('R/Gabe/clown_go.rds')
  
}
multiome_object <- readRDS('/Volumes/jrhodes/Fish Lab/Experiments/sex change single nuc POA/R/RNA Object.rds')
Idents(multiome_object) <- "harmony.wnn_res0.4_clusters"


#### Find anchorsets and transfer labels #####
previous_object <- readRDS('/Volumes/jrhodes/Fish Lab/Experiments/snRNAseq POA sex differences/analysis/data/anenomefish_clustered_061821.rds')

#### Clusters 49 ####
Idents(previous_object) <- 'clusters49'
anchors <- FindTransferAnchors(reference = previous_object, query = multiome_object,
                               reference.reduction = "pca")

predictions <- TransferData(anchorset = anchors, refdata = previous_object$clusters49)
multiome_object <- AddMetaData(multiome_object, metadata = predictions)

Idents(multiome_object) <- 'predicted.id'
DimPlot(multiome_object, label = T)

prev_labels <- c(
  "0" = c('10, 25'),
  "1" = 25,
  "2" = 15,
  "3" = 39,
  "4" = 14,
  "5" = as.character('12, 11, 5, 6, 8'),
  "6" = as.character(c('31, 17, 26, 47, 29')),
  "7" = 34,
  "8" = as.character(c('25, 10')),
  "9" = 34,
  "10" = 34,
  "11" = 38,
  "12" = 37,
  "13" = 26,
  "14" = 22,
  "15" = as.character(c('0, 2, 25')),
  "16" = 33,
  "17" = 25,
  "18" = 28,
  "19" = 19,
  "20" = as.character(c('21, 36')),
  "21" = 45,
  "22" = 44,
  "23" = as.character(c('19, 25')),
  "24" = 34,
  "25" = 48,
  "26" = 40,
  "27" = 10,
  "28" = 10,
  "29" = 10,
  "30" = 35,
  "31" = 25
)

#### By marker genes #####

#convert from percula gene names to ocellaris
gene_names <- readRDS('/Volumes/jrhodes/Fish Lab/Experiments/sex change single nuc POA/R/old_nemo_ensembl_new_nemo_ncbi_gene_conversion_v2.rds')
percula_to_ocellaris <- function(gene){
converted_name <- gene_names$ncbi_final_gene_names[gene_names$old_nemo_gene_names==gene]
return(converted_name)
  }


Idents(multiome_object) <- "harmony.wnn_res0.4_clusters"
markers <- unique(c('gad2',#GABA 
             'LOC111588076', #gad1
             'LOC111584103', #vgat2.1
             'slc17a6b',
             'slc17a7a', #vglut1
             percula_to_ocellaris('CYP19A1'),
             percula_to_ocellaris('gfap'),
             percula_to_ocellaris('crocc2'),
             percula_to_ocellaris('slc6a11b'),
             percula_to_ocellaris('glula'),
             percula_to_ocellaris('mbpa'),
             percula_to_ocellaris('cspg4'),
             percula_to_ocellaris('p2ry12'),
             percula_to_ocellaris('vegfd'),
             percula_to_ocellaris('chat'),
             percula_to_ocellaris('slc18a3b'),
             percula_to_ocellaris('th'),
             percula_to_ocellaris('slc6a3'),
             percula_to_ocellaris('dbh'),
             percula_to_ocellaris('tph1'),
             percula_to_ocellaris('slc6a4b'),
             percula_to_ocellaris('gnrh1'),
             percula_to_ocellaris('kiss1'),
             percula_to_ocellaris('galn'),
             percula_to_ocellaris('npy'),
             percula_to_ocellaris('oxt'),
             percula_to_ocellaris('oxtr'),
             percula_to_ocellaris('avp'),
             percula_to_ocellaris('avpr'),
             percula_to_ocellaris('pdyn'),
             percula_to_ocellaris('penka'),
             percula_to_ocellaris('tac1'),
             percula_to_ocellaris('tac3a'),
             percula_to_ocellaris('ccka'),
             percula_to_ocellaris('cckb'),
             percula_to_ocellaris('scg2b'),
             percula_to_ocellaris('adcyap1b'),
             percula_to_ocellaris('hmx3a'),#POA
             percula_to_ocellaris('hmx2'),#POA,
             percula_to_ocellaris('crhb'),#POA,
             percula_to_ocellaris('trh'),#POA,
              percula_to_ocellaris('crh'),#POA,
             percula_to_ocellaris('ar'),#POA,
             percula_to_ocellaris('esr2a'),#POA,
             percula_to_ocellaris('esr2b'),#POA,
             percula_to_ocellaris('pgr'),#POA,
             percula_to_ocellaris('nr3c1'),#POA,
          'otpa',
          'ccka',
          'cckbra',
          'elavl3',
                       percula_to_ocellaris('ptprc'),
                    'rbm47',
          'gdpd5a',
          'col15a1b',
          'flvcr2b',
          'slc13a4',
                    'slc13a2',
          'cbln2',
          'tent5aa',
          'kiss1',
          'tac1',
          'tac3a',
          'ar',
          'galr'

) )


plot  <- DotPlot(object = multiome_object, 
                 group.by = "harmony.wnn_res0.4_clusters", 
                 features = markers,
                 cols = c("#D2B4DE", "#8E44AD", "#6C3483")
) + 
  coord_flip()

plot

cichlid_clusters <- c(
  "0" = '15.4,15.3,15.5,15.2,15.1',
  "1" = '15.4,15.3,15.5,15.2,15.1',
  "2" = '1.1,1.2',
  "3" = 'NA',
  "4" = '2.2',
  "5" = '4.2,4.1, 8.9,8.6, 9.1,9.6,9.3, 15.4,15.2,15.1',
  "6" = '9.6, 9.2, 11.1, 8-9, 8.7,8.10',
  "7" = 'NA',
  "8" = '15.2,15.1, 15.4',
  "9" = 'NA',
  "10" = 'NA',
  "11" = 'NA',
  "12" = 'NA',
  "13" = '11.1, 8.7',
  "14" = '1.3',
  "15" = '8.2,8.1,8.8,8.3, 10.2,10.1, 15.2,15.1',
  "16" = '14',
  "17" = '15.2,15.1',
  "18" = '2.1',
  "19" = '15.5',
  "20" = '6, 15.3',
  "21" = 'NA',
  "22" = '1.2,1.1',
  "23" = '15.2,15.1, 15.5',
  "24" = 'NA',
  "25" = '15.4,15.3',
  "26" = '1.3',
  "27" = '15.4,15.2,15.5',
  "28" = '15.4,15.2,15.5',
  "29" = '15.4,15.2,15.5',
  "30" = '5.2,5.1',
  "31" = '15.2,15.1'
)

#cichlid.spatial <- read.csv('/Volumes/jrhodes/Fish Lab/Experiments/sex change single nuc POA/spatial mapping cichlid/bb_clusters_on_annotated_brain_regions_in_spatial.csv')

cluster.ids <- data.frame(cluster = c(0:31),
                          previous_object_49 = prev_labels,
                          cichlid_cluster  = cichlid_clusters,
                          region = NA
)

#im just gona do it manually

cluster.ids$region <- c(
  "0" = 'Vi, Vl, Vv, SP-u',
  "1" = 'Vi, Vl, Vv, SP-u',
  "2" = 'Radial_Glia',
  "3" = 'NA',
  "4" = 'Oligodendrocyte',
  "5" = 'Vs, Dl-d, Dc-1/2, Dd, Vv',
  "6" = 'Dm-1,Dd, Dp, Dc-4, Dp',
  "7" = 'NA',
  "8" = 'Vv',
  "9" = 'NA',
  "10" = 'NA',
  "11" = 'NA',
  "12" = 'NA',
  "13" = 'Dp',
  "14" = 'Microglia',
  "15" = 'Dm-3, Di-g, Vv',
  "16" = 'tract',
  "17" = 'Vv',
  "18" = 'OPC',
  "19" = 'Vv',
  "20" = 'Dm-2c, Vi',
  "21" = 'NA',
  "22" = 'Radial_Glia',
  "23" = 'Vv',
  "24" = 'NA',
  "25" = 'Vi, Vl, Vv',
  "26" = 'Radial_Glia',
  "27" = 'Vv, Vc',
  "28" = 'Vv, Vc', #
  "29" = 'Vv, Vc',
  "30" = 'OB',
  "31" = 'Vv')

mpoa_markers <- c( #### based on hypomap
                 'cbln2b',
                 'pcsk1',
                 'eomesb',
                 'lhx9',
                 'hmx2',
                 'hmx3a',
                 'stc1',
                 'lmo2',
                 'zic1',
                 'zic5',
                 'zic4',
                 'nfixa',
                 'reln',
                 'nfia',
                 'slc32a1',
                 'sp9',
                 'dlx1a',
                 'ebf1a',
                 'LOC111583095',
                 'otx2b',
                 'mia',
                 'megf11',
                 'moxd1l',
                 'cplx3b',
                 'gpc2',
                 'crhb',
                 'nts',
                 'lhx6a',
                 'zeb2b',
                 'nr2e1',
                 'pbx3b',
                 'galn',
                 'rgs1',
                 'pnoc',
                 'csta2',
                 'isl1',
                 'kcnk9',
                 'serpine2',
                 'ecel1',
                 'amigo2',
                 'bub3',
                 'rasgrp1',
                 'ngfra'
                  )

DotPlot(object = multiome_object, 
                 group.by = "harmony.wnn_res0.4_clusters", 
                 features = mpoa_markers,
                 cols = c("#D2B4DE", "#8E44AD", "#6C3483")
) + 
  coord_flip()



###

pvn_markers <- c('trh', #### based on hypomap
                 'ebf1a',
                 'asic4a',
                 'cbln2b',
                 'zic4',
                 'neurod1',
                 'zic1',
                 'onecut2',
                 'uncx',
                 'cbln1',
                 'zic5',
                 'islr2',
                 'adcyap1a',
                 'zic2a',
                 'sncga',
                 'kcnk9',
                 'foxg1a'
                 )

pvn_2 <- c('shox2',
           'trh',
           'cbln2b',
            'sox14',
                      'lef1',
           'uncx',
           'oxt',
           'tbx19',
           'col11a1',
           'ebf3',
           'crhb',
           'npsr1',
           'pbx3',
           'lhx6a',
           'sparc',
           'npy'


)
           
 DotPlot(object = multiome_object, 
                 group.by = "harmony.wnn_res0.4_clusters", 
                 features = pvn_2,
                 cols = c("#D2B4DE", "#8E44AD", "#6C3483")
) + 
  coord_flip()



cluster.ids$final_annotation <- c(
    "0" = 'POA_1_Mixed',
  "1" = 'POA_2_GABA',
  "2" = 'Radial_Glia',
  "3" = 'gal+_GLUT', #potentially pvn
  "4" = 'Oligodendrocyte',
  "5" = 'Vs, Dl-d, Dc-1/2, Dd, Vv_npy+_cckbra+_GABA',
  "6" = 'Dm-1,Dd, Dp, Dc-4, Dp_cckb+_GLUT',
  "7" = 'tac1+_GLUT',
  "8" = 'Vv_cckb+_GLUT_1',
  "9" = 'tac1+_GLUT_2',
  "10" = 'ccka+_GLUT_1',
  "11" = 'ccka+_GLUT_2',
  "12" = 'ACH+_GLUT_1',
  "13" = 'Dp_GLUT',
  "14" = 'Microglia',
  "15" = 'Dm-3, Di-g, Vv_pdyn+_GLUT',
  "16" = 'tract_GLUT',
  "17" = 'POA_3_GLUT',
  "18" = 'OPC',
  "19" = 'POA_4_Mixed',
  "20" = 'Dm-2c, Vi_npy+_GABA', #probably also hypothamaus
  "21" = 'ACH+_GLUT_2',
  "22" = 'Ependymal',
  "23" = 'POA_5_PVN_GLUT',
  "24" = 'GLUT_1',
  "25" = 'Vi, Vl, Vv_GABA',
  "26" = 'Leukocytes',
  "27" = 'POA_6_PVN_GLUT', #PVN or hypothalamus
  "28" = 'Vv, Vc_ccka+_cckbra+_GABA', ### actually I want to say hypothalamus now 
  "29" = 'Fibroblasts', 
  "30" = 'POA_6_PVN_GABA', #PVN
  "31" = 'Vv_GABA')

cluster_29 <- FindMarkers(multiome_object, 29)
head(cluster_29)  

cluster_27 <- FindMarkers(multiome_object, 27)
head(cluster_27,10)  #hard to get more specific than POA

write.csv(cluster.ids, '/Volumes/jrhodes/Fish Lab/Experiments/sex change single nuc POA/spatial mapping cichlid/mapping multiome regions using previous object and cichlid data.csv')


