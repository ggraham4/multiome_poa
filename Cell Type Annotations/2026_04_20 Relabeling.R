# renaming and relabeling celltypes

library(Seurat)
library(ggplot2)

obj = readRDS("C:/Users/Gabe/Desktop/nemo.orig_harmony.integration_all_testd_clusters.rds")
DotPlot(obj, 
        c('eomesa'),
        group.by = 'res0.8_50nn_40PC_45LSI')

Marks_26 =FindMarkers(obj, 26, group.by = 'res0.8_50nn_40PC_45LSI')

new_labs = list(
  '0' = 'POA_1', # mixed
  '5' = 'POA_3', #gab 
  '25' = 'POA_2', #gaba
  '4' = 'POA_4', #glut
  '7' = 'POA_5', #gab
  '8' = 'POA_6', #glut
  '17' = 'POA_7', #glut
  '18' = 'POA_8', #glut
  '1'='RGC',
  '2'='OL',
  '3' ='VD/VS',
  '6' ='nPO',
  '9' = 'IMN_1',
  '10' = 'ACh_Kiss',
  '11' ='MG',
  '12' = 'IMN_2',
  '13' = 'OPC',
  '14'= 'IMN_3',
  '15' = 'EC',
  '16' = 'ACh_Tac',
  '19' = 'Dl',
  '20' = 'Leuko',
  '21' = 'Pallidum',
  '23' = 'IMN_5',
  '22' ='IMN_4',
  '24' = 'Mixed',
  '26' = 'DG'
)



renamed = unname(unlist(new_labs[as.character(obj$res0.8_50nn_40PC_45LSI)]))

obj$labs = renamed
obj$labs = factor(obj$labs, levels =
                    c(
                      'nPO',
                      paste0('POA_', 1:8),
                      'ACh_Kiss',
                      'ACh_Tac',
                      'DA',
                      'Mixed',
                      'Dl',
                      'VD/VS',
                      'Pallidum',
                      paste0("IMN_",1:5),
                      'RGC',
                      'EC',
                      'DG',
                      'OL',
                      'OPC',
                      'MG',
                      'Leuko'
                    ))

DotPlot(obj, group.by = 'labs',
        features = c(
          'elavl3',
          'gfap',
          'LOC111577263',
          'crocc2',
          'top2a',
          'mbpa',
          'cspg4',
          'ptprc',
          'p2ry12',
          'LOC111588076', #gad1
          'gad2',
          'slc17a6a',
          'slc17a6b',
          'slc17a7a',
          'slc18a3b',
          'th2',
          'hmx2',
          'hmx3a',
          'hmx3b',
          'dlx2a',
          'dlx5a',
          'nes',
          'sox2',
          'sox4b',
          'sox11a',
          'mex3a',
          'pax6b',
          'eomesa',
          'cux2b',
          'zic2a',
          'tac1',
          'tacr1a',
          'kiss1',
          'kiss1rb',
          'sst1.1'
          
        ),
        dot.min = 0.1
        )+
  coord_flip()+
  theme(axis.text.x = element_text(angle =-90,
                                   vjust =.5))+
  scale_size(range = c(1,4))



  