library(Seurat)
library(ggplot2)
library(tidyverse)

obj = readRDS("~/Desktop/optimal_clustering_rna_only.rds")

marks =FindAllMarkers(obj, only.pos = T)%>%
  subset(p_val_adj<0.05)

genes = c(
  # cell types
  'elavl3',
  
  
  'slc17a6b',
  'slc17a7a', #glut
  'LOC111588076', #gad1b
  'gad2',#gaba 
  
  #0
  'soul5l', 
  'glra1', #glycine receptor subunit
  
  #1
  'LOC111577263', #arob
  'gfap', 
  
  #2
  'mbpa',
  
  #3
  'kcnt1a',
  'scn5lab',
  
  #4
  'tac1',
  'tshz1',
  
  #5
  'lhx9',
  'sema3h',
  
  #POA
  'prkd1',
  'hmx2',
  'hmx3a',
  'hmx3b',
  'avp',
  'oxt',
  'cckb', 
  'esr2b',
  'ar',
  
  #7
  'bhlhe22',
  'znf536',
  #8
  'barhl1b',
  'prdm8b',
  
  #9
  'mex3a',
  
  #10
  'nr4a2a',
  'prkcq',
  
  #11
  'ptprc',
  'p2ry12',
  
  #12
  'susd5',
  'tbr1b',
  
  #13
  'cspg4',
  
  # vsst
  'npy',
  'sst1.1',
  
  #15
  'crocc2',
  
  #16
  'rap1gapl',
  'nptx2b',
  
  #17
  'nr5a2',
  'calb2a',
  
  #18
  'nmbr',
  'LOC111577964',
  
  #19
  'gata2b',
  'gata3',
  
  #21
  'lhx6a',
  'LOC111581922', #lhx8
  
  #22
  'ghrhrb',
  'esrrb',
  
  #23
  'pde9ac',
  'tmc2b',
  
  #24
  'fat2',
  'arhgap46a',
  
  #25
  'LOC111585833',
  'cart2',
  
  #26
  'cenpp',
  'top2a'
)

p=DotPlot(obj, genes, 
          dot.min = 0.1,
          dot.scale = 3)+
  coord_flip()+
  scale_x_discrete(labels = c(
  # cell types
  'elavl3',
  
  
  'slc17a6b',
  'slc17a7a', #glut
  'gad1b-like', #gad1b
  'gad2',#gaba 
  
  #0
  'soul5l', 
  'glra1', #glycine receptor subunit
  
  #1
  'cyp19a1b', #arob
  'gfap', 
  
  #2
  'mbpa',
  
  #3
  'kcnt1a',
  'scn5lab',
  
  #4
  'tac1',
  'tshz1',
  
  #5
  'lhx9',
  'sema3h',
  
  #POA
  'prkd1',
  'hmx2',
  'hmx3a',
  'hmx3b',
  'avp',
  'oxt',
  'cckb', 
  'esr2b',
  'ar',
  
  #7
  'bhlhe22',
  'znf536',
  #8
  'barhl1b',
  'prdm8b',
  
  #9
  'mex3a',
  
  #10
  'nr4a2a',
  'prkcq',
  
  #11
  'ptprc',
  'p2ry12',
  
  #12
  'susd5',
  'tbr1b',
  
  #13
  'cspg4',
  
  # vsst
  'npy',
  'sst1.1',
  
  #15
  'crocc2',
  
  #16
  'rap1gapl',
  'nptx2b',
  
  #17
  'nr5a2',
  'calb2a',
  
  #18
  'nmbr',
  'cb1b-like', #cb1b-like
  
  #19
  'gata2b',
  'gata3',
  
  #21
  'lhx6a',
  'lhx8-like', #lhx8
  
  #22
  'ghrhrb',
  'esrrb',
  
  #23
  'pde9ac',
  'tmc2b',
  
  #24
  'fat2',
  'arhgap46a',
  
  #25
  'dmbx1-like',
  'cart2',
  
  #26
  'cenpp',
  'top2a'
))

p     

ggsave(plot = p,
       file = 'marker dotplot.svg',
       device = "svg",
       units = "in",
       width = 5.5,
       height = 6.5,
       path = '/Users/ggraham/Desktop/multiome_poa/Manuscript/Plots/Manuscript v1.2.3')


library(scCustomize)

p2=Clustered_DotPlot(obj, 
                     rev(genes), 
                     cluster_feature = F,
                     colors_use_exp = c('lightgrey','blue'),
                     plot_km_elbow=F,
                     show_annotation_name=F,
                     )
p2    

library(ComplexHeatmap)

# Open an SVG device, draw, then close
svg(filename = '/Users/ggraham/Desktop/multiome_poa/Manuscript/Plots/Manuscript v1.2.3/marker dotplot clustered.svg',
    width = 6.5,
    height = 6.5)
draw(p2)
dev.off()

#ggsave(plot = p2,
 #      file = 'marker dotplot clustered.svg',
  ##     device = "svg",
    #   units = "in",
     #  width = 5.5,
      # height = 6.5,
       #path = '/Users/ggraham/Desktop/multiome_poa/Manuscript/Plots/Manuscript v1.2.3')



