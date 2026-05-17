library(Seurat)
library(patchwork)
library(tidyverse)
library(ggplot2)
library(Polychrome)
library(emmeans)
library(ggsignif)
clown_go = readRDS("Functions/clown_go2")  
library(clusterProfiler)
library(CytoTRACE)
library(lme4)

obj  = readRDS("~/Desktop/optimal_clustering_rna_only.rds")

sub_6 = FindSubCluster(obj, 6, graph.name='harmony.wsnn')
Idents(sub_6) <- 'sub.cluster'
sub_6 = subset(sub_6, final_clusters ==6)
sub_6$sub.cluster = factor(sub_6$sub.cluster, levels = c(paste0('6_',0:3)))
sub_6$Status = factor(sub_6$Status, levels = c('NRM','M',"D",'E','NF','F'))

### proportion of 6 #####
total_cells = sub_6@meta.data%>%
  group_by(Status, individual)%>%
  subset(Status !='NRM')%>%
  summarize(total_cells = n())

######degs
degs = read.csv('/Users/ggraham/Desktop/multiome_poa/DEG Outputs/FINAL degs classified w singular.csv')
degs_6 = subset(degs, cluster == 6)

p_value_annotation <- function(p) {
  if (p < 0.001) {
    return("***")
  } else if (p < 0.01) {
    return("**")
  } else if (p < 0.05) {
    return("*")
  } else {
    return(paste0("p = ", round(p, 3)))
  }
}

obj@meta.data$Phase = obj@meta.data$Status
obj@meta.data$Phase = ifelse(obj@meta.data$Phase == 'D', 'I', obj@meta.data$Phase)
obj@meta.data$Phase = ifelse(obj@meta.data$Phase == 'E', 'LI', obj@meta.data$Phase)
obj@meta.data$Phase = factor(obj@meta.data$Phase, levels = c('NRM','M','I','LI','NF','F'))

make_summary = function(dat, numerator_col, denominator_col) {
  dat %>%
    filter(Phase != "NF") %>%
    group_by(Phase) %>%
    summarise(
      mean_prop = mean(.data[[numerator_col]] / .data[[denominator_col]], na.rm = TRUE),
      se = sd(.data[[numerator_col]] / .data[[denominator_col]], na.rm = TRUE) / sqrt(n()),
      .groups = 'drop'
    )
}

### cyto by degs ####
degs =read.csv('/Users/ggraham/Desktop/multiome_poa/DEG Outputs/FINAL degs classified w singular.csv')

degs_plasticity = c(
  'LOC111588913', # nxph1
  'cntn4',
  'LOC111567620', #cntnapl4
  'pcdh10b',
  'sdc2',
  'LOC111585095', # dsel
  'bcan'
)

degs_neurotransmission = c(
  'gabra3',
  'nmbr',
  'npy7r',
  'pgr',
  'LOC111577100', # cd87
  'grm5a',
  'LOC111568069', #pgrl
  'tacr3a',
  'drd3'
)


sub_6_plasticity = sub_6@assays$RNA$data[degs_plasticity,]
sub_6_nt = sub_6@assays$RNA$data[degs_neurotransmission,]

plasticity_positive_cells = colSums(sub_6_plasticity)>0
nt_positive_cells = colSums(sub_6_nt)>0


sub_6$plasticity = ifelse(plasticity_positive_cells, T, F)
sub_6$nt = ifelse(nt_positive_cells, T, F)

FeaturePlot(sub_6, 'plasticity')
FeaturePlot(sub_6, 'nt')

positives = sub_6@meta.data%>%
  group_by(individual, Status)%>%
  summarize(n_plas = sum(plasticity==T),
            n_nt = sum(nt ==T))%>%
  right_join(total_cells, by = 'individual')%>%
  mutate(prop_nt = n_nt/total_cells,
         prop_plas = n_plas/total_cells)

ggplot(positives, aes(x = Status.x, y = prop_plas))+
  geom_boxplot()+
  geom_point() # only plasticity is true

ggplot(positives, aes(x = Status.x, y = prop_nt))+
  geom_boxplot()+
  geom_point()

sub_6$tacr3a = ifelse(sub_6@assays$RNA$data['tacr3a',]>0, T, F)
sub_6$cckb = ifelse(sub_6@assays$RNA$data['cckb',]>0, T, F)

