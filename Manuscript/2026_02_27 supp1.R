
#supp fig 1
#a - Dotplot 
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
  library(emmeans)
  library(ggsignif)
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

}


multiome_object <- readRDS("~/Desktop/nemo.orig_harmony.integration_all_testd_clusters.rds") 
Idents(multiome_object) = 'res0.8_50nn_40PC_45LSI'
`%notin%` <- Negate(`%in%`)
obj = multiome_object

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

#ggsave(plot = marker_gene_plot,
 #      file = "marker_gene_plot.svg",
  #     device = "svg",
   #    units = "in",
    #   width = 6,
     #  height = 4.5,
      # path = "Manuscript/Plots/Supp 1")

### proportions####
## b - cell proportions####

status_to_phase = c(
  "D"='I',
  'M' = 'M',
  'F' ='F',
  'NF' ='NF',
  'NRM' = 'NRM',
  'E' = 'LI'
)

proportion_plotter = function(data_frame, cluster_id, signif_data){
  dag = subset(data_frame, res0.8_50nn_40PC_45LSI == cluster_id &
                 Status != 'NRM')
  
  deg_data = subset(signif_data, cluster == cluster_id)  # fixed scoping bug
  d_f_p.value = deg_data$d_f_p.value
  d_m_p.value = deg_data$d_m_p.value
  f_m_p.value = deg_data$f_m_p.value
  
  d_f_annotation = p_value_annotation(d_f_p.value)
  d_m_annotation = p_value_annotation(d_m_p.value)
  f_m_annotation = p_value_annotation(f_m_p.value)

  dag$Phase = status_to_phase[dag$Status]
  dag$Phase = factor(dag$Phase, levels = c('M', 'I', 'LI', 'NF', "F"))
  
  dag$prop = dag$ncells / dag$total_cells  # create prop column explicitly
  
  r = ggplot(dag, aes(x = Phase, y = prop))+
    geom_boxplot(data = subset(dag, Phase != "NF"),
                 aes(x = Phase, y = prop,
                     fill = Phase), outlier.shape = NA)+
    geom_point(size = 0.5)+
    theme_classic()+
    scale_x_discrete(drop = FALSE)+  
    scale_y_continuous(labels = scales::percent)+
    labs(x = 'Phase', y = 'Proportion of cells')+
    geom_signif(xmin = c(1),
                xmax = c(1.9),
                y_position = c(max(dag$prop * 1.1)),
                annotation = c(d_m_annotation),
                color = "black",
                tip_length = c(0, 0), textsize = 3)+
    geom_signif(xmin = c(1),
                xmax = c(5),
                y_position = c(max(dag$prop * 1.3)),
                annotation = c(f_m_annotation),
                color = "black",
                tip_length = c(0, 0), textsize = 3)+
    geom_signif(xmin = c(2.1),
                xmax = c(5),
                y_position = c(max(dag$prop * 1.1)),
                annotation = c(d_f_annotation),
                color = "black",
                tip_length = c(0, 0), textsize = 3)+
    theme(legend.position = 'none', axis.title.y = element_blank())
  
  return(r)
}
total_cells = obj@meta.data%>%
  group_by(individual, Status)%>%
  summarize(total_cells = n())

n_cells=obj@meta.data%>%
  group_by(individual, res0.8_50nn_40PC_45LSI)%>%
  summarize(ncells = n())

joint = total_cells%>%
  left_join(n_cells, by = 'individual')

joint$Status = as.character(joint$Status)

dat = data.frame()
for(i in 0:26){
  sub= subset(joint, res0.8_50nn_40PC_45LSI ==i & Status %in% c('M','D','F'))
  mat = cbind(sub$ncells, sub$total_cells-sub$ncells)
  
  mod = glm(mat ~ Status, data = sub, family = 'binomial')
  p = anova(mod, 'III')
  p.value = p$`Pr(>Chi)`[2]
  
  pair = pairs(emmeans(mod, 'Status'), adjust = 'none')%>%as.data.frame()
  
  newd = data.frame(cluster = i,
                    p = p.value,
                    d_m_p.value = pair$p.value[pair$contrast == 'D - M'],
                    f_m_p.value = pair$p.value[pair$contrast == 'F - M'],
                    d_f_p.value = pair$p.value[pair$contrast == 'D - F'],
                    d_m_estimate = pair$estimate[pair$contrast == 'D - M'],
                    f_m_estimate = pair$estimate[pair$contrast == 'F - M'],
                    d_f_estimate = pair$estimate[pair$contrast == 'D - F']
)
  dat = rbind(newd, dat)

  
}

dat$av_q.value = p.adjust(dat$p, 'fdr', 27)
dat$sig = ifelse(dat$av_q.value < 0.05, '*', NA)
dat$singular = F

textsize = 3

proportion_plotter(joint_neuron, 0, dat)

# systematically write all plots 
differentially_present_cells = dat$cluster[!is.na(dat$sig )]

for(i in differentially_present_cells){
  p = proportion_plotter(joint, i, dat)+
      labs(title = i)+
    theme(axis.title.x = element_blank())
  
    ggsave(plot = p,
       file = paste0('cluster_', i, '_cell.svg'),
       device = "svg",
       units = "in",
       width = 1.35,
       height = 1.45,
       path = "Manuscript/Plots/Supp 1/cell_proportions")
  
}

write.csv(dat, '/Users/ggraham/Desktop/multiome_poa/Manuscript/Supplementary Tables/proportion_of_cells_stats.csv')
write.csv(joint, '/Users/ggraham/Desktop/multiome_poa/Manuscript/Supplementary Tables/proportion_of_cells_raw.csv')


## c - neuron proportions ####

neurons_only <- subset(obj, 
                     #oligos
                     res0.8_50nn_40PC_45LSI!=2&
                     #microglia
                    res0.8_50nn_40PC_45LSI!=11&
                    #opcs
                    res0.8_50nn_40PC_45LSI!=13&
                    #dividing glia
                    res0.8_50nn_40PC_45LSI!=26&
                    #leuko
                    res0.8_50nn_40PC_45LSI!=20&
                    #ependymal
                    res0.8_50nn_40PC_45LSI!=15
                    &  res0.8_50nn_40PC_45LSI!=1
                    )

total_cells_neuron = neurons_only@meta.data%>%
  group_by(individual, Status)%>%
  summarize(total_cells = n())

n_cells_neuron=neurons_only@meta.data%>%
  group_by(individual, res0.8_50nn_40PC_45LSI)%>%
  summarize(ncells = n())

joint_neuron = total_cells_neuron%>%
  left_join(n_cells_neuron, by = 'individual')

joint_neuron$Status = as.character(joint_neuron$Status)

dat_neuron = data.frame()
for(i in unique(joint_neuron$res0.8_50nn_40PC_45LSI)){
  sub= subset(joint_neuron, res0.8_50nn_40PC_45LSI ==i & Status %in% c('M','D','F'))
  mat = cbind(sub$ncells, sub$total_cells-sub$ncells)
  
  mod = glm(mat ~ Status, data = sub, family = 'binomial')
  p = anova(mod, 'III')
  p.value = p$`Pr(>Chi)`[2]
  
  pair = pairs(emmeans(mod, 'Status'), adjust = 'none')%>%as.data.frame()
  
  newd = data.frame(cluster = i,
                    p = p.value,
                    d_m_p.value = pair$p.value[pair$contrast == 'D - M'],
                    f_m_p.value = pair$p.value[pair$contrast == 'F - M'],
                    d_f_p.value = pair$p.value[pair$contrast == 'D - F'],
                    d_m_estimate = pair$estimate[pair$contrast == 'D - M'],
                    f_m_estimate = pair$estimate[pair$contrast == 'F - M'],
                    d_f_estimate = pair$estimate[pair$contrast == 'D - F']
)
  dat_neuron = rbind(newd, dat_neuron)

  
}

dat_neuron$av_q.value = p.adjust(dat_neuron$p, 'fdr',nrow(dat_neuron))
dat_neuron$sig = ifelse(dat_neuron$av_q.value < 0.05, '*', NA)

differentially_present_neurons = dat_neuron$cluster[!is.na(dat_neuron$sig )]

for(i in differentially_present_neurons){
  p = proportion_plotter(joint_neuron, i, dat)+
    labs(title = i)+
    theme(axis.title.x = element_blank())
  
    ggsave(plot = p,
       file = paste0('cluster_', i, '_neuron.svg'),
       device = "svg",
       units = "in",
       width = 1.35,
       height = 1.45,
       path = "Manuscript/Plots/Supp 1/neuron_proportions")
  
}

write.csv(dat_neuron, '/Users/ggraham/Desktop/multiome_poa/Manuscript/Supplementary Tables/proportion_of_neurons_stats.csv')
write.csv(joint_neuron, '/Users/ggraham/Desktop/multiome_poa/Manuscript/Supplementary Tables/proportion_of_neurons_raw.csv')




