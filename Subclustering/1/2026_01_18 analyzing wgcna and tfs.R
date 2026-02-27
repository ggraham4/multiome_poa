# analyzing wgcna in a separate file so I never have to touch that other one again

#remotes::install_version("Seurat", version = "5.0.3")
#remotes::install_version("SeuratObject", version = "5.0.1") 

{
library(SeuratObject)
library(Seurat)
library(hdWGCNA)
library(patchwork)
library(tidyverse)
  clown_go = readRDS("Functions/clown_go2")  

  
  prop_plotter =function(obj, clustering){
    n_cells = obj@meta.data%>%
      group_by(individual, Status)%>%
    summarise(total_cells = n())
    
    n_cells_cluster = obj@meta.data%>%
      group_by(individual, .data[[clustering]])%>%
    summarise(n_cells = n())
    
    joint = n_cells%>%
      left_join(n_cells_cluster, by = 'individual')%>%
      mutate(prop = n_cells / total_cells)
    
    p = ggplot(joint, aes(x = Status, y = prop))+
      geom_boxplot()+
      geom_point()+
      facet_wrap(~sub.cluster, scales = 'free')
    return(p)
    
  }
mecp  = readRDS('Functions/mean_expression_cluster_plot.rds')          

  
  plot_module <-  function(module){
  #ur supposed to use ME not hME
  me_subset=MEs%>%
    group_by(individual, Status)%>%
    summarize(mean_module = mean(.data[[module]]),
              se_module = sd(.data[[module]])/sqrt(n()))
  
  me_subset$Status = factor(me_subset$Status, levels = c('NRM','M',"D",'E','NF','F'))
  
  plot <- ggplot(me_subset, aes(x = Status, y = mean_module, color = Status))+
    geom_boxplot(aes(group = Status, fill = Status), alpha = 0.25, outlier.shape = NA)+
    geom_pointrange(aes(x = Status, y = mean_module,
                        ymin = mean_module - se_module,
                        ymax = mean_module +se_module),
                    position = position_jitterdodge(1))+
    theme_classic()+
    labs(x  ='Status', y = module)+
    theme(axis.text.x = element_text(angle = -45, vjust = 1, hjust=0))
  return(plot)
  
}

plot_module_cluster <-  function(module){
  #ur supposed to use ME not hME
  me_subset=MEs%>%
    group_by(individual, sub.cluster)%>%
    summarize(mean_module = mean(.data[[module]]),
              se_module = sd(.data[[module]])/sqrt(n()))
  
  plot <- ggplot(me_subset, aes(x = sub.cluster, y = mean_module, color = sub.cluster))+
    geom_boxplot(alpha = 0.25, outlier.shape = NA)+
    geom_pointrange(aes(x = sub.cluster, y = mean_module,
                        ymin = mean_module - se_module,
                        ymax = mean_module +se_module),
                    position = position_jitterdodge(1))+
    theme_classic()+
    labs(x  ='sub.cluster', y = module)+
    theme(axis.text.x = element_text(angle = -45, vjust = 1, hjust=0))
  return(plot)
  
}
expression_by_cluster = function(gene){
  dat = rgc@assays$RNA$data[gene,]%>%as.data.frame()
  dat$sub.cluster = rgc$sub.cluster
  dat$individual = rgc$individual    
  dat$Status = rgc$Status
  
  plot_dat = dat%>%
    group_by(sub.cluster, individual, Status)%>%  
    summarize(mean_exp = mean(.))
  
 p =  ggplot(plot_dat, aes(x = sub.cluster, y = mean_exp, color= Status))+
    geom_boxplot()+
    geom_point(  position = position_jitterdodge(.2))
  
 return(p)
    }

}
rgc = readRDS('/Users/ggraham/WGCNA_rgc/rgc_seurat_object.rds')

#modules
modules = GetModules(rgc) %>% subset(module != 'grey')
#hubs
hub_df <- GetHubGenes(rgc, n_hubs = 20)
#MEs
MEs <- GetMEs(rgc, harmonized=FALSE)
MEs$individual = rgc$individual
MEs$Status = rgc$Status
MEs$sub.cluster = rgc$sub.cluster
### ok time to explore my transcription factors of interest ###
deg_tfs = c('skida1',
            'gli1',
            'LOC111563338',# coup
            'LOC111568198', # pbnox2
            'LOC111587547', # dbx
            'znf710b',
            'sox6',
            'foxp2',
            'nkx2.2a'
            )

deg_tfs%in%hub_df$gene_name # no
deg_tfs%in%modules$gene_name # mix

modules[modules$gene_name %in% deg_tfs,] 
# so gli, pbnox, and coup are coexpressed, and sox6 is in its own camp

plot_module('blue')
# elevated in dominants, which is strange because if you look at the actual degs

mecp(rgc, 'gli1', '1')
mecp(rgc, 'LOC111568198', '1') # pbnox
mecp(rgc, 'LOC111563338', '1') # coup
# they are all more or less elevated in males and dominants

# I think its time to see the distribution
FeaturePlot(rgc, c('gli1',
                   'LOC111568198',
                   'LOC111563338'))+
  DimPlot(rgc)
# these all seem like 1_0 biased
rgc$sub.cluster= factor(rgc$sub.cluster, levels = c(paste0('1_', 0:5)))

#im an idiot there is a better way
DotPlot(rgc, deg_tfs, group.by = 'sub.cluster')+
  coord_flip()

mecp(rgc, 'nkx2.2a', '1_2', 'sub.cluster')
mecp(rgc, 'sox6', '1_2', 'sub.cluster')
mecp(rgc, 'foxp2', '1_2', 'sub.cluster')

prop_plotter(rgc, 'sub.cluster')

markers = FindAllMarkers(rgc, only.pos = T)
marks_1_0 = subset(markers, cluster =='1_0')

# perplexity suggests 1_0 is glutamatergic neurogenic progenitors
"The co-expression of foxg1a, pax6b, and zic3 establishes unambiguous dorsal telencephalic (cortical) identity"
# interestingly, this population is decreasd in females but more or less the same in dominants and males

# now what about 1_2
marks_1_2 = subset(markers, cluster =='1_2')
"The co-expression of nkx2.1, shha, and six6a immediately establishes ventral telencephalic and/or hypothalamic origin. "
"[nkx2.1] specifies cell type identity in cells that are already post-mitotic or transitioning from proliferation to differentiation."
# so we can call 1_2 maturing neurons
"asic2 marks synaptic maturation and functional integration. ASIC2 (acid-sensing ion channel 2) localizes to dendritic spines, physically associates with ASIC1a and the scaffolding protein PSD-95, and is critical for targeting ASIC1a to synapses"
"Scenario 1: Developing GABAergic Interneurons (Most Likely)
Scenario 2: Hypothalamic Neuroendocrine Neurons
Scenario 3: Mature Functional Astrocyte Cluster (least likely)
"
"GPT thinks this is more likely astrocytic, based on my previously constricted lineage I agree with perplexity more
"


# and 1_4
marks_1_4 = subset(markers, cluster =='1_4')
"neurogenic midbrain domain actively specifying dopaminergic and other mesencephalic neurons"
# very prominent wnt signaling
"Dopaminergic Neurons (Primary Output)"
"GPT Importantly, these cells are not primarily producing neurons themselves.,
Act as a signaling organizer, 
more like roof plate"

marks_1_5 = subset(markers, cluster =='1_5')
"Gene Expression Profile of Fourth Cluster: Mature/Maturing Dopaminergic Neurons
The Dopaminergic Identity Triad: Pitx3, Lmx1a, and Foxa1
"
# so the idea here is 1_4 and 1_5 are related 
"Strong dopaminergic likelihood, Post-mitotic dopaminergic neurons, or Late DA progenitors"

DimPlot(rgc, label = T)
# I am not convinced

# the 1_1 - 1_3 - 1_2 axis I think is real

quiescence =read.csv("~/Desktop/UIUC/Zack Code/Radial Glia Quiescence Score/Clownfish Genes.csv")

rgc = AddModuleScore(rgc, features = list(
  'cycling' =c( quiescence$gene[quiescence$rg_state == 'cycling']),
  'neuronal differentiation' =c(quiescence$gene[quiescence$rg_state == 'neuronal differentiation']),
  'quiescent' =c(quiescence$gene[quiescence$rg_state == 'quiescent'])
), name = 'quiescence')

DotPlot(rgc, c(paste0('quiescence', 1:3)))+
  coord_flip()

DotPlot(rgc, c('setd1a',
               'dazap2',
               'fbxl5',
               'ccnd1',
               'fgfbp3',
               's100b', # quiescent begin ----
               'LOC111575074', # fabp7
               'LOC111570141' # glula
               ))+
  coord_flip()


DotPlot(rgc, c('elavl3',
               'sox11a',
               'elavl2'
               ))+
  coord_flip()


###
hdWGCNA::ModuleRadarPlot(rgc)

clown_go(modules$gene_name[modules$color=='blue'])%>%dotplot() # translation 
clown_go(modules$gene_name[modules$color=='yellow'])%>%dotplot()  # neurogenesis
clown_go(modules$gene_name[modules$color=='red'])%>%dotplot() # probably also neurogenesis
clown_go(modules$gene_name[modules$color=='green'])%>%dotplot() # this is why I think its a tanycyte

