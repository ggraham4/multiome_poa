### DNA DAMAGE


{
  library(patchwork)
  library(ggplot2)
  library(Seurat)
  library(dplyr)
  library(tidyverse)
  library(CytoTRACE)
  library(BiocManager)
  library(lme4)
  library(enrichR)
  `%notin%` = Negate(`%in%`)
  library(car)
  library(emmeans)
  library(glmGamPoi)
  library(scran)
  library(parallel)
    library(parallel)

  library(Seurat)
  library(tidyr)
  library(lme4)
  library(dplyr)
  library(MASS)
  library(Signac)
  library('glmGamPoi')
  library(scran)
  library(emmeans)
  library(openxlsx)
  library(ggplot2)
  library(stringr)
  library(forcats)
  library(clusterProfiler)
  library(biomaRt)
  library(hdWGCNA)
  
  clown_go = readRDS("Functions/clown_go2")  
  mecd = readRDS("Functions/mean_expression_cluster_data.rds")
    mecp = readRDS("Functions/mean_expression_cluster_plot.rds")
    
geneNamer = function(gene){
  names = read.csv('Reference/genes updated.csv')
  
  name = names$NIH_description[names$NIH_accession==gene][1]
  return(name)
}

plotSeuratModule = function(module, cluster){
    me_subset=rgc@meta.data%>%
    group_by(individual, Status, sub.cluster)%>%
    summarize(mean_module = mean(.data[[module]]),
              se_module = sd(.data[[module]])/sqrt(n()))%>%
      subset(sub.cluster == cluster)
    
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
  `%notin%` <- Negate(`%in%`)
  
  #define functions

}
#read in data
obj = readRDS("~/Desktop/optimal_clustering_rna_only.rds")

Idents(obj) = 'final_clusters'
rgc <- FindSubCluster(obj, 1, 'harmony.wsnn')
rgc = subset(rgc, final_clusters == 1)

Idents(rgc) <- 'sub.cluster'

DimPlot(rgc, label = T, reduction = 'harmony_wnn.umap')

rgc$Status = factor(rgc$Status, levels = c('NRM', 'M','D','E','NF','F'))


#### make a DNA damage genes list
real_degs = read.csv( 'Subclustering/degs_1_defined_09_12_2025.csv')

real_degs$name =c(sapply(X=real_degs$gene, FUN = geneNamer))
sparse =(real_degs[,c('gene', 'cluster', 'short_label', 'name')])

#rip the GO term
ensembl <- useEnsembl(biomart = "genes", 
                      dataset = "aocellaris_gene_ensembl")


biomart_basic <-getBM(
    mart = ensembl, #working mart 
    attributes = c("entrezgene_accession",
                   'entrezgene_description',
                   'go_id',
                   'name_1006',
                   'namespace_1003'))

dna_damage_genes = unique(
  biomart_basic$entrezgene_accession[biomart_basic$name_1006 ==
                                       grepl('break repair', biomart_basic$name_1006) |
                                       grepl('DNA repair', biomart_basic$name_1006) |
                                       grepl('DNA damage', biomart_basic$name_1006) |
                                       grepl('DNA damage', biomart_basic$entrezgene_description)
                                     ]
  )
dna_damage_genes

rgc = AddModuleScore(rgc, list(dna_damage_genes),name = 'Damage')

DotPlot(rgc, 'Damage1')

plotSeuratModule('Damage1', '1_0')
plotSeuratModule('Damage1', '1_1')
plotSeuratModule('Damage1', '1_2') #woooah
plotSeuratModule('Damage1', '1_3')
plotSeuratModule('Damage1', '1_4')
plotSeuratModule('Damage1', '1_5')

data_1_1 = subset(rgc@meta.data, sub.cluster == '1_1')
damage_model = lmer(Damage1 ~ Status+(1|individual), data =data_1_1 )
car::Anova(damage_model, 3)

data_1_2 = subset(rgc@meta.data, sub.cluster == '1_2')
damage_model2 = lmer(Damage1 ~ Status+(1|individual), data =data_1_2 )
car::Anova(damage_model2, 3) # big money but what does it mean


# what does the rest of the object look like
obj = AddModuleScore(obj, list(dna_damage_genes),name = 'Damage')
DotPlot(obj, 'Damage1')

DotPlot(subset(obj, final_clusters!=26), 'Damage1')

rgc$LOC111577260 = rgc@assays$RNA$data['LOC111577260',]
plotSeuratModule('LOC111577260', '1_1')

### is this due to increased txn in males and females in 1_2 ----
rgc$txn = colSums(rgc@assays$RNA$data)
plotSeuratModule('txn', '1_2') 

rgc$txn_raw = colSums(rgc@assays$RNA$counts)
plotSeuratModule('txn_raw', '1_2') 
mod_txn_raw =  lmer(txn_raw ~ Status+(1|individual), data =subset(rgc@meta.data, sub.cluster == '1_2' & Status %in%c('M','D',
                                                                                                             'F')) )
car::Anova(mod_txn_raw, 3)


### is this due to increased DNA synthesis in males and females ----
dna_synthesis_genes = unique(
  biomart_basic$entrezgene_accession[biomart_basic$name_1006 ==
                                       grepl('DNA replication', biomart_basic$name_1006) |
                                       grepl('DNA polymerase', biomart_basic$name_1006) |
                                       grepl('DNA clamp', biomart_basic$name_1006) |
                                       grepl('DNA primase', biomart_basic$entrezgene_description)
                                     ]
  )
dna_synthesis_genes

rgc = AddModuleScore(rgc, list(dna_synthesis_genes), name = 'Replication')
DotPlot(rgc, 'Replication1')

plotSeuratModule('Replication1', '1_2') 
# I dont think so 

mod_rep =  lmer(Replication1 ~ Status+(1|individual), data =subset(rgc@meta.data, sub.cluster == '1_2' & Status %in%c('M','D', 'F')) )
car::Anova(mod_rep, 3)

plotSeuratModule('Replication1', '1_1') 

RNA_synthesis_genes = unique(
  biomart_basic$entrezgene_accession[biomart_basic$name_1006 ==
                                       grepl('transcription', biomart_basic$name_1006) |
                                       grepl('RNA polymerase', biomart_basic$name_1006) 
                                     ]
  )
rgc = AddModuleScore(rgc, list(RNA_synthesis_genes), name = 'transcription')
DotPlot(rgc, 'transcription1')
plotSeuratModule('transcription1', '1_2') 
plotSeuratModule('transcription1', '1_0') # that could be something
plotSeuratModule('transcription1', '1_4') 





