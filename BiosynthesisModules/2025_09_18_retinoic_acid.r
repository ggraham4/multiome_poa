
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

plotSeuratModule = function(obj , module, cluster, clustering){
    me_subset=obj@meta.data%>%
    group_by(individual, Status, .data[[clustering]])%>%
    summarize(mean_module = mean(.data[[module]]),
              se_module = sd(.data[[module]])/sqrt(n()))%>%
      filter(.data[[clustering]] == cluster)    
    
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


#### make a retinoic acid module
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

retinoic_go_genes = unique(
  biomart_basic$entrezgene_accession[biomart_basic$name_1006 ==
                                       grepl('retinoic acid', biomart_basic$name_1006) |
                                         grepl('retinoic acid', biomart_basic$entrezgene_description) 

                                     ]
  )
retinoic_go_genes

rgc = AddModuleScore(rgc, list(retinoic_go_genes),name = 'RetGo')

DotPlot(rgc, 'RetGo1')

plotSeuratModule(rgc, 'RetGo1', '1_1', 'sub.cluster')

data_1_1 = subset(rgc@meta.data, sub.cluster == '1_1')
retgo_model = lmer(RetGo1 ~ Status+(1|individual), data =data_1_1 )
car::Anova(retgo_model, 3) 




curated_retinoic_acid = c('LOC111582092',
                'rdh8a',
                'rdh12',
                'rdh20',
                'rdh1',
                'LOC111572826',                
                'LOC111587579',
                'LOC111587581',
                'LOC111579690',
                'si:dkey-23o4.6',
                'LOC111575778',
                'LOC111587582',
                'LOC111575523',
                'LOC111566272',
                'LOC111569321',
                'rdh14b',
                'rdh5',
                'rdh10a',
                'LOC111570850')

obj <- AddModuleScore(obj, list(retinoic_go_genes),
                      name = 'RetGo')

DotPlot(obj, 'RetGo1')

plotSeuratModule(obj, 'RetGo1', '4', 'final_clusters')
plotSeuratModule(obj, 'RetGo1', '13', 'final_clusters')
plotSeuratModule(obj, 'RetGo1', '7', 'final_clusters')
plotSeuratModule(obj, 'RetGo1', '22', 'final_clusters')
plotSeuratModule(obj, 'RetGo1', '8', 'final_clusters')
plotSeuratModule(obj, 'RetGo1', '5', 'final_clusters')
plotSeuratModule(obj, 'RetGo1', '6', 'final_clusters')
plotSeuratModule(obj, 'RetGo1', '19', 'final_clusters')
plotSeuratModule(obj, 'RetGo1', '12', 'final_clusters')
plotSeuratModule(obj, 'RetGo1', '0', 'final_clusters')
plotSeuratModule(obj, 'RetGo1', '24', 'final_clusters')
plotSeuratModule(obj, 'RetGo1', '14', 'final_clusters')
plotSeuratModule(obj, 'RetGo1', '13', 'final_clusters')

data_13 = subset(obj@meta.data, final_clusters == '13' & Status %in%c('M',"D",'F'))
retgo_model3 = lmer(RetGo1 ~ Status+(1|individual), data =data_13 )
car::Anova(retgo_model3, 3) 
pairs(emmeans(retgo_model3,'Status'), adjust ='none') # dominants significantly higher than males

for(i in 0:26){
  print(plotSeuratModule(obj, 'RetGo1', i, 'final_clusters')+
          labs(title = i))
}
# interesting clusters 
#> 26
#> 20 maybe
#> 17 maybe
#> 16 
#> 15
#> 2!

plotSeuratModule(obj, 'RetGo1', '2', 'final_clusters')

data_2 = subset(obj@meta.data, final_clusters == '2' & Status %in%c('M',"D",'F'))
retgo_model3 = lmer(RetGo1 ~ Status+(1|individual), data =data_2 )
car::Anova(retgo_model3, 3)  #how???


###curated list ####
obj <- AddModuleScore(obj, list(curated_retinoic_acid),
                      name = 'Retinoic')
DotPlot(obj, 'Retinoic1')



retinoic_data = obj@meta.data%>%
  group_by(individual, Status, final_clusters)%>%
  summarize(mean_retinoic = mean(Retinoic1),
            se_retinoic = sd(Retinoic1)/sqrt(n()))

retinoic_data$Status <- factor(retinoic_data$Status, levels= c('NRM','M','D','E','NF','F'))

ggplot(subset(retinoic_data, final_clusters ==8), aes(x = Status, y = mean_retinoic))+
  geom_boxplot()+
    geom_point()

ggplot(subset(retinoic_data, final_clusters ==22), aes(x = Status, y = mean_retinoic))+
  geom_boxplot()+
    geom_point()

ggplot(subset(retinoic_data, final_clusters ==24), aes(x = Status, y = mean_retinoic))+
  geom_boxplot()+
    geom_point()

ggplot(subset(retinoic_data, final_clusters ==23), aes(x = Status, y = mean_retinoic))+
  geom_boxplot()+
    geom_point()

dme_data = data.frame()
for(i in 0:26){
  dat = subset(obj@meta.data, final_clusters ==i & Status %in% c('M','D',"F"))
  mod = lmer(Retinoic1~Status+(1|individual), data = dat)
  av = car::Anova(mod, 3)
  p = av$`Pr(>Chisq)`[2]
  dme_data = rbind(dme_data, data.frame(cluster = i, 
                                        p = p,
                                        singular = isSingular(model)))
  
}

ggplot(subset(retinoic_data, final_clusters ==14), aes(x = Status, y = mean_retinoic))+
  geom_boxplot()+
    geom_point()

