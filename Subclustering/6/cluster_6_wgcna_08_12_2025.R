#> 6 subcluster WGCNA 08_11_2025 GJG
#> 
#> Having seen different GEX DEG patterns in cluster 6 hormone responsive genes
#> pgr vs pgr-like, I want to see if they belong to different modules, and
#> if the sexes differ in those modules, if there are interesting TFs in those
#> modules. I want to see which modules seem to occur prior to gonadal sex change
#> and which occur after
#> 

#libraries
{
  library(ggridges)
  library(ggplot2)
  library(Seurat)
  library(dplyr)
  library(tidyverse)
  library(CytoTRACE)
  library(BiocManager)
  library(lme4)
  library(enrichR)
  library(car)
  library(emmeans)
  library(hdWGCNA)
  library(WGCNA)
  library(patchwork)
  
  `%notin%` = Negate(`%in%`)
  clown_go = readRDS("Functions/clown_go2")  
  mecd = readRDS("Functions/mean_expression_cluster_data.rds")
}
#define functions
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
                    position = position_jitterdodge(3))+
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
#read in object
obj = readRDS("C:/Users/Gabe/Desktop/nemo.orig_harmony.integration_all_testd_clusters.rds")

#subcluster
obj_subclustered_6 <- FindSubCluster(obj, 6, 'harmony.wsnn')
#reduce size
obj_6_only <- subset(obj_subclustered_6, sub.cluster %in% c('6_0','6_1','6_2','6_3'))
rm(obj)

### Begin WGCNA ####
clusters_list <- c('6_0',
                   '6_1',
                   '6_2',
                   '6_3')

obj_6_only <- SetupForWGCNA(
  obj_6_only,
  gene_select = "fraction",
  fraction = 0.03, 
  wgcna_name = "clust_6" 
)

### Construct metacells
obj_6_only <- MetacellsByGroups(
  seurat_obj = obj_6_only,
  group.by = c("sub.cluster"),
  reduction = 'harmony_wnn.umap', 
  k = 25, 
  max_shared = 10, # again, not sure if these should be altered or not
  ident.group = 'sub.cluster'
)

# normalize metacell expression matrix:
obj_6_only <- NormalizeMetacells(obj_6_only)

### Coexpression Network analysis
obj_6_only <- SetDatExpr(
  obj_6_only,
  group_name = paste0(clusters_list),
  group.by='sub.cluster', 
  assay = 'RNA', 
  layer = 'data' 
)

### Select soft power threshold
obj_6_only <- TestSoftPowers(
  obj_6_only,
  networkType = 'signed' #not sure if this should be altered either
)
# plot the results:
plot_list <- PlotSoftPowers(obj_6_only)
wrap_plots(plot_list, ncol=2)

#guidance is to use greater than or equal to 0.8, which corresponds to
plot_list[[1]][["data"]]
#softpower of 6

obj_6_only <- ConstructNetwork(
  obj_6_only,
  tom_name = 'clust_6', 
  tom_outdir = 'A:/WGCNA_6/',
  overwrite_tom = T
)

PlotDendrogram(obj_6_only, main='clust_6 hdWGCNA Dendrogram')

#model eigengenes
obj_6_only <- ModuleEigengenes(
  obj_6_only,
  group.by.vars="sub.cluster"
)


# harmonized module eigengenes:
hMEs <- GetMEs(obj_6_only)
# module eigengenes:
MEs <- GetMEs(obj_6_only, harmonized=FALSE)

### have to do this because the default path is changed
obj_6_only@misc[["clust_6"]][["wgcna_net"]][["TOMFiles"]] <-'A:/WGCNA_6/clust_6_TOM.rda'
tom = load('A:/WGCNA_6/clust_6_TOM.rda')

# compute eigengene-based connectivity (kME):
obj_6_only <- ModuleConnectivity(
  obj_6_only,
  group.by = 'sub.cluster', group_name =  paste0(clusters_list)
  
)

p <- PlotKMEs(obj_6_only, ncol=5)
p

# get the module assignment table:
modules <- GetModules(obj_6_only) %>% subset(module != 'grey')

#write.csv(modules, 'A:/WGCNA_6/modules_6_subclusters_08_12_2025.csv')

### Get hub genes
hub_df <- GetHubGenes(obj_6_only, n_hubs = 10)
#write.csv(hub_df, 'A:/WGCNA_6/hubs_6_subclusters_08_12_2025.csv')

#saveRDS(obj_6_only,  'A:/WGCNA_6/wgcna_seurat_obj_08_12_2025.rds')
#### Analysis ####
obj_6_only <- readRDS( 'A:/WGCNA_6/wgcna_seurat_obj_08_12_2025.rds')
hub_df = read.csv('A:/WGCNA_6/hubs_6_subclusters_08_12_2025.csv')
modules = read.csv('A:/WGCNA_6/modules_6_subclusters_08_12_2025.csv')
MEs <- GetMEs(obj_6_only, harmonized=FALSE)

### First, what are the modules ###
## lets do GO of the hub genes ##
enrich_df = data.frame()
for(i in unique(hub_df$module)){
  #print(i)
  genes = hub_df$gene_name[hub_df$module==i]
  
  go =clown_go(genes)
  if(length(go$Description)>1){
  message(paste0(i),paste0('_', go$Description)) 
  newd = data.frame(module = i,
                    description = go$Description)
  enrich_df = rbind(newd, enrich_df)
  }
  else{
    message(paste0(i,' no enrichment'))
  }
}

#red is definitely an IEG module its hub genes are fos
#red response to progesterone
#lot of rna pol II stuff not surprising
#>response to steroid hormone
#>response to glucocorticoid

#cyan
#> glial cell-derived neurotrophic factor receptor signaling pathway
#> neuron differentiation
#> neuron dev
#> cell migration
#> positive regulation of neural precursor cell proliferation
#my priors tell me this should be 6_2

#yellow
#> metabolism

#salmon
#> cell signaling type shi

#purple
#> other misc cell stuff

#> greenyellow
#> hormone binidng 
#> semaphorin receptor complex

#blue
#> more cell signaling

# magenta
#> cell migration
#>histone deacetylase binding
#>canonical Wnt signaling pathway
#>chromatin DNA binding

#turquoise
#> development

#tan
#> kinase

#magenta 
#> neuropeptide hormone activity
#> response to testosterone
#>positive regulation of stem cell proliferation


#> So, modules of interest - red, magenta, magenta, cyan
MEs$individual = obj_6_only$individual
MEs$Status = obj_6_only$Status
MEs$sub.cluster = obj_6_only$sub.cluster

plot_module('brown')

brown_test = lmer(brown~Status+(1|individual), data = MEs)
car::Anova(brown_test, type = 'III')

#the one im actually interested in 
plot_module('blue')# seems elevated in D, is it real?
blue_t = lmer(blue~Status+(1|individual), data = MEs)
car::Anova(blue_t, type = 'III')#nope

plot_module('red')
red_t = lmer(red~Status+(1|individual), data = MEs)
car::Anova(red_t, type = 'III')### very!
#woah, what do you make of that 

plot_module('cyan')
cyan_t = lmer(cyan~Status+(1|individual), data = MEs)
car::Anova(cyan_t, type = 'III')### slightly, and this is oly cause of the expandeds which I think is a fluke

#see below
plot_module('turquoise') #meh
turquoise_t = lmer(turquoise~Status+(1|individual), data = MEs)
car::Anova(turquoise_t, type = 'III')
#yeah that seems right

### I am now fully interested in red, what clusters express it the most ####
plot_module_cluster('red')
plot_module_cluster('blue')
plot_module_cluster('turquoise') # this is why it is interesting because females are enriched for this cluster

#### Which modules are enriched for DEGs ####
deg_data = read_csv("DEG Outputs/FINAL_degs_classified_08_11_2025.csv")
deg_data_6 = subset(deg_data, cluster == 6)

#mark the degs
modules$isDEG = ifelse(modules$gene_name %in% deg_data_6$gene, TRUE, FALSE)

#summarize
module_deg_table = table(modules$module, modules$isDEG)%>%as.data.frame.matrix()
module_deg_table$diff = module_deg_table$`FALSE`-  module_deg_table$`TRUE`
module_deg_table$module = rownames(module_deg_table)

#perform fishers exact test
fish_data =data.frame()
for(module in module_deg_table$module){
  fish_matrix = matrix(NA, 2, 2)
  fish_matrix[1,1] = sum(module_deg_table$`TRUE`[module_deg_table$module==module])
  fish_matrix[1,2] = sum(module_deg_table$`FALSE`[module_deg_table$module==module])
  
  fish_matrix[2,1] = sum(module_deg_table$`TRUE`[module_deg_table$module!=module])
  fish_matrix[2,2] = sum(module_deg_table$`FALSE`[module_deg_table$module!=module])
  
  test = fisher.test(fish_matrix)
  newd = data.frame(module =module,
                    p = test$p.value)
  fish_data = rbind(fish_data, newd)
}
fish_data$isSignif = ifelse(fish_data$p<0.05, '*', NA)
## again, blue is the important one..., also turquoise?

#here can I go from here, blue is enriched for DEGs but itself is not differentially
#expressed, what does GO of the whole module have to say? Does it matter?

clown_go(modules$gene_name[modules$module=='blue'])%>%dotplot()
clown_go(hub_df$gene_name[hub_df$module=='blue'])%>%dotplot()

clown_go(modules$gene_name[modules$module=='turquoise'])%>%dotplot()
clown_go(hub_df$gene_name[hub_df$module=='turquoise'])%>%dotplot()

### is there a bimodal distribution of cytotrace scores in 6_3 #####
library(CytoTRACE)
cyto= CytoTRACE(as.matrix(obj_6_only@assays$RNA$data))
obj_6_only$cyto = cyto$CytoTRACE
obj_6_only@meta.data$Status = factor(obj_6_only@meta.data$Status, levels = c('NRM',"M",'D','E','NF',"F"))

ggplot(obj_6_only@meta.data, aes(x = sub.cluster, y = cyto, color = Status))+
  geom_violin()

ggplot(obj_6_only@meta.data[obj_6_only@meta.data$Status %in% c('M','D','F'),], aes(x = sub.cluster, y = cyto))+
  geom_violin()+
  stat_summary(geom = 'crossbar', fun = 'mean', aes(color = Status))

ggplot(obj_6_only@meta.data[obj_6_only@meta.data$sub.cluster=='6_3',], aes(x = Status, y = cyto, color = Status))+
  geom_violin()+
  geom_point()



