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
obj  = readRDS("~/Desktop/optimal_clustering_rna_only.rds")

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
  tom_outdir = '/Users/ggraham/Desktop/Misc/',
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
obj_6_only@misc[["clust_6"]][["wgcna_net"]][["TOMFiles"]] <-'/Users/ggraham/Desktop/Misc/clust_6_TOM.rda'
tom = load('/Users/ggraham/Desktop/Misc/clust_6_TOM.rda')

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

MEs$individual = obj_6_only$individual
MEs$Status = obj_6_only$Status
MEs$sub.cluster = obj_6_only$sub.cluster

for(i in unique(colnames(MEs))){
 print( plot_module(i))
}

# are any modules correlated with cyto #####
library(CytoTRACE)
cyto= CytoTRACE(as.matrix(obj_6_only@assays$RNA$data))
obj_6_only$cyto = cyto$CytoTRACE
obj_6_only@meta.data$Status = factor(obj_6_only@meta.data$Status, levels = c('NRM',"M",'D','E','NF',"F"))

# combine CytoTRACE + MEs into one dataframe
df <- cbind(MEs, cyto = obj_6_only$cyto)

# pivot longer (modules into one column)
df_long <- df %>%
  pivot_longer(
    cols = colnames(MEs)[1:13],
    names_to = "module",
    values_to = "ME_value"
  )

# plot
ggplot(df_long, aes(x = ME_value, y = cyto, color = module)) +
  theme_classic()+
  geom_smooth()

ggplot(MEs, aes(x = purple, y = obj_6_only$cyto)) +
  geom_point(alpha = 0.6) +
  theme_classic()+
  geom_smooth()

# yellow, salmon, pink, blue, greenyellow, tan, purple

enrich_df$description[enrich_df$module=='yellow']
# some kind of metabolism

enrich_df$description[enrich_df$module=='salmon']
# unclear

enrich_df$description[enrich_df$module=='pink'] #*
# pretty obv

enrich_df$description[enrich_df$module=='blue']
# signaling

enrich_df$description[enrich_df$module=='greenyellow']
# this could be something

enrich_df$description[enrich_df$module=='tan']
# methylation / epigenome

enrich_df$description[enrich_df$module=='purple']
# maybe ecm


plot_module('pink')
pink_t = lmer(pink~Status+(1|individual), data = MEs)
car::Anova(pink_t, type = 'III') # enriched for nervous system development


plot_module('blue')
blue_t = lmer(blue~Status+(1|individual), data = MEs)
car::Anova(blue_t, type = 'III')




enrich_df2 = data.frame()
for(i in unique(hub_df$module)){
  #print(i)
  genes = modules$gene_name[modules$module==i]
  
  go =clown_go(genes)
  if(length(go$Description)>1){
  message(paste0(i),paste0('_', go$Description)) 
  newd = data.frame(module = i,
                    description = go$Description)
  enrich_df2 = rbind(newd, enrich_df2)
  }
  else{
    message(paste0(i,' no enrichment'))
  }
}

ggplot(MEs, aes(x = pink, y = obj_6_only$cyto)) +
  geom_point(alpha = 0.6) +
  theme_classic()+
  geom_smooth()

ggplot(MEs, aes(x = purple, y = obj_6_only$cyto)) +
  geom_point(alpha = 0.6) +
  theme_classic()+
  geom_smooth() # ok first of all woah, super interesting and worth a look later


ggplot(MEs, aes(x = grey, y = obj_6_only$cyto)) +
  geom_point(alpha = 0.6) +
  theme_classic()+
  geom_smooth() # ok I think there are some detection effects going on here

ggplot(MEs, aes(x = grey, y = obj_6_only$nCount_RNA)) +
  geom_point(alpha = 0.6) +
  theme_classic()+
  scale_y_log10()+
  geom_smooth()  # yeah there is 



# how to avoid a "linear regression is my passion"

passion_pink = lmer(pink~obj_6_only$cyto+obj_6_only$nCount_RNA+(1|individual),
                    data = MEs,
                    control = lmerControl(autoscale=T))
car::Anova(passion_pink) # it says yes here


# compare to a module I dont think is


passion_grey =lmer(grey~obj_6_only$cyto*obj_6_only$nCount_RNA+(1|individual),
                    data = MEs,
                    control = lmerControl(autoscale=T))
car::Anova(passion_grey, 3) # it says yes here






ggplot(MEs, aes(x = blue, y = obj_6_only$cyto)) +
  geom_point(alpha = 0.6) +
  theme_classic()+
  geom_smooth()


###  I think there is something with pink ###
# synapse assembly and brain development, specific genes?
# hmx tfs part of this module
#foxg1a
#lhx5


# blue seems like a more mature neuron, should do fishers exact test to see if any modules enriched for degs


sofiya_data = data.frame(cyto = obj_6_only$cyto, 
                    nCount_RNA = obj_6_only$nCount_RNA, 
                    pink = MEs$pink,
                    individual = obj_6_only$individual)

write.csv(sofiya_data, 'sofiya_data.csv')
