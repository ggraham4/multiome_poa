#> Subclustering oligos
#> 
#> GJG 09/09/2025
#> 
#libraries

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

  `%notin%` <- Negate(`%in%`)
  
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

prop_cluster_plot=
function(object, gene, cluster, clustering = 'final_clusters'){
    library(stringr)
    library(forcats)
      options(dplyr.summarise.inform = FALSE)

    counts <- t(as.matrix(object@assays$RNA$data[,object@meta.data[[clustering]] == cluster]))
  Counts_of_interest <- as.data.frame(counts[,gene]>0) #binarize
    Counts_of_interest$expression <- Counts_of_interest[,1]
  Counts_of_interest$individual <- object@meta.data$individual[object@meta.data[[clustering]] == cluster]
    Counts_of_interest$Status <- object@meta.data$individual[object@meta.data[[clustering]] == cluster]

  results <-Counts_of_interest%>%
    group_by(individual, Status)%>%
    summarize(mean = mean(expression),
              se = sd(expression)/sqrt(n()))
  results$Sex <- results$Status
  
  
  results$Sex <- str_sub(results$individual, -1)
  results$Sex[results$individual == 'T17D'] = 'NF'
  results$Sex[results$individual == 'A12D'] = 'E'
  results$Sex[results$individual == 'T11D'] = 'E'
  results$Sex[results$individual == 'GH'] = 'NRM'
  
  results$factor <- ifelse(results$Sex == "NRM", 0, NA)
  results$factor <- ifelse(results$Sex == "M", 1, results$factor)
  results$factor <- ifelse(results$Sex == "D", 2, results$factor)
  results$factor <- ifelse(results$Sex == "E", 3, results$factor)
  results$factor <- ifelse(results$Sex == "NF", 4, results$factor)
  results$factor <- ifelse(results$Sex == "F", 5, results$factor)
  results$individual <- fct_reorder(results$individual, results$factor)
  
  results$Sex <- factor(results$Sex, levels = c('NRM','M','D','E','NF','F'))
  plot <- ggplot(results, aes(x = individual, y = mean, color = Sex))+
    geom_boxplot(aes(group = Sex, fill = Sex), alpha = 0.25, outlier.shape = NA)+
    geom_point()+
    geom_pointrange(aes(x = individual, y = mean, ymin = mean - se, ymax = mean+se))+
    theme_classic()+
    labs(x  ='FishID', y = '% Expressing', title = paste0(gene,'_cluster_',cluster))+
    theme(axis.text.x = element_text(angle = -45, vjust = 1, hjust=0))
  return(plot)
}


}
#read in data
obj = readRDS("~/Desktop/optimal_clustering_rna_only.rds")

### first of all, are DEGs enriched for anything ----
deg_data = read_csv("DEG Outputs/FINAL_degs_classified_08_11_2025.csv")
deg_data_2 = subset(deg_data, cluster == 2)
clown_go(deg_data_2$gene[deg_data_2$second_word =='Downregulated'])%>%
  dotplot() 

clown_go(deg_data_2$gene[deg_data_2$second_word =='Upregulated'])%>%
  dotplot()

clown_go(deg_data_2$gene[deg_data_2$first_word =='Early'])%>%
  dotplot()

clown_go(deg_data_2$gene[deg_data_2$first_word =='Late'])%>%
  dotplot()

clown_go(deg_data_2$gene[deg_data_2$first_word =='Transiently'])%>%
  dotplot()

#clown_go(deg_data_2$gene[deg_data_2$first_word =='Progressively'])%>%
#  dotplot()

####  subcluster and subset ----
Idents(obj) = 'final_clusters'
olig <- FindSubCluster(obj, 2, 'harmony.wsnn')
olig = subset(olig, final_clusters == 2)

Idents(olig) <- 'sub.cluster'

DimPlot(olig, label = T, reduction = 'harmony_wnn.umap')
#2 subclusters

olig$Status = factor(olig$Status, levels = c('NRM', 'M','D','E','NF','F'))
### My first question, do they differ in proportion ----
cells_ind = olig@meta.data%>%
  group_by(individual)%>%
  summarize(n_cells = n())

cells_sub_ind = olig@meta.data%>%
  group_by(individual, Status, sub.cluster)%>%
  summarize(n_cells_in = n())

cells_total = cells_ind%>%
  right_join(cells_sub_ind, by = 'individual')
cells_total$prop = cells_total$n_cells_in/cells_total$n_cells

ggplot(subset(cells_total, Status %in% c('M','D',"F")), aes(x = sub.cluster, y = prop, color = Status))+
  geom_boxplot(alpha = 0)+
  geom_point(position = position_jitterdodge(1))
#WOAH that is a clear signal

## fewer dominant cells in 1_3
cells_2_0 = subset(cells_total, sub.cluster == '2_0' & Status %in% c('M','D','F'))
matrix_2_0 = cbind(cells_2_0$n_cells_in,cells_2_0$n_cells-cells_2_0$n_cells_in)

model_2_0 = glm(matrix_2_0~Status, family = 'binomial', data = cells_2_0)
anova(model_2_0, test ='Chisq')

pairs(emmeans(model_2_0, 'Status'), adjust ='none')


### What markers differentiate the subclusters ----
olig_markers = FindAllMarkers(olig, only.pos = T)
olig_markers_good = subset(olig_markers,  p_val_adj < 0.05)

clown_go(olig_markers_good$gene[olig_markers_good$cluster=='2_0'])%>%dotplot()
clown_go(olig_markers_good$gene[olig_markers_good$cluster=='2_1'])%>%dotplot()


#### CytoTRACE, I bet 4 is the highest ----
cyto = CytoTRACE(olig@assays$RNA$data%>%as.matrix())
olig$cyto = cyto$CytoTRACE
FeaturePlot(olig, 'cyto')
DimPlot(olig)

cyto_plot = olig@meta.data%>%
  group_by(individual, Status, sub.cluster)%>%
  summarize(mean_cyto = mean(cyto),
            se_cyto = sd(cyto)/sqrt(n()))

ggplot(cyto_plot, aes(x = sub.cluster,y = mean_cyto, color = Status))+
  geom_boxplot()

### WGCNA ----
clusters_list <- c('2_0',
                   '2_1')

olig <- SetupForWGCNA(
  olig,
  gene_select = "fraction",
  fraction = 0.05, 
  wgcna_name = "olig" 
)

### Construct metacells
olig <- MetacellsByGroups(
  seurat_obj = olig,
  group.by = c("sub.cluster"),
  reduction = 'harmony_wnn.umap', 
  k = 25, 
  max_shared = 10,
  ident.group = 'sub.cluster',
  min_cells = 90 
)

# normalize metacell expression matrix:
olig <- NormalizeMetacells(olig)

### Coexpression Network analysis
olig <- SetDatExpr(
  olig,
  group_name = clusters_list,
  group.by='sub.cluster', 
  assay = 'RNA', 
  layer = 'data' 
)

### Select soft power threshold
olig <- TestSoftPowers(
  olig,
  networkType = 'signed' #not sure if this should be altered either
)

# plot the results:
plot_list <- PlotSoftPowers(olig)
wrap_plots(plot_list, ncol=2)

olig <- ConstructNetwork(
  olig,
  tom_name = 'olig', 
  tom_outdir = '/Users/ggraham/WGCNA_olig/',
  overwrite_tom = T
)

PlotDendrogram(olig, main='clust_2 hdWGCNA Dendrogram')

#model eigengenes
olig <- ModuleEigengenes(
  olig,
  group.by.vars="sub.cluster"
)

# harmonized module eigengenes:
hMEs <- GetMEs(olig)
# module eigengenes:
MEs <- GetMEs(olig, harmonized=FALSE)

olig@misc[["olig"]][["wgcna_net"]][["TOMFiles"]] <-'/Users/ggraham/WGCNA_olig/olig_TOM.rda'

# compute eigengene-based connectivity (kME):
olig <- ModuleConnectivity(
  olig,
  group.by = 'sub.cluster', group_name =  clusters_list
  
)

p <- PlotKMEs(olig, ncol=5)
p

# get the module assignment table:
modules <- GetModules(olig) %>% subset(module != 'grey')

### Get hub genes
hub_df <- GetHubGenes(olig, n_hubs = 10)

#enrich hub genes
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

ModuleRadarPlot(
  olig,
  group.by = 'sub.cluster',
  axis.label.size=4,
  grid.label.size=4
)

ModuleRadarPlot(
  olig,
  group.by = 'Status',
  axis.label.size=4,
  grid.label.size=4
)

###DEGs ----
fake_degs = FindMarkers(subset(olig, sub.cluster=='2_0'), ident.1 = 'F', group.by = 'Status')
fake_degs2 = FindMarkers(subset(olig, sub.cluster=='2_0'), ident.1 = 'D', group.by = 'Status')

#######I'm not even convinced this cluster should be subclustered looking at the marker genes====
#> for the whole cluster, do the number of cells differ by sex
 cells_ind_overall = obj@meta.data%>%
  group_by(individual)%>%
  summarize(n_cells = n())

cells_sub_ind_overall = olig@meta.data%>%
  group_by(individual, Status, final_clusters)%>%
  summarize(n_cells_in = n())

cells_total_overall = cells_ind_overall%>%
  right_join(cells_sub_ind_overall, by = 'individual')

cells_total_overall$prop = cells_total_overall$n_cells_in/cells_total_overall$n_cells

ggplot(subset(cells_total_overall, Status %in% c('M','D',"F")&final_clusters == 2), aes(x = final_clusters, y = prop, color = Status))+
  geom_boxplot(alpha = 0)+
  geom_point(position = position_jitterdodge(.2))

cells_2 = subset(cells_total_overall, final_clusters == 2 & Status %in% c('M','D','F'))
matrix_2 = cbind(cells_2$n_cells_in,cells_2$n_cells-cells_2$n_cells_in)

model_2 = glm(matrix_2~Status, family = 'binomial', data = cells_2)
anova(model_2, test ='Chisq')

hist(model_2$residuals)

pairs(emmeans(model_2, 'Status'), adjust ='none')

hist(matrix_2[,1]/matrix_2[,2])

library(glmmTMB)
beta_model = glmmTMB(matrix_2~Status, data = cells_2, family = betabinomial())
car::Anova(beta_model, 3)

