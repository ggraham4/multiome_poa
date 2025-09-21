#> Subclustering radial Glia
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
deg_data_1 = subset(deg_data, cluster == 1)
clown_go(deg_data_1$gene[deg_data_1$second_word =='Downregulated'])%>%
  dotplot()
clown_go(deg_data_1$gene[deg_data_1$second_word =='Upregulated'])%>%
  dotplot()
clown_go(deg_data_1$gene[deg_data_1$first_word =='Early'])%>%
  dotplot()
clown_go(deg_data_1$gene[deg_data_1$first_word =='Late'])%>%
  dotplot()
clown_go(deg_data_1$gene[deg_data_1$first_word =='Transiently'])%>%
  dotplot()
clown_go(deg_data_1$gene[deg_data_1$first_word =='Progressively'])%>%
  dotplot()

####  subcluster and subset ----
Idents(obj) = 'final_clusters'
rgc <- FindSubCluster(obj, 1, 'harmony.wsnn')
rgc = subset(rgc, final_clusters == 1)

Idents(rgc) <- 'sub.cluster'

DimPlot(rgc, label = T, reduction = 'harmony_wnn.umap')

rgc$Status = factor(rgc$Status, levels = c('NRM', 'M','D','E','NF','F'))
### My first question, do they differ in proportion ----
cells_ind = rgc@meta.data%>%
  group_by(individual)%>%
  summarize(n_cells = n())

cells_sub_ind = rgc@meta.data%>%
  group_by(individual, Status, sub.cluster)%>%
  summarize(n_cells_in = n())

cells_total = cells_ind%>%
  right_join(cells_sub_ind, by = 'individual')
cells_total$prop = cells_total$n_cells_in/cells_total$n_cells

ggplot(subset(cells_total, Status %in% c('M','D',"F")), aes(x = sub.cluster, y = prop, color = Status))+
  geom_boxplot(alpha = 0)+
  geom_point(position = position_jitterdodge(1))

ggplot(subset(cells_total, sub.cluster =='1_0'), aes(x = sub.cluster, y = prop, color = Status))+
  geom_boxplot(alpha = 0)+
  geom_point(position = position_jitterdodge(1))
## i would bet there is a difference here between m and f

ggplot(subset(cells_total, sub.cluster =='1_3'), aes(x = sub.cluster, y = prop, color = Status))+
  geom_boxplot(alpha = 0)+
  geom_point(position = position_jitterdodge(1))
#def something with dominants here

### my only two hypotheses for these are there are fewer female cells in 1_0, and
## fewer dominant cells in 1_3
cells_1_0 = subset(cells_total, sub.cluster == '1_0' & Status %in% c('M','D','F'))
matrix_1_0 = cbind(cells_1_0$n_cells_in,cells_1_0$n_cells-cells_1_0$n_cells_in)

model_1_0 = glmer(matrix_1_0~Status +(1|individual), family = 'binomial', data = cells_1_0)
car::Anova(model_1_0, 3)
#nope

cells_1_3 = subset(cells_total, sub.cluster == '1_3' & Status %in% c('M','D','F'))
matrix_1_3 = cbind(cells_1_3$n_cells_in,cells_1_3$n_cells-cells_1_3$n_cells_in)

model_1_3 = glmer(matrix_1_3~Status +(1|individual), family = 'binomial', data = cells_1_3)
car::Anova(model_1_3, 3)
#YES
pairs(emmeans(model_1_3, 'Status'), adjust ='none') # dominants both


### What markers differentiate the subclusters ----
rgc_markers = FindAllMarkers(rgc)
rgc_markers_good = subset(rgc_markers, pct.1>pct.2 & p_val_adj < 0.05)

clown_go(rgc_markers_good$gene[rgc_markers_good$cluster=='1_0'])%>%dotplot()
clown_go(rgc_markers_good$gene[rgc_markers_good$cluster=='1_1'])%>%dotplot()
clown_go(rgc_markers_good$gene[rgc_markers_good$cluster=='1_2'])%>%dotplot()
clown_go(rgc_markers_good$gene[rgc_markers_good$cluster=='1_3'])%>%dotplot()
clown_go(rgc_markers_good$gene[rgc_markers_good$cluster=='1_4'])%>%dotplot() ### likely proliferating
clown_go(rgc_markers_good$gene[rgc_markers_good$cluster=='1_5'])%>%dotplot()

DotPlot(rgc, 'LOC111577263')+coord_flip()#### strongly in 1_3, 1_2, and 1_0

  clown_go = readRDS("Functions/clown_go2")  

go_1_0 = FindMarkers(rgc, '1_0')%>%
  subset(p_val_adj<0.05)%>%
  rownames()%>%
  clown_go%>%
  dotplot()
go_1_0
go_1_1 = FindMarkers(rgc, '1_1')%>%
  subset(p_val_adj<0.05)%>%
  rownames()%>%
  clown_go%>%
  dotplot()
go_1_1
go_1_2 = FindMarkers(rgc, '1_2')%>%
  subset(p_val_adj<0.05)%>%
  rownames()%>%
  clown_go%>%
  dotplot()
go_1_2

go_1_3 = FindMarkers(rgc, '1_3')%>%
  subset(p_val_adj<0.05)%>%
  rownames()%>%
  clown_go%>%
  dotplot()
go_1_3

go_1_4 = FindMarkers(rgc, '1_4')%>%
  subset(p_val_adj<0.05)%>%
  rownames()%>%
  clown_go%>%
  dotplot()
go_1_4

go_1_5 = FindMarkers(rgc, '1_5')%>%
  subset(p_val_adj<0.05)%>%
  rownames()%>%
  clown_go%>%
  dotplot()
go_1_5
### difference between 3 and 1 ----
  clown_go = readRDS("Functions/clown_go2")  

markers_1_3 <- FindMarkers(rgc, '1_3', '1_1')
clown_go(rownames(markers_1_3[markers_1_3$p_val_adj<0.001 & markers_1_3$pct.1>markers_1_3$pct.2 ,]))%>%
  dotplot()

clown_go(rownames(markers_1_3[markers_1_3$p_val_adj<0.001 & markers_1_3$pct.1<markers_1_3$pct.2 ,]))%>%
  dotplot()


## can I ID the ependymal lineage
ep_markers <- FindMarkers(obj, '15', '1')

#### CytoTRACE, I bet 4 is the highest ----
cyto = CytoTRACE(rgc@assays$RNA$data%>%as.matrix())
rgc$cyto = cyto$CytoTRACE
FeaturePlot(rgc, 'cyto')
DimPlot(rgc) #woah 1 and 3

cyto_plot = rgc@meta.data%>%
  group_by(individual, Status, sub.cluster)%>%
  summarize(mean_cyto = mean(cyto),
            se_cyto = sd(cyto)/sqrt(n()))

ggplot(cyto_plot, aes(x = sub.cluster,y = mean_cyto, color = Status))+
  geom_boxplot()
# 3 seems more immature in doms, 5 may be more mature

cyto_1_3 = lmer(cyto~Status+(1|individual), data = subset(rgc@meta.data, sub.cluster =='1_3' ))
car::Anova(cyto_1_3, 3)
#nope

cyto_1_5 = lmer(cyto~Status+(1|individual), data = subset(rgc@meta.data, sub.cluster =='1_5'))
car::Anova(cyto_1_5, 3)

ggplot(cyto_plot, aes(x = Status,y = mean_cyto, color = Status))+
  geom_boxplot()
### WGCNA ----
clusters_list <- c('1_0',
                   '1_1',
                   '1_2',
                   '1_3',
                   '1_4',
                   '1_5')

rgc <- SetupForWGCNA(
  rgc,
  gene_select = "fraction",
  fraction = 0.05, 
  wgcna_name = "rgc" 
)

### Construct metacells
rgc <- MetacellsByGroups(
  seurat_obj = rgc,
  group.by = c("sub.cluster"),
  reduction = 'harmony_wnn.umap', 
  k = 25, 
  max_shared = 10,
  ident.group = 'sub.cluster',
  min_cells = 90 
)

# normalize metacell expression matrix:
rgc <- NormalizeMetacells(rgc)

### Coexpression Network analysis
rgc <- SetDatExpr(
  rgc,
  group_name = clusters_list,
  group.by='sub.cluster', 
  assay = 'RNA', 
  layer = 'data' 
)

### Select soft power threshold
rgc <- TestSoftPowers(
  rgc,
  networkType = 'signed' #not sure if this should be altered either
)

# plot the results:
plot_list <- PlotSoftPowers(rgc)
wrap_plots(plot_list, ncol=2)

rgc <- ConstructNetwork(
  rgc,
  tom_name = 'rgc', 
  tom_outdir = '/Users/ggraham/WGCNA_rgc/',
  overwrite_tom = T
)

PlotDendrogram(rgc, main='clust_1 hdWGCNA Dendrogram')

#model eigengenes
rgc <- ModuleEigengenes(
  rgc,
  group.by.vars="sub.cluster"
)

# harmonized module eigengenes:
hMEs <- GetMEs(rgc)
# module eigengenes:
MEs <- GetMEs(rgc, harmonized=FALSE)

rgc@misc[["rgc"]][["wgcna_net"]][["TOMFiles"]] <-'/Users/ggraham/WGCNA_rgc/rgc_TOM.rda'

# compute eigengene-based connectivity (kME):
rgc <- ModuleConnectivity(
  rgc,
  group.by = 'sub.cluster', group_name =  clusters_list
  
)

p <- PlotKMEs(rgc, ncol=5)
p

# get the module assignment table:
modules <- GetModules(rgc) %>% subset(module != 'grey')

### Get hub genes
hub_df <- GetHubGenes(rgc, n_hubs = 10)

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
  rgc,
  group.by = 'sub.cluster',
  axis.label.size=4,
  grid.label.size=4
)

ModuleRadarPlot(
  rgc,
  group.by = 'Status',
  axis.label.size=4,
  grid.label.size=4
)

# 1_1 and 1_3 are both enriched for the turquose, 1_3 is also enriched for yellow

clown_go(modules$gene_name[modules$color=='turquoise'])%>%dotplot()
# i swear to god if I get the dendrite GO term one more time

clown_go(modules$gene_name[modules$color=='red'])%>%dotplot()
clown_go(hub_df$gene_name[hub_df$module=='red'])%>%dotplot() # yeah I mean cmom


clown_go(modules$gene_name[modules$color=='yellow'])%>%dotplot()
# seems like yellow is more proliferating cells?

clown_go(modules$gene_name[modules$color=='brown'])%>%dotplot()

clown_go(modules$gene_name[modules$color=='green'])%>%dotplot()

# 3 expresses everything so its enriched in all modules basically with strongest
# expression of turquoise along with 1
MEs$Status = rgc$Status
MEs$individual = rgc$individual
MEs$sub.cluster = rgc$sub.cluster

plot_module('turquoise')
plot_module('yellow')

plot_module_cluster('turquoise')
### def immature 

plot_module_cluster('red') # sure enough this looks like my neuronal? lineage
plot_module('red') #something there perhaps, it shows the opposite pattern of arom
#a hub df of red


plot_module('blue')


dmes = data.frame()
for(i in unique(modules$color)){
  formula_string <- paste0(i, " ~ Status +(1|individual)")%>%as.formula()
  model = lmer(formula_string, data = subset(MEs, Status %in% c('M','D','F')))
  p =car::Anova(model, 3)
  newd = data.frame(module = i,
                    pval = p$`Pr(>Chisq)`[2])
  dmes = rbind(dmes, newd)
  
}
#nope

#arom is part of the red module

### 3 and 2 (and maybe 0) seem to be related in the same lineage, can I prove? ----
ModuleRadarPlot(
  rgc,
  group.by = 'sub.cluster',
  axis.label.size=4,
  grid.label.size=4
)

clown_go(modules$gene_name[modules$color=='red'])%>%
  dotplot() # neurogenesis
clown_go(modules$gene_name[modules$color=='yellow'])%>%
  dotplot() # cell diff
clown_go(modules$gene_name[modules$color=='blue'])%>%
  dotplot() # transcription

#1_3 and 1_1
clown_go(modules$gene_name[modules$color=='turquoise'])%>%
  dotplot()
clown_go(modules$gene_name[modules$color=='brown'])%>%
  dotplot()

clown_go(modules$gene_name[modules$color=='green'])%>%
  dotplot()

### do the sexes differ in red or yellow 

#### I want to make a radial glial index ---

""
'
https://www.cell.com/fulltext/S0092-8674(15)01124-1

radial glia
sox2,
pax6,
gli3,
slc1a3,
PDGFD,
vim,
hes1

intermediate progenitors
eomes,
elavl4,
neurog1,
neurod1,
neurod4, 
PPP1R17,
penk

'
""


#geneList = c( #just using the ieg list as a test
#  'LOC111583367',
#  'egr1',
#  'npas4a'
#) #goddamn this really shows that zacks method works for IEGs becuase
# whenever I try for other sets of genes it does not give me shit

calculateCoexpressedGenes = function(seurat_obj2, geneList, clustering = 'sub.cluster'){
  ### the goal is to split cells in a cluster into + and negative for the gene list,
  ### then find other genes coexpressed with those
  
  `%notin%` = Negate(`%in%`)
  seurat_obj <- seurat_obj2
  seurat_obj$clustering = seurat_obj[[clustering]]
  
  seurat_obj$positive = ifelse(seurat_obj@assays$RNA$data[geneList[1],]>0,T, NA)
  for(i in 2:length(geneList)){
    if(geneList[i] %in% rownames(seurat_obj@assays$RNA$data)){
      seurat_obj$positive = ifelse(seurat_obj@assays$RNA$data[geneList[i],]>0,T, seurat_obj$positive)
    }}
    seurat_obj$positive = ifelse(is.na(seurat_obj$positive), F, seurat_obj$positive)

  
  marker_all<- data.frame()
  for(cluster in unique(seurat_obj$clustering)){
    temp_obj = subset(seurat_obj, clustering == cluster)
    Idents(temp_obj) = 'positive'
    
    if(sum(temp_obj$positive)>3){
    
    markers <- FindAllMarkers(temp_obj,
                              assay = "RNA",
                              group.by = 'positive',
                              logfc.threshold = 0,
                              min.pct = 1/nrow(temp_obj@meta.data))
    markers$sub.cluster= cluster
    marker_all = rbind(markers, marker_all)
    }
  }
  #extract significant   genes
  marker_all_signif <- marker_all[marker_all$p_val_adj<0.05 &marker_all$cluster==T ,]
  
  #genes in over half of clusters
  marker_genes_half = marker_all_signif%>%
    group_by(gene)%>%
    summarize(n_clusters = n())%>%
  subset(n_clusters >= length(unique(seurat_obj$clustering))/2
  ) 
  
  message('Found ', nrow(marker_genes_half) - length(geneList), ' new markers')
  print(marker_genes_half$gene[marker_genes_half$gene %notin% geneList])
  
  all_markers = marker_genes_half$gene
  #module_score time
  seurat_obj = AddModuleScore(seurat_obj, 
                              features = list(all_markers),
                              name = 'coexModule')
  
  return(seurat_obj)
}

ieg = calculateCoexpressedGenes(rgc,c( #just using the ieg list as a test
  'LOC111583367',
  'egr1',
  'npas4a'
))

DotPlot(ieg, 'coexModule1')+
  coord_flip()
#why tf does 1_3 of course have the highest coexpression of IEGs


progenitor = calculateCoexpressedGenes(rgc,c( 
  'neurod4',
  'LOC111575776', #penk
  'eomesb',
  'neurod1'
))

DotPlot(progenitor, 'coexModule1')+
  coord_flip()

radial = calculateCoexpressedGenes(rgc,c( 
  'sox2',
  'pax6b', 
  'slc1a3b',
  'slc1a3a',
  'gli3',
  'pdgfd',
  'gfap',
  'LOC111579533' #vim
  
))

DotPlot(radial, 'coexModule1')+
  coord_flip()

# no idea 


### pseudotime ----
  library(monocle3)
gene_meta_data <- data.frame(row.names = rownames(rgc@assays$RNA$data),
                                                  gene_short_name=
                                                    rownames(rgc@assays$RNA$data)
)
                                                  
cds <- new_cell_data_set(rgc@assays$RNA$data,
                         cell_metadata = rgc@meta.data,
                         gene_metadata =gene_meta_data)
## Step 1: Normalize and pre-process the data
cds <- preprocess_cds(cds, num_dim = 60)

## Step 3: Reduce the dimensions using UMAP
cds <- reduce_dimension(cds)

## Step 4: Cluster the cells
cds <- cluster_cells(cds)
                    
## Step 5: Learn a graph
cds <- learn_graph(cds)

plot_cells(cds,
           color_cells_by='sub.cluster',
              group_label_size = 4)+
  theme(legend.position = 'top')

cds <- order_cells(cds) # choose the bottom one of the big island that is 1_1

plot_cells(cds,
           #color_cells_by='pseudotime',
           color_cells_by  = 'sub.cluster',
             group_label_size = 4
)+
  theme(legend.position = 'top')
#suggests 1-1 is basal to 1_3 which is basal to 1_2 and 1_0, though I would like to know best
# practices for monocle

rgc$pseudotime = pseudotime(cds)
rgc$pseudotime  = ifelse(rgc$pseudotime == Inf, max(rgc$pseudotime[is.finite(rgc$pseudotime)]), rgc$pseudotime )

FeaturePlot(rgc, 'pseudotime', label =T)
DotPlot(rgc,'pseudotime')+
  coord_flip()

ggplot(rgc@meta.data, aes(x = sub.cluster, y = pseudotime, color = Status) )+
  geom_boxplot()

ggplot(rgc@meta.data, aes(x = Status, y = pseudotime, color = Status) )+
  geom_violin()

#IDK pseudotime is sketchy tbh I'm not gonna read too much into this

# does pseudotime correlate with cytotrace
ggplot(subset(rgc@meta.data, pseudotime != max(rgc$pseudotime)), aes(x = cyto, y = pseudotime))+
    geom_point(aes(color = sub.cluster))+
  geom_smooth(method = 'lm')

coef = lmer(pseudotime~cyto+(1|individual), data =subset(rgc@meta.data, pseudotime != max(rgc$pseudotime)))
summary(coef)
car::Anova(coef)




####proportion of cells expressing aromatase in a cluster----

arom_data = ifelse(rgc@assays$RNA$data['LOC111577263',]>0, T, F)
rgc$arom_ex = arom_data

prop_arompos = rgc@meta.data%>%
  group_by(individual, Status, sub.cluster)%>%
  summarize(n_pos = sum(arom_ex ==T),
            n_cells = n())%>%
  mutate(diff = n_cells - n_pos)

ggplot(prop_arompos, aes(x = Status, y = (n_pos/n_cells), color = Status))+
  geom_point()+
  geom_boxplot(alpha = 0)+
  facet_wrap(~sub.cluster)

arom_prop_data = data.frame()
for(i in 0:5){
  cluster = paste0('1_',i)
  temp_data = subset(prop_arompos, sub.cluster==cluster)
  matrix = cbind(temp_data$n_pos, temp_data$diff)
  model = glmer(matrix~Status+(1|individual), data = temp_data, family = 'binomial') # 4 throws an error
  av_ = car::Anova(model, 3)
  av_p = av_$`Pr(>Chisq)`[2]
  newd = data.frame(cluster =cluster, 
                    p = av_p,
                    singular = isSingular(model))
  arom_prop_data = rbind(newd, arom_prop_data)
  print(cluster)
  print(pairs(emmeans(model, 'Status'), adjust = 'none'))
}

#1_3 is stem cell
#does aromatase expression associate with IEG expression?

radial$ieg = ieg$coexModule1
ggplot(radial@meta.data,aes(x = ieg, y = radial@assays$RNA$data['LOC111577263',], color =sub.cluster))+
  geom_point()
#i think nothing

dat = data.frame(x = radial$ieg, y = radial@assays$RNA$data['LOC111577263',])%>%
  subset(y>0)
ggplot(dat,aes(x = x, y = y))+
  geom_point()+
  coord_flip()

### what do aromatase+ cells also express ----
radial$arom = ifelse(radial@assays$RNA$data['LOC111577263',]>0, T, F)
arom_markers = FindAllMarkers(radial, group.by = 'arom')
#notable that serotoninr receptor 2a is a marker
  clown_go = readRDS("Functions/clown_go")  

clown_go(arom_markers$gene[arom_markers$cluster==TRUE &
                             arom_markers$p_val_adj<0.05&
                            arom_markers$pct.1>arom_markers$pct.2])%>%dotplot()

clown_go(arom_markers$gene[arom_markers$cluster==F &
                             arom_markers$p_val_adj<0.05&
                            arom_markers$pct.1>arom_markers$pct.2])%>%dotplot()

FeaturePlot(radial, 'cyto')+
  facet_wrap(~radial@meta.data$arom)

DimPlot(radial)+
  facet_wrap(~radial@meta.data$arom)
#so maybe youd argue 1_1 is the basal and 1_3 is a step down

ggplot(radial@meta.data, aes(x = sub.cluster, y= cyto, color = arom))+
  geom_boxplot()

#### ok what is the difference between 1_2 and 1_3
marks_1_2_3 <- FindMarkers(rgc, '1_2','1_3',latent.vars =  "nuc_prep_batch", test.use = 'LR')
clown_go(rownames(marks_1_2_3[marks_1_2_3$p_val_adj<0.05 & marks_1_2_3$avg_log2FC<(-2),]))%>%dotplot()
# ok so all of these processes are turned down, are they migrating or somehting 

marks_1_2_5 <- FindMarkers(rgc, '1_2','1_5',latent.vars =  "nuc_prep_batch", test.use = 'LR')
#sorcs1 may suggest neural progenitor
AverageExpression(rgc, features = 'LOC111562803') ### strongest in 1_2
# also lhx2, six6

marks_1_0_5 <- FindMarkers(rgc, '1_0','1_5',latent.vars =  "nuc_prep_batch", test.use = 'LR')
clown_go(rownames(marks_1_2_3[marks_1_0_5$p_val_adj<0.05 & marks_1_0_5$avg_log2FC>(1),]))%>%dotplot()
#so I guess not neurogenic?
clown_go(rownames(marks_1_2_3[marks_1_0_5$p_val_adj<0.05 & marks_1_0_5$avg_log2FC>(1.5),]))%>%dotplot()
clown_go(rownames(marks_1_2_3[marks_1_0_5$p_val_adj<0.05 & marks_1_0_5$avg_log2FC<(-2),]))%>%dotplot()

###DEGs ----
fake_degs = FindMarkers(subset(rgc, sub.cluster=='1_3'), ident.1 = 'F', group.by = 'Status')
fake_degs2 = FindMarkers(subset(rgc, sub.cluster=='1_3'), ident.1 = 'D', group.by = 'Status')
fake_degs3 = FindMarkers(subset(rgc, sub.cluster=='1_3'), ident.1 = 'M', group.by = 'Status')

fake_degs4 = FindMarkers(rgc[,rgc@meta.data$sub.cluster=='1_3' & rgc@assays$RNA$data['LOC111577263',]>0], ident.1 = 'D', group.by = 'Status')
fake_degs5 = FindMarkers(rgc[,rgc@meta.data$sub.cluster=='1_3' & rgc@assays$RNA$data['LOC111577263',]>0], ident.1 = 'F', group.by = 'Status')
fake_degs6 = FindMarkers(rgc[,rgc@meta.data$sub.cluster=='1_3' & rgc@assays$RNA$data['LOC111577263',]>0], ident.1 = 'M', group.by = 'Status')

#real degs
real_degs = read.csv( 'Subclustering/degs_1_defined_09_12_2025.csv')

# are they enriched for anything
clown_go(real_degs$gene)%>%dotplot()
clown_go(real_degs$gene[real_degs$second_word =='Upregulated'])%>%dotplot()
clown_go(real_degs$gene[real_degs$second_word =='Downregulated'])%>%dotplot()


plot_real = real_degs%>%
  group_by(cluster, short_label)%>%
  summarize(n_degs = n())

ggplot(plot_real, aes(x = cluster, y=n_degs, fill = short_label))+
  geom_bar(stat= 'identity', position ='stack')

plot_real2 = real_degs%>%
  group_by(cluster, first_word)%>%
  summarize(n_degs = n())

ggplot(plot_real2, aes(x = cluster, y=n_degs, fill = first_word))+
  geom_bar(stat= 'identity', position ='stack')

  real_degs$arom = ifelse(real_degs$gene == 'LOC111577263', 'arom', 'other')

  plot_real3 = real_degs%>%
  group_by(cluster, short_label, arom)%>%
  summarize(n_degs = n())

ggplot(plot_real3, aes(x = cluster, y=n_degs, fill = short_label, color = arom))+
  geom_bar(stat= 'identity', position ='stack', linewidth = 2)+
  scale_color_manual(values = c('black', 'gray'))
#I wonder if arom has organizing effects on 1_3 or it is a module with the other gees

mecp(rgc, 'LOC111577263', '1_0', 'sub.cluster')
mecp(rgc, 'LOC111577263', '1_1', 'sub.cluster')
mecp(rgc, 'LOC111577263', '1_3', 'sub.cluster')
mecp(rgc, 'LOC111577263', '1_4', 'sub.cluster') ### 1_4 gains aromatase basically absent prior to sex change

AverageExpression(rgc, feature = 'LOC111577263')
 #               g1-1     g1-0     g1-2     g1-5     g1-3     g1-4
#LOC111577263 2.114407 17.19239 28.98238 3.073186 22.02101 2.830006
### highest expression in 1_2 and 1_3

mecp(rgc, 'LOC111577263', '1_2', 'sub.cluster') # interesting that the primary cluster it's expressed in is not a degs


### is aromatase DEG due to a global reduction in mRNA txn ----
rgc$arom = (rgc@assays$RNA$data['LOC111577263',])

arom_means = rgc@meta.data%>%
  subset(sub.cluster=='1_3')%>%
  group_by(Status)%>%
    summarize(mean_mean = mean(colMeans))

  rgc$colMeans = colMeans(rgc@assays$RNA$data[1:1000,])
ggplot(subset(rgc@meta.data,sub.cluster == '1_3'), aes(x = Status, y = colMeans))+
  geom_boxplot()+
  geom_point(data = arom_means, aes(x = Status, y = mean_mean))
#ok no its normal

### is there an empirical p-value for the fold change in arom vs every other gene -----
Idents(rgc)<- 'sub.cluster'
fold_dom_v_mf <- function(gene, cluster){
  # Create data frame with gene expression values
  rgc_data = data.frame(Status = rgc@meta.data$Status, 
                        sub.cluster = rgc@meta.data$sub.cluster,
                        gene_expression = as.numeric(rgc@assays$RNA$data[gene,]))
  
  # Filter for specific cluster and calculate means by Status
  av <- rgc_data %>%
    subset(sub.cluster == cluster) %>%
    group_by(Status) %>%
    summarize(mean_ex = mean(gene_expression, na.rm = TRUE), .groups = 'drop')
  
  # Calculate fold change: D vs mean of M and F
  fold <- av$mean_ex[av$Status == 'D'] / mean(av$mean_ex[av$Status %in% c('M','F')], na.rm = TRUE)
  
  return(fold)
}
empirical_p = data.frame()
for(marker in unique(rgc_markers$gene)){
  print(marker)
    fold =fold_dom_v_mf(marker, '1_3')
    newd =data.frame(fold = fold, 
                     gene = marker)
    empirical_p = rbind(empirical_p, newd)
}
hist(empirical_p$fold, 100)
ggplot(empirical_p, aes(x = fold))+
  geom_histogram()+
  geom_vline(xintercept =empirical_p$fold[empirical_p$gene== 'LOC111577263'] ,color ='red')

empirical_p2 <- subset(empirical_p,!is.na(empirical_p$fold) & empirical_p$fold!=Inf & empirical_p$fold!=0)

sum(na.omit(empirical_p2$fold<empirical_p2$fold[empirical_p2$gene== 'LOC111577263']))/nrow(empirical_p2)
#p.value = 0.09

#clown_go(empirical_p2$gene[empirical_p2$fold<0.6066228])%>%dotplot()
clown_go(empirical_p2$gene[empirical_p2$fold<1])%>%dotplot()
clown_go(empirical_p2$gene[empirical_p2$fold>1])%>%dotplot()

### jaccard index -----

all_markers = FindAllMarkers(obj)
good_markers = subset(all_markers, p_val_adj < 0.05)

bored = data.frame()
for(j in 0:5){
  print(j)
for(i in 0:25){
rgc_mat = as.matrix((rgc@assays$RNA$data[,rgc$sub.cluster == paste0('1_',j)]>0))[which(rownames(rgc@assays$RNA$data)%in%good_markers$gene),1:200]  
mat_6 =as.matrix((obj@assays$RNA$data[,obj$final_clusters ==i]>0))[which(rownames(rgc@assays$RNA$data)%in%good_markers$gene),1:200] 
rgc_mat = ifelse(rgc_mat==T, 1, 0)
mat_6 = ifelse(mat_6==T, 1, 0)

d <-Jaccard(rgc_mat, mat_6 )

score =  mean(d[!is.nan(d)])
bored = rbind(data.frame(background = j, 
                         cluster = i, 
                         score = score), bored)
}
}

#1_2 is closest to 22 which is ependymal but so is everything I think, next best is 24 which seems right to me

### which subclusters are most similar

sub_bored = data.frame()
for(j in 0:5){
  print(j)
for(i in 0:5){
  if(i ==j){next}
rgc_mat = as.matrix((rgc@assays$RNA$data[,rgc$sub.cluster == paste0('1_',j)]>0))[which(rownames(rgc@assays$RNA$data)%in%rgc_markers$gene),1:200]  
rgc_mat2 = as.matrix((rgc@assays$RNA$data[,rgc$sub.cluster == paste0('1_',i)]>0))[which(rownames(rgc@assays$RNA$data)%in%rgc_markers$gene),1:200]  
rgc_mat = ifelse(rgc_mat==T, 1, 0)
rgc_mat2 = ifelse(rgc_mat2==T, 1, 0)

d <-Jaccard(rgc_mat, rgc_mat2 )

score =  mean(d[!is.nan(d)])
sub_bored = rbind(data.frame(background = j, 
                         cluster = i, 
                         score = score), sub_bored)
}
}

### markers 24

DotPlot(rgc, c(
                'sox2',
                'pax6b',
               'dclk1a',
               'dclk2a',
               'stmn2a',
               'LOC111585510', #hes1
               'elavl3',
               'gap43',
               'LOC111575074', #fabp7b
               'pax7a',
               'efnb1',
               'gfap',
               'pcna',
               's100b',
               'slc1a3a',
               'LOC111570141', #glutamine synthetase 
               'crocc2',
               'shha',
               'wnt8b',
               'sfrp1a',
               'LOC111582800', #cdnk1c
               'mki67',
               'LOC111574306', #msi1
               'LOC111562473', #bcl2
               'bcl2l1',
               'LOC111581681', #baxa
               'LOC111562459', #bax like
               'LOC111577263',#arom
               'LOC111579695' #cyp17a1
               ))+
  coord_flip()

DotPlot(obj, c('LOC111563640',
               'LOC111579695'), group.by = 'final_clusters') #strongly expressed in
#> ependymal cells, ach+ kiss+ cells, 6, 19 (also POA)
mecp(obj,'LOC111579695',22 )
marks_22 = FindMarkers(obj, 22, only.pos = T) #interesting I need to look more at this cluster

DimPlot(obj, label = T)

DotPlot(obj, c('gper1',
        'LOC111578888'))
### is there obvious evidence for neurogenesis?? -----
DotPlot(obj, 'dclk2a')+
  coord_flip()

DotPlot(obj, c('jun',
        'fosb',
        'egr1',
        'npas4a'))

cyto_whole <- CytoTRACE(obj@assays$RNA$data%>%as.matrix())
obj$cyto= cyto_whole$CytoTRACE

ggplot(obj@meta.data, aes(x = final_clusters ,y = cyto))+
  geom_boxplot()
# converging evidence suggests 0 might be the immature population? The jaccard + cyto
#potentially also 8 and 5. I wonder if thats why 0 is mixed?, 19?

# Cellular process modules ----
### there was an interesting stategy wei et al had to come up with different 
#modules for processes and then see enrichment

BiocManager::install('biomaRt', force = T)
library(biomaRt)

human_mart <- useEnsembl(biomart = 'genes',
                         dataset = 'hsapiens_gene_ensembl')
att = listAttributes(human_mart)

#human_go <- getBM(mart = human_mart, 
#                  attributes = c('entrezgene_accession', 
#                                 'go_id',
#                                 'name_1006' 
               #   ))
#human_to_ocellaris <- read.csv("Reference/hsapiens_to_aocellaris.csv")
"
right_join(human_go, by =join_by('hsapiens_name' == 'entrezgene_accession'))%>% 
  subset(!is.na(hsapiens_name)& 
           hsapiens_name != ''&
           !is.na(aocellaris_name)&
           aocellaris_name != '')


cell_cycle
"


##### how many positive marker genes does 1_1 have in common with the others ----
rgc_markers_good

rgc_markers_good_1 = subset(rgc_markers_good, cluster=='1_1')
  n_markers_1 = nrow(rgc_markers_good_1)
  
overlaps = data.frame()
for(i in 0:5){
  if(i ==1){next}
  markers = subset(rgc_markers_good, cluster == paste0('1_',i))
  
  n_markers = nrow(markers)
  
  markers_in_1 = sum(markers$gene %in% rgc_markers_good_1$gene)
  
  perc_in_1 = markers_in_1/n_markers
  
  perc_in_query = sum(rgc_markers_good_1$gene %in% markers$gene)/n_markers_1
  
  newd = data.frame(cluster =i,
                    n_markers = n_markers,
                    perc_in_1 = perc_in_1,
                    perc_1_in_query = perc_in_query)
  
  overlaps= rbind(overlaps, newd)
  
}

similarityScore = function(queryCluster){
  
    query_markers = subset(rgc_markers_good, cluster == queryCluster)
    n_query_markers = nrow(query_markers)
    
    out_data = data.frame()
    for(i in 0:5){
      testCluster = paste0('1_',i)
      if(testCluster == queryCluster){next}
      testMarkers = subset(rgc_markers_good, cluster == testCluster)
      n_test_markers = nrow(testMarkers)
      
      query_in_test = sum(query_markers$gene %in% testMarkers$gene)
      
      test_in_query = sum(testMarkers$gene %in% query_markers$gene)
      
      percent_query_in_test =query_in_test / n_query_markers
      percent_test_in_query = test_in_query/n_test_markers
      
      newd = data.frame(queryCluster = queryCluster,
                        testCluster = testCluster,
                        n_query_markers = n_query_markers,
                        n_test_markers = n_test_markers, 
                        percent_query_in_test = percent_query_in_test,
                        percent_test_in_query = percent_test_in_query)
      out_data = rbind(out_data, newd)
            
    }
return(out_data)
  
}

similarityScore('1_5')
similarityScore('1_3')


##### network analysis of marker gene overlap ----
# Required libraries
library(dplyr)
library(ggplot2)
library(igraph)
library(tidyr)

# Your existing function (assuming it's already defined)
# similarityScore = function(queryCluster){...}

# Generate similarity matrix using your function
clusters <- paste0('1_', 0:5)
all_similarities <- data.frame()

for(cluster in clusters) {
  sim_result <- similarityScore(cluster)
  all_similarities <- rbind(all_similarities, sim_result)
}

# Create symmetric similarity matrix
# First, create a complete pairwise dataset
complete_pairs <- expand.grid(queryCluster = clusters, testCluster = clusters, 
                             stringsAsFactors = FALSE)

# Add similarities from your data
similarity_data <- complete_pairs %>%
  left_join(
    all_similarities %>% 
      select(queryCluster, testCluster, percent_query_in_test, percent_test_in_query),
    by = c("queryCluster", "testCluster")
  ) %>%
  # For diagonal (self-similarity), set to 1
  mutate(
    similarity = case_when(
      queryCluster == testCluster ~ 1,
      !is.na(percent_query_in_test) ~ (percent_query_in_test + percent_test_in_query) / 2,
      TRUE ~ 0  # This shouldn't happen with complete data
    )
  )

# Convert to wide format (matrix)
similarity_matrix <- similarity_data %>%
  select(queryCluster, testCluster, similarity) %>%
  pivot_wider(names_from = testCluster, values_from = similarity, values_fill = list(similarity = 0))

# Convert to matrix format
sim_matrix <- as.matrix(similarity_matrix[, -1])
rownames(sim_matrix) <- similarity_matrix$queryCluster
colnames(sim_matrix) <- colnames(similarity_matrix)[-1]

# Heatmap visualization
library(pheatmap)
pheatmap(sim_matrix, 
         cluster_rows = TRUE, 
         cluster_cols = TRUE,
         display_numbers = TRUE,
         number_format = "%.2f",
         main = "Cluster Similarity Matrix\n(Average Bidirectional Overlap)")

# Network analysis
# Set threshold for edges (adjust as needed)
threshold <- 0.1  # 10% similarity threshold

# Create adjacency matrix
adj_matrix <- sim_matrix
adj_matrix[adj_matrix < threshold] <- 0
diag(adj_matrix) <- 0  # Remove self-loops

# Create igraph object
g <- graph_from_adjacency_matrix(adj_matrix, 
                                mode = "undirected", 
                                weighted = TRUE, 
                                diag = FALSE)

# Add node attributes (number of markers)
marker_counts <- all_similarities %>%
  select(queryCluster, n_query_markers) %>%
  distinct() %>%
  arrange(queryCluster)

V(g)$marker_count <- marker_counts$n_query_markers[match(V(g)$name, marker_counts$queryCluster)]

# Plot network
par(mfrow = c(1, 2))

# Layout 1: Force-directed
plot(g, 
     vertex.size = sqrt(V(g)$marker_count) / 10,  # Size by marker count
     vertex.label = V(g)$name,
     vertex.label.cex = 0.8,
     vertex.color = "lightblue",
     edge.width = E(g)$weight * 10,  # Edge width by similarity
     edge.label = round(E(g)$weight, 2),
     edge.label.cex = 0.6,
     layout = layout_with_fr(g),
     main = "Similarity Network\n(Force-directed layout)")

# Layout 2: Hierarchical
plot(g, 
     vertex.size = sqrt(V(g)$marker_count) / 10,
     vertex.label = V(g)$name,
     vertex.label.cex = 0.8,
     vertex.color = "lightcoral",
     edge.width = E(g)$weight * 10,
     edge.label = round(E(g)$weight, 2),
     edge.label.cex = 0.6,
     layout = layout_as_tree(g, root = which(V(g)$name == "1_1")),
     main = "Similarity Network\n(Hierarchical layout, root = 1_1)")

# Print network statistics
cat("Network Statistics:\n")
cat("Number of edges:", ecount(g), "\n")
cat("Network density:", edge_density(g), "\n")
cat("Clusters by marker count:\n")
print(marker_counts[order(-marker_counts$n_query_markers), ])

# Alternative: Minimum spanning tree for clearest lineage structure
mst <- mst(g)
plot(mst,
     vertex.size = sqrt(V(mst)$marker_count) / 8,
     vertex.label = V(mst)$name,
     vertex.label.cex = 1,
     vertex.color = "lightgreen",
     edge.width = E(mst)$weight * 15,
     edge.label = round(E(mst)$weight, 2),
     edge.label.cex = 0.8,
     layout = layout_as_tree(mst, root = which(V(mst)$name == "1_1")),
     main = "Minimum Spanning Tree\n(Clearest lineage relationships)")

par(mfrow = c(1, 1))

# Print detailed similarity table for interpretation
cat("\nDetailed similarity patterns:\n")
detailed_sim <- all_similarities %>%
  arrange(queryCluster, -percent_query_in_test) %>%
  mutate(
    relationship = case_when(
      percent_query_in_test > 0.3 & percent_test_in_query > 0.3 ~ "Sister clusters",
      percent_query_in_test > 0.1 & percent_test_in_query < 0.1 ~ "Query → Test differentiation",
      percent_query_in_test < 0.1 & percent_test_in_query > 0.1 ~ "Test → Query differentiation",
      TRUE ~ "Distant relationship"
    )
  )

print(detailed_sim %>% select(queryCluster, testCluster, percent_query_in_test, 
                             percent_test_in_query, relationship))



###good markers 1_4
marks_4 = FindMarkers(rgc, '1_4')
clown_go(rownames(marks_4[marks_4$p_val_adj<0.01 &marks_4$avg_log2FC>1, ]))%>%dotplot()

marks_0 = FindMarkers(rgc, '1_0')
clown_go(rownames(marks_0[marks_0$p_val_adj<0.01 &marks_0$avg_log2FC>1, ]))%>%dotplot()


marks_4_0  = FindMarkers(rgc, '1_4','1_0' )

DotPlot(rgc, c('lmx1a',
               'lmx1al',
               'LOC111582483',#1ba
               'lmx1bb'
               ))

#### examine rgc degs ----
geneNamer = function(gene){
  names = read.csv('Reference/genes updated.csv')
  
  name = names$NIH_description[names$NIH_accession==gene][1]
  return(name)
}

real_degs$name =c(sapply(X=real_degs$gene, FUN = geneNamer))
sparse =(real_degs[,c('gene', 'cluster', 'short_label', 'name')])

DotPlot(rgc, c('nkx2.1',
               'LOC111587547' #dbx1b like
               ,'reln',
               'olig2',
               'htr3a'
               ))



####fucking around with hub genes -----
ggplot(data = rgc@meta.data, aes(x = rgc@assays$RNA$data['LOC111577263',],
                                y = rgc@assays$RNA$data['LOC111577260',]))+
  geom_point(aes(color = Status), alpha = 0.2)+
  geom_smooth(aes(color = Status),method='lm', se = F)+
  xlim(c(0.1,6))+
  ylim(c(0.1,5))


mod = cor(rgc@assays$RNA$data['LOC111577260',][rgc@assays$RNA$data['LOC111577260',]>0 & rgc@assays$RNA$data['LOC111577263',]>0],
    rgc@assays$RNA$data['LOC111577263',][rgc@assays$RNA$data['LOC111577260',]>0 & rgc@assays$RNA$data['LOC111577263',]>0])
mod

red_genes= hub_df$gene_name[hub_df$module=='red']
red_exp = rgc@assays$RNA$data[red_genes,]%>%as.matrix()%>%t()
colnames(red_exp) = rownames( rgc@assays$RNA$data[red_genes,])
cor_matrix = cor(red_exp)

library(pheatmap)
heatmap(cor_matrix)

all_red_genes = modules$gene_name[modules$module=='red']
red_exp_all = rgc@assays$RNA$data[all_red_genes,]%>%as.matrix()%>%t()
colnames(red_exp_all) = rownames( rgc@assays$RNA$data[all_red_genes,])
cor_matrix_big = cor(red_exp_all)
diag(cor_matrix_big) <-0
heatmap(cor_matrix_big)

view(cor_matrix_big[,'LOC111577263'])

hist(cor_matrix_big)
range(cor_matrix_big)

which(cor_matrix_big == max(cor_matrix_big), arr.ind = TRUE)

means = colMeans(cor_matrix_big) # this is how you get the hub genes


