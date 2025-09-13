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

ggplot(cells_total, aes(x = sub.cluster, y = prop, color = Status))+
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

### difference between 3 and 1 ----
  clown_go = readRDS("Functions/clown_go2")  

markers_1_3 <- FindMarkers(rgc, '1_3', '1_1')
clown_go(rownames(markers_1_3[markers_1_3$p_val_adj<0.001 & markers_1_3$pct.1>markers_1_3$pct.2 ,]))%>%
  dotplot()

clown_go(rownames(markers_1_3[markers_1_3$p_val_adj<0.001 & markers_1_3$pct.1<markers_1_3$pct.2 ,]))%>%
  dotplot()

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

cyto_1_3 = lmer(cyto~Status+(1|individual), data = subset(rgc@meta.data, sub.cluster =='1_3'))
car::Anova(cyto_1_3, 3)
#nope

cyto_1_5 = lmer(cyto~Status+(1|individual), data = subset(rgc@meta.data, sub.cluster =='1_5'))
car::Anova(cyto_1_5, 3)
#YUP
pairs(emmeans(cyto_1_5, 'Status'), adjust ='none')
# but no intereting diffs

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

####proportion of cells expressing aromatase in a cluster----

DotPlot(radial, "LOC111577263")+
  coord_flip()
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
