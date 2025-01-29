#WGCNA
#https://smorabit.github.io/hdWGCNA/articles/basic_tutorial.html
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
  library(glmGamPoi)
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
  library(CytoTRACE)
  library(ggrepel)
  library(tidyverse)
  library(cowplot)
  library(patchwork)
  library(Polychrome)
P40 <- createPalette(40, c("#FF0000", "#00FF00", "#0000FF"), range = c(30, 80))
swatch(P40)
names(P40) <- NULL

mean_expression_cluster_plot<- readRDS('Functions/mean_expression_cluster_plot.rds')
prop_cluster_plot<- readRDS( 'Functions/prop_cluster_plot.rds')
define_degs_prop<- readRDS('Functions/define_degs_prop.rds')
mean_expression_cluster_data<- readRDS('Functions/mean_expression_cluster_data.rds')
prop_deg_function.rds<- readRDS('Functions/DEG_functions/prop_deg_function.rds')
define_behavior_degs<- readRDS('Functions/define_behavior_degs')
clown_go<- readRDS('Functions/clown_go')

}


#### God the support for this package is awful I have to install all this stuff manually###
library(devtools)
BiocManager::install('WGCNA', force = T)
library(WGCNA)

BiocManager::install()
BiocManager::install(c('UCell',
                       'GeneOverlap'))
devtools::install_github('wjawaid/enrichR', force = T)
devtools::install_github('smorabit/hdWGCNA', ref='dev')

library(hdWGCNA)
###FINALLY I THINK IT WORKED FUCK YOU MR SMORTABIT

####-- Tutorial --####
seurat_obj <- readRDS('/Users/ggraham/Desktop/snRNA-seq R Files 122524/RNA Object.rds')

#Yup thats my object
DimPlot(seurat_obj, group.by='harmony.wnn_res0.4_clusters', label=TRUE) +
   umap_theme() 

### Set up seurat object
seurat_obj <- SetupForWGCNA(
  seurat_obj,
  gene_select = "fraction", # the gene selection approach
  fraction = 0.02, # fraction of cells that a gene needs to be expressed in order to be included
  ### I lowered fraction cause I want to include tacr3a which is in 0.02
  wgcna_name = "tutorial" # the name of the hdWGCNA experiment
)

### Construct metacells
seurat_obj <- MetacellsByGroups(
  seurat_obj = seurat_obj,
  group.by = c("harmony.wnn_res0.4_clusters"), # specify the columns in seurat_obj@meta.data to group by
  reduction = 'harmony_wnn.umap', # select the dimensionality reduction to perform KNN on
  k = 25, # nearest-neighbors parameter
  max_shared = 10, # maximum number of shared cells between two metacells
  ident.group = 'harmony.wnn_res0.4_clusters' # set the Idents of the metacell seurat object
)

# normalize metacell expression matrix:
seurat_obj <- NormalizeMetacells(seurat_obj)

### Coexpression Network analysis
seurat_obj <- SetDatExpr(
  seurat_obj,
  group_name = "19", # the name of the group of interest in the group.by column
  group.by='harmony.wnn_res0.4_clusters', # the metadata column containing the cell type info. This same column should have also been used in MetacellsByGroups
  assay = 'RNA', # using RNA assay
  layer = 'data' # using normalized data
)

### Select soft power threshold

seurat_obj <- TestSoftPowers(
  seurat_obj,
  networkType = 'signed' # you can also use "unsigned" or "signed hybrid"
)

# plot the results:
plot_list <- PlotSoftPowers(seurat_obj)

# assemble with patchwork
wrap_plots(plot_list, ncol=2)

power_table <- GetPowerTable(seurat_obj)
head(power_table) ###ok so the guidance is treater than or equal to 0.8, so power 9 is the closest to that

###Construct co-expression network
seurat_obj <- ConstructNetwork(
  seurat_obj,
  tom_name = '19' # name of the topoligical overlap matrix written to disk
)

PlotDendrogram(seurat_obj, main='19 hdWGCNA Dendrogram')

#inspect topological overlap matrix
TOM <- GetTOM(seurat_obj)

### Module Eigengenes and Connectivity
# need to run ScaleData first or else harmony throws an error:
#seurat_obj <- ScaleData(seurat_obj, features=VariableFeatures(seurat_obj))

# compute all MEs in the full single-cell dataset
seurat_obj <- ModuleEigengenes(
 seurat_obj,
 group.by.vars="harmony.wnn_res0.4_clusters"
)

# harmonized module eigengenes:
hMEs <- GetMEs(seurat_obj)

# module eigengenes:
MEs <- GetMEs(seurat_obj, harmonized=FALSE)

# compute eigengene-based connectivity (kME):
seurat_obj <- ModuleConnectivity(
  seurat_obj,
  group.by = 'harmony.wnn_res0.4_clusters', group_name = '19'
)

# rename the modules
seurat_obj <- ResetModuleNames(
  seurat_obj,
  new_name = "19-M"
)

# plot genes ranked by kME for each module
p <- PlotKMEs(seurat_obj, ncol=5)

p


# get the module assignment table:
modules <- GetModules(seurat_obj) %>% subset(module != 'grey')

#saveRDS(modules, '~/Desktop/multiome_poa/WGCNA/Outputs/modules_012825.rds')

# show the first 6 columns:
head(modules[,1:6])

### Get hub genes
hub_df <- GetHubGenes(seurat_obj, n_hubs = 10)
#saveRDS(hub_df, '~/Desktop/multiome_poa/WGCNA/Outputs/hub_genes_012825.rds')

head(hub_df)

#saveRDS(seurat_obj, file='~/Desktop/snRNA-seq R Files 122524/hdWGCNA_object.rds')
#seurat_obj <- readRDS('~/Desktop/snRNA-seq R Files 122524/hdWGCNA_object.rds')

#### Compute hub signature score
seurat_obj <- ModuleExprScore(
  seurat_obj,
  n_genes = 25,
  method='UCell'
)

# make a featureplot of hMEs for each module
plot_list <- ModuleFeaturePlot(
  seurat_obj,
  reduction = 'harmony_wnn.umap',
  features='hMEs', # plot the hMEs
  order=TRUE # order so the points with highest hMEs are on top
)

# stitch together with patchwork
wrap_plots(plot_list, ncol=6)

subset_19_obj <-   subset(seurat_obj,harmony.wnn_res0.4_clusters==19) 
#this doesnt rly work tbh
ModuleRadarPlot(
  subset_19_obj,
  barcodes =  rownames(subset_19_obj@meta.data),
  axis.label.size=4,
  grid.label.size=4
)

###Modules <- 
MEs <- GetMEs(seurat_obj, harmonized=TRUE)
modules <- GetModules(seurat_obj)
mods <- levels(modules$module); mods <- mods[mods != 'grey']

# add hMEs to Seurat meta-data:
seurat_obj@meta.data <- cbind(seurat_obj@meta.data, MEs)

DotPlot(seurat_obj, features=mods, group.by = 'harmony.wnn_res0.4_clusters')+
  coord_flip()

####module 19-m7 contains tacr3a
genes_19_7 <- modules$gene_name[modules$module=='19-M7']
clown_go(genes_19_7)%>%dotplot()

module_go_function <- function(module_name){
  genes <-  modules$gene_name[modules$module==module_name]
 return(clown_go(genes)%>%dotplot())
  }
module_go_function('19-M3')
module_go_function('19-M19')
module_go_function('19-M2')
module_go_function('19-M7')
module_go_function('19-M18')


modules%>%
  group_by(module)%>%
  subset(module!='grey')%>%
  summarize(n = n())%>%
  plot()

#m7 is the largest by far, 

### More plotting ####
seurat_obj@misc$active_wgcna
seurat_obj@misc[["tutorial"]][["wgcna_net"]][["TOMFiles"]] <-'/Users/ggraham/Desktop/snRNA-seq R Files 122524/TOM/19_TOM.rda'

ModuleNetworkPlot(
  seurat_obj,
  outdir = '/Users/ggraham/Desktop/snRNA-seq R Files 122524/TOM/',
   n_inner = 20, # number of genes in inner ring
    n_outer = 30, # number of genes in outer ring
    n_conns = Inf, # show all of the connections
    plot_size=c(10,10), # larger plotting area
    vertex.label.cex=1 # font size
)

options('future.globals.maxSize'=Inf)
# hubgene network
HubGeneNetworkPlot(
  seurat_obj,
  n_hubs = 3, n_other=5,
  edge_prop = 0.75,
  mods = 'all',
  return_graph = T
)

# get the list of modules:
modules <- GetModules(seurat_obj)
mods <- levels(modules$module); mods <- mods[mods != 'grey']

# hubgene network
HubGeneNetworkPlot(
  seurat_obj,
  n_hubs = 10, n_other=20,
  edge_prop = 0.75,
  mods = mods[c(3,7,23)] # only select 5 modules
)

###aplying umap
seurat_obj <- RunModuleUMAP(
  seurat_obj,
  n_hubs = length(unique(modules$module))-1, # number of hub genes to include for the UMAP embedding
  n_neighbors=15, # neighbors parameter for UMAP
  min_dist=0.1 # min distance between points in UMAP space
)

umap_df <- GetModuleUMAP(seurat_obj)

ggplot(umap_df, aes(x=UMAP1, y=UMAP2)) +
  geom_point(
   color=umap_df$color, # color each point by WGCNA module
   size=umap_df$kME*2 # size of each point based on intramodular connectivity
  ) +
  umap_theme()

ModuleUMAPPlot(
  seurat_obj,
  edge.alpha=0.25,
  sample_edges=TRUE,
  edge_prop=0.1, # proportion of edges to sample (20% here)
  label_hubs=2 ,# how many hub genes to plot per module?
  keep_grey_edges=FALSE

)

#https://smorabit.github.io/hdWGCNA/articles/network_visualizations.html

#saveRDS(seurat_obj, file='~/Desktop/snRNA-seq R Files 122524/hdWGCNA_object.rds')
