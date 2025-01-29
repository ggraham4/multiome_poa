#WGCNA on POA Cells Only
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
  library(scCustomize)

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
library(hdWGCNA)
library(WGCNA)

}

### Read in RNA only obj
obj <- readRDS('/Users/ggraham/Desktop/snRNA-seq R Files 122524/RNA Object.rds')

###subset to POA only clusters via label transfer from old project
### exclude 15 cause it is missing in some
obj <- subset(obj, harmony.wnn_res0.4_clusters %in%
                c(0,
                 1,
                 8,
                 17,
                 19,
                 21,
                 23,
                 25,
                 27,
                 28,
                 31
  )
)

clusters_list <-   c(0,
                 1,
                 8,
                 17,
                 19,
                 21,
                 23,
                 25,
                 27,
                 28,
                 31
  )
#confirm only the ones we want are in there
unique(obj$harmony.wnn_res0.4_clusters)
DimPlot(obj)

###Find expression level of tacr3a in this group, I want to make sure it is included
tac_expression <- scCustomize::Percent_Expressing(obj, 'tacr3a')%>%
  t()
mean(tac_expression)
##set cutoff for 3%

### set up for WGCNA
obj <- SetupForWGCNA(
  obj,
  gene_select = "fraction",
  fraction = 0.03, 
  wgcna_name = "putative_poa_cells" 
)

### Construct metacells
obj <- MetacellsByGroups(
  seurat_obj = obj,
  group.by = c("harmony.wnn_res0.4_clusters"),
  reduction = 'harmony_wnn.umap', 
  k = 25, 
  max_shared = 10, # again, not sure if these should be altered or not
  ident.group = 'harmony.wnn_res0.4_clusters'
)

# normalize metacell expression matrix:
obj <- NormalizeMetacells(obj)

### Coexpression Network analysis
obj <- SetDatExpr(
  obj,
  group_name = paste0(clusters_list),
  #changing this to be all clusters
  group.by='harmony.wnn_res0.4_clusters', 
  assay = 'RNA', 
  layer = 'data' 
)

### Select soft power threshold
obj <- TestSoftPowers(
  obj,
  networkType = 'signed' #not sure if this should be altered either
)
# plot the results:
plot_list <- PlotSoftPowers(obj)
wrap_plots(plot_list, ncol=2)

#guidance is to use greater than or equal to 0.8, which corresponds to
plot_list[[1]][["data"]]
#a power of 6

###Construct co-expression network
obj <- ConstructNetwork(
  obj,
  tom_name = 'putative_poa_cells', # name of the topoligical overlap matrix written to disk
tom_outdir = '/Users/ggraham/Desktop/snRNA-seq R Files 122524/TOM/',
overwrite_tom = T
  )

PlotDendrogram(obj, main='putative_poa_cells hdWGCNA Dendrogram')

#model eigengenes
obj <- ModuleEigengenes(
 obj,
 group.by.vars="harmony.wnn_res0.4_clusters"
)


# harmonized module eigengenes:
hMEs <- GetMEs(obj)
# module eigengenes:
MEs <- GetMEs(obj, harmonized=FALSE)
### have to do this because the default path is changed
obj@misc[["putative_poa_cells"]][["wgcna_net"]][["TOMFiles"]] <-'/Users/ggraham/Desktop/snRNA-seq R Files 122524/TOM/putative_poa_cells_TOM.rda'

# compute eigengene-based connectivity (kME):
obj <- ModuleConnectivity(
  obj,
  group.by = 'harmony.wnn_res0.4_clusters', group_name =  paste0(clusters_list)
  
)

#Ok turquoise has oligo markers, 
p <- PlotKMEs(obj, ncol=5)
p

# get the module assignment table:
modules <- GetModules(obj) %>% subset(module != 'grey')

#saveRDS(modules, '~/Desktop/multiome_poa/WGCNA/Outputs/modules_putative_poa_cells_012925.rds')

### Get hub genes
hub_df <- GetHubGenes(obj, n_hubs = 10)
#saveRDS(hub_df, '~/Desktop/multiome_poa/WGCNA/Outputs/hub_genes_putative_poa_cells_012925.rds')

head(hub_df)

#### Compute hub signature score
obj <- ModuleExprScore(
  obj,
  n_genes = 25,
  method='UCell'
)

# make a featureplot of hMEs for each module
plot_list <- ModuleFeaturePlot(
  obj,
  reduction = 'harmony_wnn.umap',
  features='hMEs', # plot the hMEs
  order=TRUE 
)
wrap_plots(plot_list, ncol=6)

ModuleRadarPlot(
  obj,
  barcodes =  obj@meta.data %>% rownames(),
  axis.label.size=4,
  grid.label.size=4
)

###Modules <- 
MEs <- GetMEs(obj, harmonized=TRUE)
modules <- GetModules(obj)
mods <- levels(modules$module); mods <- mods[mods != 'grey']

obj@meta.data <- cbind(obj@meta.data, MEs)

DotPlot(obj, features=mods, group.by = 'harmony.wnn_res0.4_clusters')+
  coord_flip()

module_go_function <- function(module_name){
  genes <-  modules$gene_name[modules$module==module_name]
 return(clown_go(genes)%>%dotplot())
  }
module_go_function('brown')
module_go_function('turquoise')
module_go_function('red')
module_go_function('yellow')
module_go_function('black')
module_go_function('green')
module_go_function('magenta')
module_go_function('purple')
module_go_function('blue')


modules%>%
  group_by(module)%>%
  subset(module!='grey')%>%
  summarize(n = n())%>%
  plot()

ModuleNetworkPlot(
  obj,
  outdir = '/Users/ggraham/Desktop/snRNA-seq R Files 122524/TOM/all_cells/',
   n_inner = 20, # number of genes in inner ring
    n_outer = 30, # number of genes in outer ring
    n_conns = Inf, # show all of the connections
    plot_size=c(10,10), # larger plotting area
    vertex.label.cex=1 # font size
)

options('future.globals.maxSize'=Inf)

### still doesnt work
HubGeneNetworkPlot(
  obj,
  n_hubs = 3, n_other=5,
  edge_prop = 0.75,
  mods = 'all'
)

###aplying umap
obj <- RunModuleUMAP(
  obj,
  n_hubs = length(unique(modules$module))-1, 
  n_neighbors=15, 
  min_dist=0.1 
)
umap_df <- GetModuleUMAP(obj)

ggplot(umap_df, aes(x=UMAP1, y=UMAP2)) +
  geom_point(
   color=umap_df$color, # color each point by WGCNA module
   size=umap_df$kME*2 # size of each point based on intramodular connectivity
  ) +
  umap_theme()

ModuleUMAPPlot(
  obj,
  edge.alpha=0.25,
  sample_edges=TRUE,
  edge_prop=0.1, 
  label_hubs=2 ,
  keep_grey_edges=FALSE,
    return_graph = F
)
#saveRDS(obj, file='~/Desktop/snRNA-seq R Files 122524/hdWGCNA_object_putative_poa_cells.rds')









