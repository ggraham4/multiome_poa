{
  #> OK FUTURE GABE HERE IS THE LORE
  #> This prick rahul satija messed with GetAssayData
  #> to do who the fuck knows what, but anyway it changed an argument from 
  #> slot to layer, smorabito of hdWGCNA fame fixed it, but not for the most
  #> recent versions at time of writing 01/18/2026, 
  #> SO, you have to force install these versions below and make sure you 
  #> dont update anything to get it to work. If I ever meet him there
  #> will be violence
  
  # Rahul satija is an asshole 
  remotes::install_version("Seurat", version = "5.0.3")
  remotes::install_version("SeuratObject", version = "5.0.1") 
  
  library(SeuratObject)
  library(Seurat)
  
  packageVersion("Seurat")
  packageVersion("SeuratObject")
  
  library(hdWGCNA)
  
  library(patchwork)
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
}
#read in data
obj = readRDS("~/Desktop/optimal_clustering_rna_only.rds")


####  subcluster and subset ----
Idents(obj) = 'final_clusters'
rgc <- FindSubCluster(obj, 1, 'harmony.wsnn')
rgc = subset(rgc, final_clusters == 1)

Idents(rgc) <- 'sub.cluster'

DimPlot(rgc, label = T, reduction = 'harmony_wnn.umap')

rgc$Status = factor(rgc$Status, levels = c('NRM', 'M','D','E','NF','F'))


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

#radar plot
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

### OK we are finally up to date, I am going to save this seurat object
# so I dont have to deal with any more dependency stuff ever again
#saveRDS(rgc, '/Users/ggraham/WGCNA_rgc/rgc_seurat_object.rds')
rgc = readRDS('/Users/ggraham/WGCNA_rgc/rgc_seurat_object.rds')
