#subclustering microglia
#> GJG 2025 09 08 

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

calculateCoexpressedGenes = function(seurat_obj2, geneList, clustering = 'sub.cluster'){
  ### Find genes coexpressed with a gene list by splitting cells into positive/negative
  ### based on expression of ANY gene in the list
  
  `%notin%` = Negate(`%in%`)
  seurat_obj <- seurat_obj2
  seurat_obj$clustering = seurat_obj[[clustering]]
  
  # Filter gene list to only include genes present in the data
  available_genes <- geneList[geneList %in% rownames(seurat_obj@assays$RNA$data)]
  
  if(length(available_genes) == 0) {
    stop("None of the genes in geneList are found in the Seurat object")
  }
  
  message("Using ", length(available_genes), " out of ", length(geneList), " genes from the input list")
  if(length(available_genes) < length(geneList)) {
    missing_genes <- geneList[geneList %notin% rownames(seurat_obj@assays$RNA$data)]
    message("Missing genes: ", paste(missing_genes, collapse = ", "))
  }
  
  # Initialize all cells as negative
  seurat_obj$positive <- FALSE
  
  # Set cells to positive if ANY of the available genes is expressed (> 0)
  for(gene in available_genes) {
    gene_expression <- seurat_obj@assays$RNA$data[gene, ]
    seurat_obj$positive <- seurat_obj$positive | (gene_expression > 0)
  }
  
  marker_all <- data.frame()
  clusters_processed <- 0
  
  for(cluster in unique(seurat_obj$clustering)) {
    temp_obj <- subset(seurat_obj, clustering == cluster)
    
    # Check if we have enough positive and negative cells
    n_positive <- sum(temp_obj$positive)
    n_negative <- sum(!temp_obj$positive)
    
    if(n_positive >= length(unique(seurat_obj$clustering)) && length(unique(seurat_obj$clustering)) >= 3) {
      Idents(temp_obj) <- 'positive'
      
      tryCatch({
        markers <- FindAllMarkers(temp_obj,
                                  assay = "RNA",
                                  group.by = 'positive',
                                  logfc.threshold = 0,
                                  min.pct = 1/nrow(temp_obj@meta.data),
                                  verbose = FALSE)
        markers$sub.cluster <- cluster
        marker_all <- rbind(markers, marker_all)
        clusters_processed <- clusters_processed + 1
      }, error = function(e) {
        message("Warning: Could not find markers for cluster ", cluster, ": ", e$message)
      })
    } else {
      message("Skipping cluster ", cluster, ": insufficient cells (positive: ", 
              n_positive, ", negative: ", n_negative, ")")
    }
  }
  
  if(nrow(marker_all) == 0) {
    return(message('No markers found. Check if there are sufficient positive/negative cells in each cluster'))
  }
  
  message("Successfully processed ", clusters_processed, " clusters")
  
  # Extract significant genes (positive markers only)
  marker_all_signif <- marker_all[marker_all$p_val_adj < 0.05 & marker_all$cluster == TRUE, ]
  
  if(nrow(marker_all_signif) == 0) {
    return(message('No significant markers found'))
  }
  
  # Genes found in at least half of the clusters
  min_clusters <- ceiling(length(unique(seurat_obj$clustering)) / 2)
  marker_genes_half <- marker_all_signif %>%
    group_by(gene) %>%
    summarize(n_clusters = n(), .groups = 'drop') %>%
    filter(n_clusters >= min_clusters)
  
  new_markers <- marker_genes_half$gene[marker_genes_half$gene %notin% available_genes]
  message('Found ', length(new_markers), ' new coexpressed markers')
  
  if(length(new_markers) > 0) {
    print(new_markers)
  }
  
  # Assign to global environment
  assign('markersCoex', marker_genes_half$gene, env = .GlobalEnv)
  
  # Calculate module score - filter markers to only include genes present in object
  all_markers <- marker_genes_half$gene
  available_markers <- all_markers[all_markers %in% rownames(seurat_obj@assays$RNA$data)]
  
  if(length(available_markers) == 0) {
    warning("No markers available for module score calculation")
    seurat_obj$coexModule1 <- 0
  } else {
    message("Calculating module score with ", length(available_markers), " markers")
    if(length(available_markers) != length(all_markers)) {
      missing_markers <- all_markers[all_markers %notin% rownames(seurat_obj@assays$RNA$data)]
      message("Note: ", length(missing_markers), " markers not found in object: ", 
              paste(head(missing_markers, 5), collapse = ", "))
    }
    
    seurat_obj <- AddModuleScore(seurat_obj, 
                                features = list(available_markers),
                                name = 'coexModule')
  }
  
  # Add summary information as metadata
  seurat_obj@misc$coexpression_analysis <- list(
    input_genes = geneList,
    available_genes = available_genes,
    clusters_processed = clusters_processed,
    total_markers = length(all_markers),
    new_markers = length(new_markers)
  )
  
  return(seurat_obj)
}
}
#read in data
obj = readRDS("~/Desktop/optimal_clustering_rna_only.rds")

####  subcluster and subset ----
Idents(obj) = 'final_clusters'
microglia <- FindSubCluster(obj, 11, 'harmony.wsnn')
microglia = subset(microglia, final_clusters == 11)

#3 clusters to be expected
DimPlot(microglia, reduction = 'harmony_wnn.umap', group.by = 'sub.cluster')

FeaturePlot(obj,'ptprc', reduction = 'harmony_wnn.umap')
FeaturePlot(obj,'p2ry12', reduction = 'harmony_wnn.umap')
#they are definitely microglia so why tf are my GO results so weird below??

####  markers ----
microglia_markers = FindAllMarkers(microglia, group.by = 'sub.cluster', test.use = 'LR')

# what are they enriched for
clown_go(microglia_markers$gene[microglia_markers$cluster=='11_0' & microglia_markers$p_val_adj <0.001])%>%dotplot()
# what
clown_go(microglia_markers$gene[microglia_markers$cluster=='11_1' & microglia_markers$p_val_adj <0.001])%>%dotplot()

clown_go(microglia_markers$gene[microglia_markers$cluster=='11_2' & microglia_markers$p_val_adj <0.001])%>%dotplot()
#this dendrite shit is stupid

# ok moving on because what
markers_11_0 = microglia_markers$gene[microglia_markers$cluster=='11_0']
markers_11_1 = microglia_markers$gene[microglia_markers$cluster=='11_1']
markers_11_2 = microglia_markers$gene[microglia_markers$cluster=='11_2']
#too much overlap, I am going to have to intersect them

good_markers_11_0 = markers_11_0[markers_11_0 %notin% c(markers_11_1, markers_11_2)]
good_markers_11_1 = markers_11_1[markers_11_1 %notin% c(markers_11_0, markers_11_2)]
good_markers_11_2 = markers_11_2[markers_11_2 %notin% c(markers_11_1, markers_11_0)]

#clown_go(good_markers_11_0)%>%dotplot()
#clown_go(good_markers_11_1)%>%dotplot()
#clown_go(good_markers_11_2)%>%dotplot()
#still confused 


#### CytoTRACE ----
cyto = CytoTRACE(microglia@assays$RNA$data%>%as.matrix())
microglia$cyto = cyto$CytoTRACE

ggplot(microglia@meta.data, aes(x = sub.cluster,y = cyto))+
  geom_violin()+
  stat_summary(geom ='crossbar', fun = 'median')
#11_2 is immature which I expected from the sheer number of marker genes

microglia@meta.data$Status = factor(microglia@meta.data$Status, levels = c('NRM','M','D','E','NF','F'))
ggplot(microglia@meta.data, aes(x = sub.cluster,y = cyto, fill = Status))+
  stat_summary(geom='bar', position = 'dodge', fun = 'mean')
#like you could maybe argue there is a difference between males and females in 11_2
# but i think that is a reach

t.test(microglia@meta.data$cyto[microglia@meta.data$Status=='M'], microglia@meta.data$cyto[microglia@meta.data$Status=='F'])
#nope

### proportion ----
props = microglia@meta.data%>%
  group_by(sub.cluster)%>%
  summarize(n_cells = n(), 
            total = nrow(microglia@meta.data),
            prop = n_cells/total)

ggplot(props, aes(x = sub.cluster, y = prop))+
  geom_point()
# ok so 0 is the most populous easily and 2 is almost nonexistant
props_ind = microglia@meta.data%>%
  group_by(sub.cluster, individual, Status)%>%
  summarize(n_cells = n())
props_total = microglia@meta.data%>%
  group_by( individual)%>%
  summarize(n = n())

props_ind_joined = props_ind%>%
  right_join(props_total, by = 'individual')
props_ind_joined$prop = props_ind_joined$n_cells/props_ind_joined$n

props_ind_joined$Status = factor(props_ind_joined$Status , levels = c('NRM', 'M','D','E','NF','F'))

ggplot(props_ind_joined, aes(x = sub.cluster, y = prop, color = Status))+
  geom_point(position =position_jitterdodge(1))+
  geom_boxplot(alpha = 0)
#reallly loooks like nothing here

### WGCNA ----
clusters_list <- c('11_0',
                   '11_1',
                   '11_2')

microglia <- SetupForWGCNA(
  microglia,
  gene_select = "fraction",
  fraction = 0.05,
  wgcna_name = "microglia" 
)

### Construct metacells
microglia <- MetacellsByGroups(
  seurat_obj = microglia,
  group.by = c("sub.cluster"),
  reduction = 'harmony_wnn.umap', 
  k = 25, 
  max_shared = 10,
  ident.group = 'sub.cluster',
  min_cells = 90 # to keep 2
)

# normalize metacell expression matrix:
microglia <- NormalizeMetacells(microglia)

### Coexpression Network analysis
microglia <- SetDatExpr(
  microglia,
  group_name = c('11_0', '11_1', '11_2'),
  group.by='sub.cluster', 
  assay = 'RNA', 
  layer = 'data' 
)

### Select soft power threshold
microglia <- TestSoftPowers(
  microglia,
  networkType = 'signed' #not sure if this should be altered either
)

# plot the results:
plot_list <- PlotSoftPowers(microglia)
wrap_plots(plot_list, ncol=2)

#guidance is to use greater than or equal to 0.8, which corresponds to
plot_list[[1]][["data"]]
#softpower of 8 

microglia <- ConstructNetwork(
  microglia,
  tom_name = 'microglia', 
  tom_outdir = '/Users/ggraham/WGCNA_microglia/',
  overwrite_tom = T
)

PlotDendrogram(microglia, main='clust_11 hdWGCNA Dendrogram')
# two main modules thats what we like to see

#model eigengenes
microglia <- ModuleEigengenes(
  microglia,
  group.by.vars="sub.cluster"
)

# harmonized module eigengenes:
hMEs <- GetMEs(microglia)
# module eigengenes:
MEs <- GetMEs(microglia, harmonized=FALSE)

microglia@misc[["microglia"]][["wgcna_net"]][["TOMFiles"]] <-'/Users/ggraham/WGCNA_microglia/microglia_TOM.rda'

# compute eigengene-based connectivity (kME):
microglia <- ModuleConnectivity(
  microglia,
  group.by = 'sub.cluster', group_name =  c('11_0','11_1', '11_2')
  
)

p <- PlotKMEs(microglia, ncol=5)
p

# get the module assignment table:
modules <- GetModules(microglia) %>% subset(module != 'grey')

### Get hub genes
hub_df <- GetHubGenes(microglia, n_hubs = 10)

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
  microglia,
  group.by = 'sub.cluster',
  axis.label.size=4,
  grid.label.size=4
) # 11_2 is enriched for turquoise

## still very confused

MEs$individual = microglia$individual
MEs$Status = microglia$Status
MEs$sub.cluster = microglia$sub.cluster

plot_module('blue') #maybe? a difference, though Im not convinced
plot_module('turquoise') # again hard to say
plot_module('brown') #no difference fs

turquoise_test = lmer(turquoise~Status+(1|individual), data = MEs)
car::Anova(turquoise_test, 3)

blue_test = lmer(blue~Status+(1|individual), data = MEs)
car::Anova(blue_test, 3)
#nothing

plot_module_cluster('blue') # strongly enriched for 0 id say
plot_module_cluster('turquoise') # VERY strongly enriched for 2
plot_module_cluster('brown')  # 1


clown_go(modules$gene_name[modules$color=='turquoise'])%>%dotplot() # turquoise seems like binding to neurons maybe
clown_go(modules$gene_name[modules$color=='blue'])%>%dotplot() # so blue seems like activated module
clown_go(modules$gene_name[modules$color=='brown'])%>%dotplot() # brown seems like homeostatic or proliferating

### so, it seems like 11_1 (brown) does not change by sex, whie blue could increase transiently though the stats dont support
##also, my guess is that turqoise is a marker of immature/ division though that is not supported by enrichment

## i am going to lower the threshold to get 11_2 back in here
# done, and results are the same

### are any modules enriched for degs ----
deg_data = read_csv("DEG Outputs/FINAL_degs_classified_08_11_2025.csv")
deg_data_11 = subset(deg_data, cluster == 11)

#mark the degs
modules$isDEG = ifelse(modules$gene_name %in% deg_data_11$gene, TRUE, FALSE)

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
#nope

#### what do the module members suggest the module is ----
clown_go(modules$gene_name[modules$module=='blue'])%>%dotplot()
#lysosome, phagocytosis
clown_go(hub_df$gene_name[hub_df$module=='blue'])%>%dotplot()

#blue expresses p2ry12, hexb suggesting it is  M0, but m2 markers are also present so I think blue is too all
# encompassing

clown_go(modules$gene_name[modules$module=='turquoise'])%>%dotplot()
#neurotransmission
clown_go(modules$gene_name[hub_df$module=='turquoise'])%>%dotplot()

clown_go(modules$gene_name[modules$module=='brown'])%>%dotplot()
# dna and rna txn, and also response to steroids?, and cell cycle is interesting
clown_go(modules$gene_name[hub_df$module=='brown'])%>%dotplot()

brown_go_df = data.frame(
gene_id= clown_go(modules$gene_name[modules$module=='brown'])$geneID,
desc = clown_go(modules$gene_name[modules$module=='brown'])$Description
)

turquoise_go_df = data.frame(
gene_id= clown_go(modules$gene_name[modules$module=='turquoise'])$geneID,
desc = clown_go(modules$gene_name[modules$module=='turquoise'])$Description
)


### im wondering if I can make like an index of inflammation or something for each nucleus
# in the same way zack made the ieg score
Idents(microglia) <-'sub.cluster'
DotPlot(microglia, c('p2ry12', #homeostatic
                     'apoeb', #disease
                     'lamp1a', #phagocytosis
                     'cd68'#cd68, phagycytosis marker
                     ))+
  coord_flip()

# i think I have to conclude that the number of clusters is wrong, these should be distinct populations


#### trying different clusters ----
microglia2 <- FindSubCluster(obj, 11, 'harmony.wsnn', resolution = 0.8)
microglia2 = subset(microglia2, final_clusters == 11)
Idents(microglia2) = 'sub.cluster'
DimPlot(microglia2)

DotPlot(microglia2, c('p2ry12', #homeostatic
                     'apoeb', #disease
                     'lamp1a', #phagocytosis
                     'cd68'#cd68, phagycytosis marker
                     ))+
  coord_flip()
# I dont really know whats going on----

### trying some meghan markers ----
microglia = AddModuleScore(microglia,
               list(c('p2ry12',
                 'lcp1',
                 'rgs18',
                 'hivep1',
                 'nr4a3',
                 'dusp5',
                 'marco',
                 'grapa',
                 'rac2',
                 'tnfaip8l2b')),
               ctrl = 10,
               name = 'Homeostatic_')

microglia = AddModuleScore(microglia,
               list(tolower(c('fxr2',
                 'Marcksl1a',
                 'Marcksa',
                 'Marcksb',
                 'Pxylp1'))),
               ctrl = 10,
               name = 'Dev')

microglia = AddModuleScore(microglia,
               list(tolower(c(
                 'cebpb',
                 'irf1b',
                 'gpnmb',
                 'socs3a',
                 'bhlhe40',
                 'jun',
                 'ninj1'))),
               ctrl = 10,
               name = 'Inflammatory')


DotPlot(microglia, c('Homeostatic_1', 'Dev1', 'Inflammatory1'))+
  coord_flip()



ggplot(microglia@meta.data, aes(x =Status, y = Inflammatory1))+
  geom_boxplot()


""
'
classical activation state

il1b
TNF
il6
il-12
ifn-gamma
ccl2
'
""
classical = c(
  #'LOC111583887', #il1Beta like
  #'tnfa',
 # 'tnfb',
  #'LOC118471830', #interferon gamma like
  #'LOC111565507', # ccl2 like
 'apoeb',
 'cst3',
#'spp1',
 'cd74b',
 'LOC111580055',
'rpl29'
)

classicalCoex = calculateCoexpressedGenes(microglia, classical) # none of them are expressed?
DotPlot(classicalCoex, 'coexModule1') # 0 seems to be activated?
DotPlot(classicalCoex, classical)+
  coord_flip()


""
'
alternative activation - neuroprotective
il4
il3
il4
il10
il13
tgf-beta

'
""
alternative = c(
  'il10',
  'tgfb3'
)
#altCoex = calculateCoexpressedGenes(microglia, alternative) 
FeaturePlot(microglia, 'il10') 
FeaturePlot(microglia, 'tgfb3') 
FeaturePlot(obj, 'il10') 
FeaturePlot(obj, 'tgfb3') 

#collasssal waste of time

### Genes coexpressed with grik4 ----
FeaturePlot(microglia, 'grik4') 

grikCoex = calculateCoexpressedGenes(microglia, c('grik4')) 
clown_go(markersCoex)%>%dotplot()

nosCoex = calculateCoexpressedGenes(microglia, c('nos1')) 
clown_go(markersCoex)%>%dotplot()

sema5baCoex = calculateCoexpressedGenes(microglia, c('sema5ba')) 
clown_go(markersCoex)%>%dotplot()

mecp(microglia, 'grik4', '11_0', 'sub.cluster') #0 and 1 drive the DEGness
mecp(microglia, 'grik4', '11_1', 'sub.cluster')
mecp(microglia, 'grik4', '11_2', 'sub.cluster')
mecp(microglia, 'grik4', '11', 'final_clusters')

DotPlot(microglia, 'grik4')+
  coord_flip()

## mbpl-----
mbplCoex = calculateCoexpressedGenes(microglia, c('LOC111580029')) 
clown_go(markersCoex)%>%dotplot()

mecp(microglia, 'LOC111580029', '11', 'final_clusters')
mecp(microglia, 'LOC111580029', '11_0', 'sub.cluster') 
mecp(microglia, 'LOC111580029', '11_1', 'sub.cluster')
mecp(microglia, 'LOC111580029', '11_2', 'sub.cluster')

DotPlot(microglia, 'LOC111580029')+
  coord_flip()


### DEG module ----

mg_degs = read_csv("DEG Outputs/FINAL_degs_classified_08_11_2025.csv")
mg_degs = subset(mg_degs, cluster =='11')

deg_module = calculateCoexpressedGenes(microglia, mg_degs$gene) 
DotPlot(deg_module, 'positive')+
  coord_flip() # degs are predominantly in 11_0 and 11_2

DotPlot(microglia, mg_degs$gene)+
  coord_flip() 

### meghan stress module ----
meghan_stress_genes = read.xlsx('/Users/ggraham/downloads/media-3.xlsx',3)
sig_stress_genes = meghan_stress_genes$gene[meghan_stress_genes$q.value_stress<0.05]

ensembl_ocellaris <- useEnsembl(biomart = "genes", 
                      dataset = "aocellaris_gene_ensembl")
att = listAttributes(ensembl_ocellaris)
biomart_names = 
  getBM(
    mart = ensembl_ocellaris, 
    attributes = c('external_gene_name', 
                   'entrezgene_accession'))

  biomart_homologs <-
  getBM(
    mart = ensembl_ocellaris, 
    attributes = c('external_gene_name', 
                   'mmusculus_homolog_associated_gene_name'))
  
  biomart_joined = biomart_names%>%
    subset(external_gene_name != '' & entrezgene_accession != '')%>%
    right_join(biomart_homologs, by = 'external_gene_name')



clown_sig_stress_genes = unique(
  biomart_joined$entrezgene_accession[biomart_joined$mmusculus_homolog_associated_gene_name %in% sig_stress_genes]) 

good_stress_genes = clown_sig_stress_genes[clown_sig_stress_genes %in% rownames(microglia@assays$RNA$data)]

microglia = AddModuleScore(microglia, 
                           features = list(good_stress_genes),
                           name ='Stress')
DotPlot(microglia, 'Stress1')+
  coord_flip()

sex_diff_stress =lmer(Stress1~Status +(1|individual), microglia@meta.data)
summary(sex_diff_stress)
car::Anova(sex_diff_stress, 3) #fuck

ggplot(microglia@meta.data, aes(x = Status, y = Stress1, color = sub.cluster))+
  geom_boxplot()

ggplot(microglia@meta.data, aes(x = sub.cluster, y = Stress1))+
  geom_violin()

# I wonder if the clustering needs to be different
Idents(microglia) <- 'final_clusters'

microglia = FindSubCluster(microglia, 11, 
                           'harmony.wsnn', 
                           'sub_2',
                           resolution = .8)
Idents(microglia) = 'sub_2'
DimPlot(microglia)

ggplot(microglia@meta.data, aes(x = sub_2, y = Stress1))+
  geom_violin()

Idents(microglia) = 'sub.cluster'

# I dont really think so 

table(microglia$sub.cluster, microglia$Stress1 > mean(microglia$Stress1))


