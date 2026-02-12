
{
  library(ggsignif)
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
  
  geneNamer = function(gene){
  names = read.csv('Reference/genes updated.csv')
  
  name = names$NIH_description[names$NIH_accession==gene][1]
  
  if(is.na(name)){name = gene}
  return(name)
}

  plot_gene_ieg = function(gene_name, seurat_obj = sub_6, statuses = c('M','D','F')) {
  
  # Check if gene exists in the data
  if (!gene_name %in% rownames(seurat_obj@assays$RNA$data)) {
    stop(paste("Gene", gene_name, "not found in the data"))
  }
  
  # Create binary gene expression variable
  seurat_obj@meta.data[[gene_name]] = ifelse(
    seurat_obj@assays$RNA$data[gene_name,] > 0, 
    TRUE, 
    FALSE
  )
  
  # Prepare plot data
  plot_data = seurat_obj@meta.data %>%
    group_by(individual, Status, !!sym(gene_name)) %>%
    summarize(mean_ieg = mean(ieg_PC1), .groups = 'drop') %>%
    subset(Status %in% statuses)

  # Create plot
  p = ggplot(plot_data, aes(x = Status, y = mean_ieg, color = !!sym(gene_name))) +
    geom_boxplot() +
    geom_point(position = position_jitterdodge(.5)) +
    labs(title = paste("IEG expression by", gene_name, "status"),
         y = "IEG PC1",
         color = gene_name) +
    theme_minimal()
  
  return(p)
}

  #define functions

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
  #seurat_obj = AddModuleScore(seurat_obj, 
  #                            features = list(all_markers),
  #                            name = 'coexModule')
  
  return(all_markers)
}




}
obj = readRDS("~/Desktop/optimal_clustering_rna_only.rds")

sub_6 = FindSubCluster(obj, 6, graph.name='harmony.wsnn')
Idents(sub_6) <- 'sub.cluster'
sub_6 = subset(sub_6, final_clusters ==6)
sub_6$sub.cluster = factor(sub_6$sub.cluster, levels = c(paste0('6_',0:3)))
sub_6$Status = factor(sub_6$Status, levels = c('NRM','M',"D",'E','NF','F'))

DimPlot(sub_6)

# Find IEGs ---
iegs <- c('npas4a', 'fosb','egr1') #using fosb instead of cfos , I think it should be more or less the same

iegs_= calculateCoexpressedGenes(sub_6, iegs) #c-fos like is included anyway so its good

### Define a function to return PC1 for the module
whole_matrix = sub_6@assays$RNA$data%>%as.matrix()%>%t()

module_pc1 = function(module){
  genes = colnames(whole_matrix)
  module_matrix = whole_matrix[,which(genes %in% module)]
  
  pc = prcomp(module_matrix, scale = T)
  scores <- pc$x[,1]

  if(mean(pc$rotation[c("fosb","egr1","npas4a"),1], na.rm=TRUE) < 0){
  scores <- -scores
}
return(scores)

}

sub_6$ieg_PC1 = module_pc1(iegs_)

# check to make sure it makes sense
DotPlot(sub_6, 'ieg_PC1')

ggplot(subset(sub_6@meta.data, Status != 'NRM'), aes(x = Status, y = ieg_PC1))+
  geom_boxplot()+
  facet_wrap(~sub.cluster, scales ='free')

ggplot(subset(sub_6@meta.data, Status != 'NRM'), aes(x = sub.cluster, y = ieg_PC1))+
  geom_boxplot()

DotPlot(sub_6, iegs)
# ok the PC score is correct I believe

plot_gene_ieg('LOC111571064')+
    facet_wrap(~LOC111571064)

plot_gene_ieg('cckb')+
  facet_wrap(~cckb)

plot_gene_ieg('LOC111562384')+
  facet_wrap(~LOC111562384)

plot_gene_ieg('tacr3a')+
  facet_wrap(~tacr3a)

plot_gene_ieg('drd3')+
  facet_wrap(~drd3)

plot_gene_ieg('npy7r')+
  facet_wrap(~npy7r)

plot_gene_ieg('kiss1ra')+
  facet_wrap(~kiss1ra)

plot_gene_ieg('nmbr')+
  facet_wrap(~nmbr)

plot_gene_ieg('sdc2')+
  facet_wrap(~sdc2)

plot_gene_ieg('cntn4')+
  facet_wrap(~cntn4)

plot_gene_ieg('bcan')+
  facet_wrap(~bcan)

plot_gene_ieg('gabra3')+
  facet_wrap(~gabra3)

plot_gene_ieg('grm5a')+
  facet_wrap(~grm5a)

plot_gene_ieg('pgr')+
  facet_wrap(~pgr)

plot_gene_ieg('cckbrb')+
  facet_wrap(~cckbrb)

plot_gene_ieg('tacr1a')+
  facet_wrap(~tacr1a)

plot_gene_ieg('esr2a')+
  facet_wrap(~esr2a)

plot_gene_ieg('esr2b')+
  facet_wrap(~esr2b)

plot_gene_ieg('ar')+
  facet_wrap(~ar)

## degs
degs =read.csv('/Users/ggraham/Desktop/multiome_poa/DEG Outputs/FINAL degs classified w singular.csv')
degs_6 = degs$gene[degs$cluster==6]

for(i in degs_6){
  p = plot_gene_ieg(paste0(i))+
  facet_wrap(~paste0(i))
  print(p)
}
