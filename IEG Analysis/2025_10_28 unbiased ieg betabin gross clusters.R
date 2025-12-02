
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
obj = subset(obj, Status %in% c('M', 'D','F' ) & final_clusters %notin%c(1, 2, 11, 15, 26, 20, 13))
Idents(obj) = 'final_clusters'
# Find IEGs ---
iegs <- c('npas4a', 'fosb','egr1') #using fosb instead of cfos , I think it should be more or less the same

iegs_= calculateCoexpressedGenes(obj, iegs, clustering = 'final_clusters') #c-fos like is included anyway so its good

RNA_data <- obj@assays$RNA$data%>%as.matrix()
ieg_rows_rna <- RNA_data[which(rownames(RNA_data) %in% iegs_),]
ieg_scores <- colSums(ieg_rows_rna>0)
obj$ieg_scores = ieg_scores


### beta binomial regression
library(glmmTMB)
library(performance)

out_dat = data.frame()
for(i in unique(obj$final_clusters)){
  print(i)
  temp_data = subset(obj@meta.data, final_clusters ==i & Status != 'NRM')
  
  model <- glmmTMB(cbind(ieg_scores,length(iegs_)-ieg_scores)  ~ Status + (1|individual),
  data = temp_data,
  family=betabinomial(link = "logit"))  
  
  av = car::Anova(model, type = 'III') 

  newd <- data.frame(cluster = i,
                   status_p.value = av$`Pr(>Chisq)`[2],
                  singular =  check_singularity(model)
  )

  out_dat = rbind(out_dat, newd)
}

out_dat$q.value = p.adjust(out_dat$status_p.value, 'fdr', nrow(out_dat))
out_dat$q.value = round(out_dat$q.value, 5)








