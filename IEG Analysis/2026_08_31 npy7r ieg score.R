
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
obj = subset(obj, Status %in% c('M', 'D','F' ) & final_clusters==6)
obj$npy7r = ifelse(obj@assays$RNA$data['npy7r',]>0, T, F)
Idents(obj) = 'npy7r'
# Find IEGs ---
iegs <- c('npas4a', 'fosb','egr1') #using fosb instead of cfos , I think it should be more or less the same

iegs_= calculateCoexpressedGenes(obj, iegs, clustering = 'npy7r') #c-fos like is included anyway so its good

RNA_data <- obj@assays$RNA$data%>%as.matrix()
ieg_rows_rna <- RNA_data[which(rownames(RNA_data) %in% iegs_),]
ieg_scores <- colSums(ieg_rows_rna>0)
obj$ieg_scores = ieg_scores


### beta binomial regression
library(glmmTMB)
library(performance)

out_dat = data.frame()
for(i in unique(obj$npy7r)){
  print(i)
  temp_data = subset(obj@meta.data, npy7r ==i & Status != 'NRM')
  
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


b=temp_data%>%
  group_by(individual, Status)%>%
  summarize(mean_ieg = mean(ieg_scores))

ggplot(b, aes(x= Status, y = mean_ieg))+
  geom_boxplot()+
  geom_point()

library(emmeans)
library(ggplot2)

# Refit (or reuse) the model for npy7r == TRUE
temp_data_T <- subset(obj@meta.data, npy7r == TRUE & Status != 'NRM')

model_T <- glmmTMB(
  cbind(ieg_scores, length(iegs_) - ieg_scores) ~ Status + (1 | individual),
  data = temp_data_T,
  family = betabinomial(link = "logit")
)

# Get estimated marginal means by Status, back-transformed to proportion scale
emm_T <- emmeans(model_T, ~ Status, type = "response")
emm_df <- as.data.frame(emm_T)

# emm_df should have columns like: Status, prob (or response), asymp.LCL, asymp.UCL
print(emm_df)

ggplot(emm_df, aes(x = Status, y = response, ymin = asymp.LCL, ymax = asymp.UCL)) +
  geom_point(size = 3) +
  geom_errorbar(width = 0.15) +
  labs(
    title = "Model-estimated IEG proportion by Status (npy7r = T)",
    y = "Estimated proportion of IEGs expressed",
    x = "Status"
  ) +
  theme_minimal(base_size = 14)



library(ggplot2)
library(dplyr)

# Extract egr1 expression from the RNA data matrix
egr1_expr <- RNA_data['egr1', ]
obj$egr1_expr <- egr1_expr
obj$egr1_pos <- ifelse(egr1_expr > 0, 1, 0)  # binary detected/not detected

# ---- Option A: continuous expression value, npy7r = T only ----
temp_data_T <- subset(obj@meta.data, npy7r == TRUE & Status != 'NRM')
temp_data_T$egr1_expr <- egr1_expr[rownames(obj@meta.data) %in% rownames(temp_data_T)]

b_egr1 <- temp_data_T %>%
  group_by(individual, Status) %>%
  summarize(mean_egr1 = mean(egr1_expr), .groups = "drop")

ggplot(b_egr1, aes(x = Status, y = mean_egr1)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position = position_jitter(width = 0.1)) +
  labs(
    title = "egr1 expression by Status (npy7r = T)",
    y = "Mean egr1 expression (per individual)",
    x = "Status"
  ) +
  theme_minimal(base_size = 14)
