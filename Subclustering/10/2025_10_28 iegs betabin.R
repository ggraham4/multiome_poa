# cluster 10 beta binomial IEG analysis


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

  plot_gene_ieg = function(gene_name, seurat_obj = sub_10, statuses = c('M','D','F')) {
  
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
    group_by(individual, Status, !!sym(gene_name), sub.cluster) %>%
    summarize(mean_ieg = mean(ieg_PC1), .groups = 'drop') %>%
    subset(Status %in% statuses)
  
  # Create plot
  p = ggplot(plot_data, aes(x = Status, y = mean_ieg, color = !!sym(gene_name))) +
    geom_boxplot() +
    geom_point(position = position_jitterdodge(.5)) +
    facet_wrap(~sub.cluster, scales = 'free') +
    labs(title = paste("IEG expression by", gene_name, "status"),
         y = "IEG PC1",
         color = gene_name) +
    theme_minimal()
  
  return(p)
}

  #define functions

prop_cluster_plot=function(object, gene, cluster, clustering = 'final_clusters'){
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

deg_plotter = function(object = rgc, 
                       gene, 
                       cluster , 
                       clustering='sub.cluster',
                       signif_dm  ,
                       signif_df ,
                       signif_mf ,
                       singular=F,
                       common_name){
  set.seed(10)
  singular = ifelse(singular == T, 'Singular', '')
  meta= object@meta.data
  meta$gene = object@assays$RNA$data[gene,]
  
  meta_grouped_and_sded = meta%>%
    filter(Phase != 'NRM' & !!sym(clustering) == cluster) %>%
    group_by(individual, Phase)%>%
    summarize(mean_gene = mean(gene),
              se_gene = sd(gene)/sqrt(n()))
  
  plot_lower_lim = min(meta_grouped_and_sded$mean_gene -meta_grouped_and_sded$se_gene )
  plot_upper_lim= max(meta_grouped_and_sded$mean_gene +meta_grouped_and_sded$se_gene ) * 1.2
  plot_signif_lower = max(meta_grouped_and_sded$mean_gene +meta_grouped_and_sded$se_gene ) * 1.05
  plot_signif_upper = max(meta_grouped_and_sded$mean_gene +meta_grouped_and_sded$se_gene ) * 1.25
  

  
    signif_dm = p_annotate(signif_dm)
    signif_df = p_annotate(signif_df)
    signif_mf = p_annotate(signif_mf)
  
  textsize_dm = ifelse(grepl("\\*", signif_dm), 10, 3)  
  textsize_df = ifelse(grepl("\\*", signif_df), 10, 3)  
  textsize_mf = ifelse(grepl("\\*", signif_mf), 10, 3)    

  
      plot_upper_lim= ifelse(signif_mf!= 'NS', max(meta_grouped_and_sded$mean_gene +meta_grouped_and_sded$se_gene ) * 1.4,plot_upper_lim )

  
plot = ggplot(meta_grouped_and_sded, aes(x = Phase, y = mean_gene,fill = Phase))+
    geom_boxplot(alpha = 0.5,  outlier.shape = NA)+
  geom_pointrange(aes(x = Phase,
                      y = mean_gene,
                      ymin = mean_gene-se_gene,
                      ymax= mean_gene+se_gene),
                  position = position_jitterdodge(1), 
                  size = 0.2
                  )+
  labs(y = 'Mean +/- SE Expression', subtitle = singular)+
  ggtitle(paste0(common_name, ': ', cluster))+
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5, size =12),
        plot.subtitle = element_text(hjust = 0.5, size =8))+
  theme(legend.position = 'none')+
    geom_signif(xmin = c(1.0),
              xmax = c(1.9),
              y_position = c(plot_signif_lower),
              annotation =c(signif_dm), 
              color = "black",
              tip_length = c(0,0),
              textsize=textsize_dm)+
  geom_signif(xmin = c(2.1),
              xmax = c(5),
              y_position = c(plot_signif_lower),
              annotation =c(signif_df), 
              color = "black",
              tip_length = c(0,0),
              textsize=textsize_df)+
  ylim(plot_lower_lim, plot_upper_lim)
   
 if(signif_mf != 'NS'){
   plot  <- plot+
      geom_signif(xmin = c(1),
              xmax = c(5),
              y_position = c(plot_signif_upper),
              annotation =c(signif_mf), 
              color = "black",
              tip_length = c(0,0),
              textsize=textsize_mf)
    
  }

return(plot)
  
}

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

sub_10 = FindSubCluster(obj, 10, graph.name='harmony.wsnn')
Idents(sub_10) <- 'sub.cluster'
sub_10 = subset(sub_10, final_clusters ==10)
sub_10$sub.cluster = factor(sub_10$sub.cluster, levels = c(paste0('10_',0:1)))
sub_10$Status = factor(sub_10$Status, levels = c('NRM','M',"D",'E','NF','F'))

DimPlot(sub_10)

# Find IEGs ---
iegs <- c('npas4a', 'fosb','egr1') #using fosb instead of cfos , I think it should be more or less the same

iegs_= calculateCoexpressedGenes(sub_10, iegs) #c-fos like is included anyway so its good

RNA_data <- sub_10@assays$RNA$data%>%as.matrix()
ieg_rows_rna <- RNA_data[which(rownames(RNA_data) %in% iegs_),]
ieg_scores <- colSums(ieg_rows_rna>0)

sub_10$ieg_scores = ieg_scores
### beta binomial regression
library(glmmTMB)
library(performance)

data_10_0 = subset(sub_10@meta.data, sub.cluster =='10_0')

model <- glmmTMB(cbind(ieg_scores,length(iegs_)-ieg_scores)  ~ Status + (1|individual),
  data = data_10_0,
  family=betabinomial(link = "logit"))  

car::Anova(model, type = 'III') # the last gasp fails

emm <- emmeans(model, ~ Status, type = "response")

# Convert emmeans to a data frame
emm_df <- summary(emm)

# Plot predicted proportions with confidence intervals
ggplot(emm_df, aes(x = Status, y = response)) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), width=0.2) +
  ylab("Predicted Proportion of IEGs Expressed") +
  theme_minimal() # ah that is a problem

# what about other sub.clusters
out_dat = data.frame()
for(i in unique(sub_10$sub.cluster)){
  
  temp_data = subset(sub_10@meta.data, sub.cluster ==i & Status != 'NRM')
  
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

out_dat$q.value = p.adjust(out_dat$status_p.value, 'fdr', 4)
out_dat$q.value = round(out_dat$q.value, 5)


i = '10_0'

