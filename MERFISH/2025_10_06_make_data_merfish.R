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
p_annotate <- function(p_value) {
  if (is.na(p_value)) {
    return("NA")
  }
  
  if (p_value < 0.001) {
    return("***")
  } else if (p_value < 0.01) {
    return("**")
  } else if (p_value < 0.05) {
    return("*")
  } else if (p_value < 0.1) {
    return(paste0("p = ", round(p_value, 3)))
  } else {
    return("NS")
  }
}

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
  
  textsize_dm = ifelse(grepl("\\*", signif_dm), 6, 3)  
  textsize_df = ifelse(grepl("\\*", signif_df), 6, 3)  
  textsize_mf = ifelse(grepl("\\*", signif_mf), 6, 3)    

  
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


coalesce_interaction = function(LRpair, CCpair, whole_sub){
  "
  the goal here is to coalesce the means of each individual for the given
  cluster cluster pair and ligand receptor pair
  LRpair : Ligand - receptor pair, string e.g., 'TAC1_TACR3'
  CCpair : Cluster - cluster pair, sring e.g., '0|10'
  whole_sub: which list to reference, should be of level 'whole' or 'sub'
  "
  
  if(whole_sub == 'whole'){
    ref_signif = whole_significant
    ref_means = whole_individual_data
  }
  if(whole_sub == 'sub'){
    ref_signif = sub_significant
    ref_means = sub_individual_data
    
  }
  
  coalesced_data = data.frame()
  for(ind in individuals){
    individual_df =  ref_means[[ind]]
    CCpair_column = which(colnames(individual_df)==CCpair)
    LRpair_row = which(individual_df$interacting_pair== LRpair)
    
    individual_mean = (individual_df[LRpair_row, CCpair_column]) %>% as.numeric()
    if(length(individual_mean) == 0 || is.na(individual_mean)){
      next  # Skip this individual if mean is NA or empty
    }    
    newd = data.frame(individual = ind,
                      mean = individual_mean)
    coalesced_data =rbind(coalesced_data, newd)
  }
  full_data = coalesced_data%>%
    right_join(measures, by = join_by('individual'=='Fish'))%>%
    subset(Status %in% c('M','D','E','NF','F'))
  full_data$Status = factor(full_data$Status, levels = c('M','D','E','NF','F'))
  return(full_data)
}

plot_interaction = function(LRpair, CCpair, whole_sub){
  if(whole_sub=='whole'){
    status_q_value = whole_significant$main_effect_q.value[whole_significant$cluster_pair == CCpair & 
                                                             whole_significant$interacting_pair == LRpair]
    status_q_value = round(status_q_value, 5)
    
  } else if(whole_sub=='sub'){
    status_q_value = sub_significant$main_effect_q.value[sub_significant$cluster_pair == CCpair & 
                                                           sub_significant$interacting_pair == LRpair]
    status_q_value = round(status_q_value, 5)
  } else {
    message('Plot ERROR: whole_sub not "whole" or "sub"')
    status_q_value = 'Invalid whole_sub parameter'
  }
  
  if(length(status_q_value) == 0 || is.na(status_q_value)){
    status_q_value = ' Interaction not relevant'
  }  
  
  data_for_plot = coalesce_interaction(LRpair, CCpair, whole_sub)
  
  data_for_plot = data_for_plot[!is.na(data_for_plot$mean),]
  plot = ggplot(data_for_plot, aes(x = Status, y = mean, 
                                   color = Status, 
                                   fill = Status, 
                                   shape = Status))+
    geom_point(position = position_jitterdodge(2))+
    geom_boxplot(alpha = 0.25, outlier.shape = NA)+
    theme_minimal()+
    labs(x = 'Status', y = 'Interaction Mean',
         title = paste0(LRpair, ' ', CCpair),
         subtitle = paste0('Status q.value =',status_q_value ))
  
  return(plot)
  
}


}
obj = readRDS("~/Desktop/optimal_clustering_rna_only.rds")

# gonna have to do markers for 6 and 1 subclusters separately

all_markers = FindAllMarkers(obj)

#subcluster
sub_6 = FindSubCluster(obj, 6, graph = 'harmony.wsnn')
sub_6_markers = FindAllMarkers(sub_6, group.by = 'sub.cluster')

sub_1 = FindSubCluster(obj, 1, graph = 'harmony.wsnn')
sub_1_markers = FindAllMarkers(sub_1, group.by = 'sub.cluster')


#join_dfs
sub_6_sub_markers = subset(sub_6_markers,
                           cluster %in% c(paste0(6,'_', 0:3)) &
                             p_val_adj < 0.05 & 
                             pct.1 > pct.2)
sub_1_sub_markers = subset(sub_1_markers, cluster %in% c(paste0('1_', 0:5))&
                             p_val_adj < 0.05 & 
                             pct.1 > pct.2)

whole_markers = subset(all_markers, cluster %notin% c(1, 6) &
                             p_val_adj < 0.05 & 
                             pct.1 > pct.2 )

joined_markers = rbind(whole_markers,sub_1_sub_markers,  sub_6_sub_markers)


## whip out biomart
ensembl <- useEnsembl(biomart = "genes", 
                      dataset = "aocellaris_gene_ensembl")

att = listAttributes(ensembl)

biomart_basic <-getBM(
    mart = ensembl, #working mart 
    attributes = c("entrezgene_accession",
                   'transcript_length'))

# add transcript length
add_transcript_length = function(gene){
  length = biomart_basic$transcript_length[biomart_basic$entrezgene_accession==gene]%>%min()
  if(length == Inf){length = 0}
  return(length)
}

joined_markers$transcript_length = sapply(joined_markers$gene, add_transcript_length)

#max expression 
counts_matrix = obj@assays$RNA$counts%>%as.matrix()
add_max_expression = function(gene){
  return(max(counts_matrix[gene,]))
  
}
joined_markers$max_expression = lapply(joined_markers$gene, add_max_expression)

## filter time
joined_markers_filtered = subset(joined_markers, max_expression <= 136 & transcript_length>= 1440)

unique_good_genes = unique(joined_markers_filtered$gene)

# output
truncated_matrix = counts_matrix[unique_good_genes, ]

sparse_mat = Matrix(truncated_matrix, sparse = T)
library(Matrix)

writeMM(sparse_mat, file = "/Users/ggraham/MERFISH/2025_10_06_genes_meeting_criteria.mtx")
test_pull =readMM("/Users/ggraham/MERFISH/2025_10_06_genes_meeting_criteria.mtx")

write.csv(colnames(sparse_mat), "/Users/ggraham/MERFISH/2025_10_06_sparse_mat_colnames.csv")
write.csv(rownames(sparse_mat), "/Users/ggraham/MERFISH/2025_10_06_sparse_mat_rownames.csv")

dict = data.frame(  clusterName = joined_markers_filtered$cluster,
                  markers = joined_markers_filtered$gene
                )

write.csv(dict, "/Users/ggraham/MERFISH/dict.csv")

### clusters list 
obj$final_clusters <- unfactor(obj$final_clusters )

obj$final_clusters2 <- ifelse(obj$final_clusters==1, sub_1$sub.cluster, obj$final_clusters)
obj$final_clusters2 <- ifelse(obj$final_clusters==6, sub_6$sub.cluster, obj$final_clusters2)

DimPlot(obj, group.by = 'final_clusters2')
DimPlot(subset(obj, final_clusters==6), group.by = 'final_clusters2')

dat = data.frame(clusters = obj$final_clusters2, 
                 cell_label = rownames(obj@meta.data))
#write.csv(dat, "/Users/ggraham/MERFISH/2025_10_06_clusters.csv")

#### top 3 marker genes from each cluster
joined_markers_filtered$transcript_length <- as.numeric(joined_markers_filtered$transcript_length)
joined_markers_filtered$max_expression <- as.numeric(joined_markers_filtered$max_expression)

joined_markers_log2fc <- joined_markers_filtered %>%
  group_by(cluster) %>%
  distinct(gene, .keep_all = TRUE) %>%
  slice_max(order_by = avg_log2FC, n = 3)

joined_markers_pval <- joined_markers_filtered %>%
  group_by(cluster) %>%
  distinct(gene, .keep_all = TRUE) %>%
  slice_max(order_by = p_val_adj, n = 3)

write.csv(joined_markers_pval, "/Users/ggraham/MERFISH/2025_10_06_joined_markers_pval.csv")
write.csv(joined_markers_log2fc, "/Users/ggraham/MERFISH/2025_10_06_joined_markers_log2fc.csv")

###
NSmarkers =read_csv("/Users/ggraham/MERFISH/NSForest Output/final_clusters2_markers.csv")
  geneNamer = function(gene){
  names = read.csv('Reference/genes updated.csv')
  
  name = names$NIH_description[names$NIH_accession==gene][1]
  
  if(is.na(name)){name = gene}
  return(name)
  }
  
  NSmarkers$longName = sapply(NSmarkers$markerGene, geneNamer)
write.csv(NSmarkers, "MERFISH/2025_10_06_NSMarkers_named.csv")

