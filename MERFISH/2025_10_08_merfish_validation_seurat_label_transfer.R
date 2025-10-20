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

add_transcript_length = function(gene){
  length = biomart_basic$transcript_length[biomart_basic$entrezgene_accession==gene]%>%min()
  if(length == Inf){length = 0}
  return(length)
}
add_max_expression = function(gene){
  print(gene)
  return(max(counts_matrix[gene,]))
  
}

## whip out biomart
ensembl <- useEnsembl(biomart = "genes", 
                      dataset = "aocellaris_gene_ensembl")

att = listAttributes(ensembl)

biomart_basic <-getBM(
    mart = ensembl, #working mart 
    attributes = c("entrezgene_accession",
                   'transcript_length'))
}
obj = readRDS("~/Desktop/optimal_clustering_rna_only.rds")
counts_matrix = obj@assays$RNA$counts%>%as.matrix()

### start -----

all_marks_6_vs_0 <- FindMarkers(obj, 
                                ident.1 = '6', 
                                ident.2 = c('0',3),
                                group.by = 'final_clusters')

#candidate genes
candidate_gene = read.xlsx("MERFISH/candidate_genes.xlsx")
candidate_genes = data.frame(gene =candidate_gene$gene%>%unique())
candidate_genes$gene[candidate_genes$gene =='fgf12a'] = 'fosb'
candidate_genes$gene[candidate_genes$gene =='grb10b'] = 'npas4a'

#candidate_genes[100:122, 'gene'] <- rownames(all_marks_6_vs_0[1:23,])


candidate_genes <- candidate_genes%>%
  distinct(gene)

candidate_genes$gene[candidate_genes$gene =='gal'] =rownames(all_marks_6_vs_0[24,])
candidate_genes$gene[candidate_genes$gene =='prlh2'] =rownames(all_marks_6_vs_0[25,])
candidate_genes$gene[candidate_genes$gene =='chgb'] =rownames(all_marks_6_vs_0[26,])
candidate_genes$gene[candidate_genes$gene =='LOC111564112'] =rownames(all_marks_6_vs_0[27,])

candidate_genes$max_expression = sapply(candidate_genes$gene, add_max_expression)
candidate_genes$transcript_length = sapply(candidate_genes$gene, add_transcript_length)

# manually filling based on NIH
candidate_genes$transcript_length[candidate_genes$gene =='LOC111568258'] = 3021
candidate_genes$transcript_length[candidate_genes$gene =='LOC111564112'] =9486
candidate_genes$transcript_length[candidate_genes$gene =='LOC129349764'] =2437
candidate_genes$transcript_length[candidate_genes$gene =='chgb'] =2365
candidate_genes$transcript_length[candidate_genes$gene =='ttn.2'] =87356

merfish_genes =candidate_genes

# test find transfer labels -- thank you claude
library(Seurat)
library(dplyr)
library(ggplot2)

# MERFISH Gene Panel Validation
# Tests if selected genes can accurately map held-out cells back to original clusters

# Parameters
holdout_fraction <- 0.2  # Fraction of cells to hold out for testing
n_neighbors <- 30        # Number of neighbors for label transfer
seed <- 42

set.seed(seed)

# Ensure merfish_genes exist in the object
merfish_genes_filtered <- merfish_genes$gene[merfish_genes$gene %in% rownames(obj)]
cat(sprintf("Using %d/%d MERFISH genes found in object\n", 
            length(merfish_genes_filtered), length(merfish_genes$gene)))

# Split cells into reference and query sets
all_cells <- colnames(obj)
n_test <- round(length(all_cells) * holdout_fraction)
test_cells <- sample(all_cells, n_test)
ref_cells <- setdiff(all_cells, test_cells)

cat(sprintf("\nReference cells: %d\nTest cells: %d\n", 
            length(ref_cells), length(test_cells)))

# Create reference and query objects
ref_obj <- subset(obj, cells = ref_cells)
query_obj <- subset(obj, cells = test_cells)

# Store original cluster labels from query
true_labels <- query_obj$final_clusters
names(true_labels) <- colnames(query_obj)

# Subset to MERFISH genes only
ref_subset <- subset(ref_obj
                     #, features = merfish_genes_filtered
                     )
query_subset <- subset(query_obj
                      # , features = merfish_genes_filtered
                       )

# Normalize and scale data
ref_subset <- NormalizeData(ref_subset, verbose = FALSE)
ref_subset <- ScaleData(ref_subset, verbose = FALSE)
query_subset <- NormalizeData(query_subset, verbose = FALSE)

# Find anchors and transfer labels
cat("\nFinding transfer anchors...\n")
anchors <- FindTransferAnchors(
  reference = ref_subset,
  query = query_subset,
  dims = 1:min(30, length(merfish_genes_filtered) - 1),
  features = merfish_genes$gene,
  verbose = FALSE
)

cat("Transferring labels...\n")
predictions <- TransferData(
  anchorset = anchors,
  refdata = ref_subset$final_clusters,
  dims = 1:min(30, length(merfish_genes_filtered) - 1),
  k.weight = n_neighbors,
  verbose = FALSE
)

# Add predictions to query object
query_obj$predicted_clusters <- predictions$predicted.id
query_obj$prediction_score <- predictions$prediction.score.max

# Calculate accuracy metrics
accuracy <- mean(query_obj$predicted_clusters == true_labels)
cat(sprintf("\n=== RESULTS ===\n"))
cat(sprintf("Overall Accuracy: %.2f%%\n", accuracy * 100))

# Per-cluster accuracy
cluster_accuracy <- query_obj@meta.data %>%
  group_by(final_clusters) %>%
  summarise(
    n_cells = n(),
    correct = sum(predicted_clusters == final_clusters),
    accuracy = mean(predicted_clusters == final_clusters),
    mean_score = mean(prediction_score)
  ) %>%
  arrange(desc(accuracy))

print(cluster_accuracy)

# Confusion matrix
conf_matrix <- table(True = true_labels, Predicted = query_obj$predicted_clusters)
cat("\nConfusion Matrix:\n")
print(conf_matrix)

# Calculate mean prediction score
mean_score <- mean(query_obj$prediction_score)
cat(sprintf("\nMean Prediction Score: %.3f\n", mean_score))

# Visualize results
p1 <- ggplot(query_obj@meta.data, aes(x = final_clusters, fill = predicted_clusters)) +
  geom_bar(position = "fill") +
  theme_minimal() +
  labs(title = "Prediction Distribution by True Cluster",
       x = "True Cluster", y = "Proportion", fill = "Predicted") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p2 <- ggplot(query_obj@meta.data, aes(x = prediction_score)) +
  geom_histogram(bins = 50, fill = "steelblue", alpha = 0.7) +
  theme_minimal() +
  labs(title = "Distribution of Prediction Scores",
       x = "Prediction Score", y = "Count") +
  geom_vline(xintercept = mean_score, linetype = "dashed", color = "red")

print(p1)
#print(p2)

# Summary assessment
cat("\n=== ASSESSMENT ===\n")
if (accuracy > 0.9 & mean_score > 0.7) {
  cat("✓ EXCELLENT: Gene panel is highly suitable for MERFISH mapping\n")
} else if (accuracy > 0.8 & mean_score > 0.6) {
  cat("✓ GOOD: Gene panel should work well for MERFISH mapping\n")
} else if (accuracy > 0.7 & mean_score > 0.5) {
  cat("⚠ MODERATE: Gene panel may need refinement\n")
} else {
  cat("✗ POOR: Consider adding more informative genes\n")
}

# Identify problematic clusters
low_acc_clusters <- cluster_accuracy %>%
  filter(accuracy < 0.7) %>%
  pull(final_clusters)

if (length(low_acc_clusters) > 0) {
  cat(sprintf("\nClusters with <70%% accuracy: %s\n", 
              paste(low_acc_clusters, collapse = ", ")))
}




### what the fuck -----

obj_sub_6 = FindSubCluster(obj, '6', graph = 'harmony.wsnn')

for(i in paste0('6_', 0:3)){
  marks = FindMarkers(obj_sub_6, ident.1 = i, ident.2 = 0, group.by='sub.cluster')
  assign(paste0('best_10_',i),head(marks,20), envir   = .GlobalEnv )
  print(head(marks))
}

## look at query object
names_6_3 = obj@meta.data[obj_sub_6$sub.cluster=='6_3',]%>%rownames()
names_6_2 = obj@meta.data[obj_sub_6$sub.cluster=='6_2',]%>%rownames()
names_6_1 = obj@meta.data[obj_sub_6$sub.cluster=='6_1',]%>%rownames()
names_6_0 = obj@meta.data[obj_sub_6$sub.cluster=='6_0',]%>%rownames()

query_obj@meta.data$predicted_clusters[rownames(query_obj@meta.data) %in% c(names_6_3)]%>%table()%>%sort()
query_obj@meta.data$predicted_clusters[rownames(query_obj@meta.data) %in% c(names_6_2)]%>%table()%>%sort()
query_obj@meta.data$predicted_clusters[rownames(query_obj@meta.data) %in% c(names_6_1)]%>%table()%>%sort()
query_obj@meta.data$predicted_clusters[rownames(query_obj@meta.data) %in% c(names_6_0)]%>%table()%>%sort()

# its literally all of them EXCEPT 6_3

# How many strong markers exist?
strong_markers <- all_marks_6_vs_0[abs(all_marks_6_vs_0$avg_log2FC) > 1 & 
                                   all_marks_6_vs_0$p_val_adj < 0.01,]

misclassified_6_cells = rownames(query_obj@meta.data[query_obj$final_clusters == 6 &query_obj$predicted_clusters!=6,])
cells_0  = rownames(query_obj@meta.data[query_obj$final_clusters == 0,])
# Assign a temporary grouping
Idents(query_obj) <- "tmp_ident"
query_obj$tmp_ident <- "other"
query_obj$tmp_ident[misclassified_6_cells] <- "misclassified6"
query_obj$tmp_ident[cells_0] <- "cells_0"

# Run FindMarkers using those groups
misclassified_markers = FindMarkers(
  query_obj,
  ident.1 = "misclassified6",
  ident.2 = "cells_0",
  group.by = 'tmp_ident'
)

head(misclassified_markers)
