#Negative binomial lower stringency 
{
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
  library(Polychrome)
  P40 <- createPalette(40, c("#FF0000", "#00FF00", "#0000FF"), range = c(30, 80))
swatch(P40)
names(P40) <- NULL

  mean_expression_cluster_plot<- readRDS('Functions/mean_expression_cluster_plot.rds')
prop_cluster_plot<- readRDS( 'Functions/prop_cluster_plot.rds')
mean_expression_cluster_data<- readRDS('Functions/mean_expression_cluster_data.rds')
clown_go<- readRDS('Functions/clown_go')
define_degs<- readRDS('Functions/define_degs')

}

obj <- readRDS('/Users/ggraham/Desktop/snRNA-seq R Files 122524/RNA Object.rds')


neg.bin.mult <- function(obj,
                         cluster, 
                         clustering = 'harmony.wnn_res0.4_clusters',
                         n_cores = detectCores() - 1) {
  start_time <- Sys.time()  # Start timing
  
  message('Extracting Counts')
  counts <- obj@assays$RNA$counts[, obj@meta.data[[clustering]] == cluster & (obj@meta.data$Status == "M" | obj@meta.data$Status == "F" | obj@meta.data$Status == "D")]
  combined_counts <- counts
  
  df_counts <- data.frame(t(combined_counts))
  colnames(df_counts) <- rownames(obj@assays$RNA)
  
  n_genes = ncol(df_counts)
  n_cells = nrow(df_counts)
  
  message("Making Counts Data Frame...")
  df_counts_meta <- data.frame(rownames(df_counts))
  df_counts_meta$id <- df_counts_meta$rownames.df_counts.
  df_counts_meta$rownames.df_counts. = NULL
  df_counts_meta$individual = obj$individual[obj@meta.data[[clustering]] == cluster & (obj@meta.data$Status == "M" | obj@meta.data$Status == "F" | obj@meta.data$Status == "D")]
  df_counts_meta$Status = obj$Status[obj@meta.data[[clustering]] == cluster & (obj@meta.data$Status == "M" | obj@meta.data$Status == "F" | obj@meta.data$Status == "D")]
  
  message("Removing Genes with 0 Counts...")
  df_counts_no_0 <- df_counts[, which(colSums(df_counts) != 0)]
  
  message("Making New Counts Data Frame Without 0s...")
  n_genes_no_0 = ncol(df_counts_no_0)
  
  
  ##### ok here is where my changes are going to have to be #####
  
  df_counts_no_0 <- cbind(df_counts_no_0, df_counts_meta)
  df_counts_no_0_split_by_subject <- split(df_counts_no_0, f = df_counts_no_0$individual)
  
  message("Finding Good Genes for Subject...")
  # REMOVE GENES WITH ZERO COUNTS IN EACH SUBJECT 
  for (l in 1:length(df_counts_no_0_split_by_subject)) {
    correct_gene_names <- colnames(df_counts_no_0)
    
    temp_subject_l <- data.frame(df_counts_no_0_split_by_subject[[l]]) ### AND HERE THEY GET FUCKED UP 
    colnames(temp_subject_l) <- correct_gene_names
    
    temp_subject_l_counts <- temp_subject_l[, 1:n_genes_no_0]
    ###
    temp_subject_l_counts_no_0 <- temp_subject_l_counts #<- temp_subject_l_counts[, which(colSums(temp_subject_l_counts) != 0)]
    #ok here I am making the real code a comment to stop the filtering without fucking up the rest of the code
    #I'm realizing there does need to be some way I test for a gene being missing in several sexes huh
    out <- data.frame(colnames(temp_subject_l_counts_no_0))
    assign(x = paste0("gene_list_subject_", l), value = get("out"))
  }
  
  # GENERATE A LIST OF GENES FOR EACH SUBJECT -- I believe everything else should still work
  good_gene_list <- gene_list_subject_1$colnames.temp_subject_l_counts_no_0.
  
  for (m in 2:length(df_counts_no_0_split_by_subject)) {
    temp_good_gene_list_m <- data.frame(value = get(paste0("gene_list_subject_", m)))
    temp_good_gene_list_m <- temp_good_gene_list_m$colnames.temp_subject_l_counts_no_0.
    good_gene_list <- intersect(good_gene_list, temp_good_gene_list_m)
  }
  
  p <- length(good_gene_list)
  
  message('Making Gene Data Frame for Each Subject...')
  valid_genes <- good_gene_list[good_gene_list %in% colnames(df_counts_no_0)]
  v <- length(valid_genes)
  
  df_counts_no_0_all_subjects <- df_counts_no_0[, valid_genes]
  count_matrix_final <- as.matrix(df_counts_no_0_all_subjects)
  count_matrix_final <- as.data.frame(t(count_matrix_final))
  df_counts_no_0_all_subjects <- cbind(df_counts_no_0_all_subjects, df_counts_meta)
  
  Status <- as.factor(df_counts_no_0_all_subjects$Status)
  subject <- as.factor(df_counts_no_0_all_subjects$individual)
  
  message("Estimating Dispersion Using Gamma-Poisson...")
  cluster_size <- ncol(count_matrix_final)
  
  size_factors <- calculateSumFactors(count_matrix_final,
                                      clusters = NULL,
                                      ref.clust = NULL,
                                      max.cluster.size = cluster_size,
                                      positive = TRUE,
                                      scaling = NULL,
                                      min.mean = NULL,
                                      subset.row = NULL)
  
  coldata <- data.frame(Status)
  fit <- glm_gp(as.matrix(count_matrix_final), col_data = coldata, size_factors = size_factors, design = ~ Status, on_disk = FALSE)
  dispersions.RAW <- fit$overdispersion_shrinkage_list$ql_disp_estimate
  log.sizeFactors.RAW <- log(size_factors)
  
  offset <- log.sizeFactors.RAW
  index <- v
  
  results <- mclapply(1:index, function(i) {
    message('Calculating Gene', paste0(i, ' of ', index, '...'))
    
    dispersion <- dispersions.RAW[i]
    outcome <- df_counts_no_0_all_subjects[, i]
    
    
     tryCatch({suppressMessages(
    glmer_model <- glmer(outcome ~ Status + (1 | subject),
                         offset = offset,
                         family = MASS::negative.binomial(theta = 1 / dispersion)))
     }, error = function(e) {
      return(NULL)
    })
    
    pairs <- pairs(emmeans(glmer_model, 'Status'), adjust = 'none')
    
    output_df <- data.frame(
      gene = valid_genes[i],
      f_m_estimate = as.data.frame(pairs[pairs@grid$contrast == 'F - M'])[, 2], # Ensure the correct reference to pairs columns
      f_m_p.value = as.data.frame(pairs[pairs@grid$contrast == 'F - M'])[, 6],
      d_m_estimate = as.data.frame(pairs[pairs@grid$contrast == 'D - M'])[, 2],
      d_m_p.value = as.data.frame(pairs[pairs@grid$contrast == 'D - M'])[, 6],
      d_f_estimate = as.data.frame(pairs[pairs@grid$contrast == 'D - F'])[, 2],
      d_f_p.value = as.data.frame(pairs[pairs@grid$contrast == 'D - F'])[, 6],
      warning = ifelse(length(glmer_model@optinfo$conv$lme4$code) != 0, substr(glmer_model@optinfo$conv$lme4$messages, 1, 50), NA),
      singular = ifelse(isSingular(glmer_model), TRUE, FALSE)
    )
    return(output_df)
  }, mc.cores = n_cores
  )
  
  results <- do.call(rbind, results)
  results <- as.data.frame(results, stringsAsFactors = FALSE)
  results$f_m_p.value <- ifelse(results$f_m_p.value == 0, 1, results$f_m_p.value)
  results$d_m_p.value <- ifelse(results$d_m_p.value == 0, 1, results$d_m_p.value)
  results$d_f_p.value <- ifelse(results$d_f_p.value == 0, 1, results$d_f_p.value)
  
  results$f_m_q.value <- ifelse(is.na(results$warning), p.adjust(as.numeric(results$f_m_p.value), method = "fdr"), "NA")
  results$d_m_q.value <- ifelse(is.na(results$warning), p.adjust(as.numeric(results$d_m_p.value), method = "fdr"), "NA")
  results$d_f_q.value <- ifelse(is.na(results$warning), p.adjust(as.numeric(results$d_f_p.value), method = "fdr"), "NA")
  
  assign(paste0("results_", "cluster", cluster), results, envir = .GlobalEnv)
  
  message('Complete')
  end_time <- Sys.time()  # End timing
  message(end_time - start_time)  # Print the time difference
  return(results)
}

for (i in 31:0) {
  print(i)
  output <- neg.bin.mult(obj,
                         i,
                         'harmony.wnn_res0.4_clusters')
  write.csv(output, paste0('/Volumes/jrhodes/Fish Lab/Experiments/sex change single nuc POA/Seurat Outputs/012425 Neg Bin w Doms Lower Stringency/cluster_', i, '.csv'))
}

### ANALYSIS ####
together_data <- data.frame()
for(i in 0:31){
  data <- read.csv(paste0('/Volumes/jrhodes/Fish Lab/Experiments/sex change single nuc POA/Seurat Outputs/012425 Neg Bin w Doms Lower Stringency/cluster_', i, '.csv'))
  #data <- get(paste0('results_cluster',i))
  data$cluster <- i
  together_data <- rbind(together_data, data)
}

together_data_defined <- define_degs(together_data)

together_data_defined_summed <- together_data_defined%>%
  subset(!is.na(class))%>%
  group_by(class, cluster)%>%
  summarize(class_count = n())

ggplot(together_data_defined_summed, aes(x = cluster, y = class_count, fill = class)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(x = "Cluster", y = "Number of DEGs") +
  geom_bar(position="stack", stat="identity")+
  theme(axis.text.x = element_text(angle = -45, vjust = 1, hjust=0))+
    scale_x_continuous(breaks = c(0:31))+
  scale_y_continuous()+
  scale_fill_manual(values = P40)

#### Chisq test for enrichment ####
together_data_chisq <-  together_data_defined%>%
  subset(!is.na(class))%>%
  group_by(cluster)%>%
  summarize(class_count = n())
total_sum <- sum(together_data_chisq$class_count)
mean_sum <- total_sum / nrow(together_data_chisq)
chisq_test <- chisq.test(together_data_chisq$class_count, p = rep(1/nrow(together_data_chisq), nrow(together_data_chisq)))
expected <- chisq_test$expected
residuals <- (together_data_chisq$class_count - expected) / sqrt(expected)
together_data_chisq$residuals <- residuals
together_data_chisq$significant <- together_data_chisq$residuals > 2
together_data_chisq$issignif <- ifelse(together_data_chisq$significant==T, '*',NA)

together_data_defined_summed_plot <- together_data_defined_summed%>%
  right_join(together_data_chisq, by = 'cluster')

ggplot(together_data_defined_summed_plot, aes(x = cluster, y = class_count.x, fill = class)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(x = "Cluster", y = "Number of DEGs") +
  geom_bar(position="stack", stat="identity")+
  theme(axis.text.x = element_text(angle = -45, vjust = 1, hjust=0))+
    scale_x_continuous(breaks = c(0:31))+
  scale_y_continuous()+
  scale_fill_manual(values = P40)+
  geom_text(aes(label = issignif, y = class_count.y), size = 10)

degs_data <- together_data_defined[!is.na(together_data_defined$class),c(2,14,15)]

ensembl <- useEnsembl(biomart = "genes", 
                      dataset = "aocellaris_gene_ensembl")
biomart_basic <-
  getBM(
    mart = ensembl, #working mart 
    attributes = c("entrezgene_accession",
                   'entrezgene_description'))

clust_2_go <- clown_go(degs_data$gene[degs_data$cluster==2])
dotplot(clust_2_go)

clust_3_go <- clown_go(degs_data$gene[degs_data$cluster==3])
dotplot(clust_3_go)

clust_6_go <- clown_go(degs_data$gene[degs_data$cluster==6])
dotplot(clust_6_go)

clust_7_go <- clown_go(degs_data$gene[degs_data$cluster==7])
dotplot(clust_7_go)
clust_7_go$geneID
clust_7_degs <- degs_data$gene[degs_data$cluster==7]
named_7 <- biomart_basic[biomart_basic$entrezgene_accession %in%clust_7_degs,]


clust_8_go <- clown_go(degs_data$gene[degs_data$cluster==8])
dotplot(clust_8_go)
clust_8_go$geneID
clust_8_degs <- degs_data$gene[degs_data$cluster==8]


named_8 <- biomart_basic[biomart_basic$entrezgene_accession %in%clust_8_degs,]


clust_14_go <- clown_go(degs_data$gene[degs_data$cluster==14])
dotplot(clust_14_go)
clust_14_degs <- degs_data$gene[degs_data$cluster==14]
clust_14_go$geneID

biomart_basic[biomart_basic$entrezgene_accession %in%clust_14_degs,]


clust_19_go <- clown_go(degs_data$gene[degs_data$cluster==19])
dotplot(clust_19_go)
clust_19_degs <- degs_data$gene[degs_data$cluster==19]
named_19 <- biomart_basic[biomart_basic$entrezgene_accession %in%clust_19_degs,]


clust_29_go <- clown_go(degs_data$gene[degs_data$cluster==29])
dotplot(clust_29_go)


### are any DEGs repeated ####
length(unique(degs_data$gene))
length(degs_data$gene)###Oh several are repeated

degs_data$gene%>%table()%>%sort()

##### wdfy2 ####
degs_data$cluster[degs_data$gene == 'wdfy2']
degs_data$class[degs_data$gene == 'wdfy2']

mean_expression_cluster_plot(obj,
                             'wdfy2',
                             3)

mean_expression_cluster_plot(obj,
                             'wdfy2',
                             6)

mean_expression_cluster_plot(obj,
                             'wdfy2',
                             15)
mean_expression_cluster_plot(obj,
                             'wdfy2',
                             19)

##### LOC111567884 #####

degs_data$cluster[degs_data$gene == 'LOC111567884']
degs_data$class[degs_data$gene == 'LOC111567884']
biomart_basic[biomart_basic$entrezgene_accession %in%'LOC111567884',]

mean_expression_cluster_plot(obj,
                             'LOC111567884',
                             0)

mean_expression_cluster_plot(obj,
                             'LOC111567884',
                             9)
mean_expression_cluster_plot(obj,
                             'LOC111567884',
                             23)

### LOC111565980 ####
degs_data$cluster[degs_data$gene == 'LOC111565980']
degs_data$class[degs_data$gene == 'LOC111565980']
biomart_basic[biomart_basic$entrezgene_accession %in%'LOC111565980',]

mean_expression_cluster_plot(obj,
                             'LOC111565980',
                             1)
mean_expression_cluster_plot(obj,
                             'LOC111565980',
                             4)

mean_expression_cluster_plot(obj,
                             'LOC111565980',
                             14)
### tcf7l2 ####
degs_data$cluster[degs_data$gene == 'tcf7l2']
degs_data$class[degs_data$gene == 'tcf7l2']
biomart_basic[biomart_basic$entrezgene_accession %in%'tcf7l2',]

mean_expression_cluster_plot(obj,
                             'tcf7l2',
                             17)
mean_expression_cluster_plot(obj,
                             'tcf7l2',
                             18)
###slc8a3 ####
degs_data$cluster[degs_data$gene == 'slc8a3']
degs_data$class[degs_data$gene == 'slc8a3']
biomart_basic[biomart_basic$entrezgene_accession %in%'slc8a3',]

mean_expression_cluster_plot(obj,
                             'slc8a3',
                             6)

mean_expression_cluster_plot(obj,
                             'slc8a3',
                             16)

### slc8a1b ###
degs_data$cluster[degs_data$gene == 'slc8a1b']
degs_data$class[degs_data$gene == 'slc8a1b']
biomart_basic[biomart_basic$entrezgene_accession %in%'slc8a1b',]

mean_expression_cluster_plot(obj,
                             'slc8a1b',
                             7)

mean_expression_cluster_plot(obj,
                             'slc8a1b',
                             8)


#### si:dkey-21e13.3 ####

degs_data$cluster[degs_data$gene == 'si:dkey-21e13.3']
degs_data$class[degs_data$gene == 'si:dkey-21e13.3']
biomart_basic[biomart_basic$entrezgene_accession %in%'si:dkey-21e13.3',]

mean_expression_cluster_plot(obj,
                             'si:dkey-21e13.3',
                             5)

mean_expression_cluster_plot(obj,
                             'si:dkey-21e13.3',
                             6)

### si:ch211-200p22.4   ####

degs_data$cluster[degs_data$gene == 'si:ch211-200p22.4']
degs_data$class[degs_data$gene == 'si:ch211-200p22.4']
biomart_basic[biomart_basic$entrezgene_accession %in%'si:ch211-200p22.4',]

mean_expression_cluster_plot(obj,
                             'si:ch211-200p22.4',
                             10)
mean_expression_cluster_plot(obj,
                             'si:ch211-200p22.4',
                             31)
###  prkcaa ####
degs_data$cluster[degs_data$gene == 'prkcaa']
degs_data$class[degs_data$gene == 'prkcaa']
biomart_basic[biomart_basic$entrezgene_accession %in%'prkcaa',]

mean_expression_cluster_plot(obj,
                             'prkcaa',
                             3)
mean_expression_cluster_plot(obj,
                             'prkcaa',
                             4)

### olfm1b ####

degs_data$cluster[degs_data$gene == 'olfm1b']
degs_data$class[degs_data$gene == 'olfm1b']
biomart_basic[biomart_basic$entrezgene_accession %in%'olfm1b',]

#ignoring this one cause cluster 30 doesnt count
 ### mpped2a ####
degs_data$cluster[degs_data$gene == 'mpped2a']
degs_data$class[degs_data$gene == 'mpped2a']
biomart_basic[biomart_basic$entrezgene_accession %in%'mpped2a',]

#ignoring cause 15 doesnt count

#### meis1b ####
degs_data$cluster[degs_data$gene == 'meis1b']
degs_data$class[degs_data$gene == 'meis1b']
biomart_basic[biomart_basic$entrezgene_accession %in%'meis1b',]

mean_expression_cluster_plot(obj,
                             'meis1b',
                             6)
mean_expression_cluster_plot(obj,
                             'meis1b',
                             22)

### LOC129350740 ####
degs_data$cluster[degs_data$gene == 'LOC129350740']
degs_data$class[degs_data$gene == 'LOC129350740']
biomart_basic[biomart_basic$entrezgene_accession %in%'LOC129350740',]

mean_expression_cluster_plot(obj,
                             'LOC129350740',
                             8)

mean_expression_cluster_plot(obj,
                             'LOC129350740',
                             9)

#### LOC129349260 #####
degs_data$cluster[degs_data$gene == 'LOC129349260']
degs_data$class[degs_data$gene == 'LOC129349260']
biomart_basic[biomart_basic$entrezgene_accession %in%'LOC129349260',]

mean_expression_cluster_plot(obj,
                             'LOC129349260',
                             6)

mean_expression_cluster_plot(obj,
                             'LOC129349260',
                             14)

### LOC129348729 ####
degs_data$cluster[degs_data$gene == 'LOC129348729']
degs_data$class[degs_data$gene == 'LOC129348729']
biomart_basic[biomart_basic$entrezgene_accession %in%'LOC129348729',]

mean_expression_cluster_plot(obj,
                             'LOC129348729',
                             2)

mean_expression_cluster_plot(obj,
                             'LOC129348729',
                             13)

### LOC129348286 ####
degs_data$cluster[degs_data$gene == 'LOC129348286']
degs_data$class[degs_data$gene == 'LOC129348286']
biomart_basic[biomart_basic$entrezgene_accession %in%'LOC129348286',]

#Im so tired of these uncharacterized LOCs I’m gonna just look at named genes now
### LOC111584813 ####
degs_data$cluster[degs_data$gene == 'LOC111584813']
degs_data$class[degs_data$gene == 'LOC111584813']
biomart_basic[biomart_basic$entrezgene_accession %in%'LOC111584813',]

#### LOC111584642 ####
degs_data$cluster[degs_data$gene == 'LOC111584642']
degs_data$class[degs_data$gene == 'LOC111584642']
biomart_basic[biomart_basic$entrezgene_accession %in%'LOC111584642',]

mean_expression_cluster_plot(obj,
                             'LOC111584642',
                             4)

mean_expression_cluster_plot(obj,
                             'LOC111584642',
                             6)
###LOC111568891 ###

degs_data$cluster[degs_data$gene == 'LOC111568891']
degs_data$class[degs_data$gene == 'LOC111568891']
biomart_basic[biomart_basic$entrezgene_accession %in%'LOC111568891',]

mean_expression_cluster_plot(obj,
                             'LOC111568891',
                             1)
mean_expression_cluster_plot(obj,
                             'LOC111568891',
                             8)

### LOC111567620 ####
degs_data$cluster[degs_data$gene == 'LOC111567620']
degs_data$class[degs_data$gene == 'LOC111567620']
biomart_basic[biomart_basic$entrezgene_accession %in%'LOC111567620',]

mean_expression_cluster_plot(obj,
                             'LOC111567620',
                             8)
mean_expression_cluster_plot(obj,
                             'LOC111567620',
                             19)
### LOC111564853 ####
degs_data$cluster[degs_data$gene == 'LOC111564853']
degs_data$class[degs_data$gene == 'LOC111564853']
biomart_basic[biomart_basic$entrezgene_accession %in%'LOC111564853',]

mean_expression_cluster_plot(obj,
                             'LOC111564853',
                             2)
mean_expression_cluster_plot(obj,
                             'LOC111564853',
                             8)
### kcnip3a ###
degs_data$cluster[degs_data$gene == 'kcnip3a']
degs_data$class[degs_data$gene == 'kcnip3a']
biomart_basic[biomart_basic$entrezgene_accession %in%'kcnip3a',]
#ignore cause 30 is fake

#### igf2bp1 ####
degs_data$cluster[degs_data$gene == 'igf2bp1']
degs_data$class[degs_data$gene == 'igf2bp1']
biomart_basic[biomart_basic$entrezgene_accession %in%'igf2bp1',]

mean_expression_cluster_plot(obj,
                             'igf2bp1',
                             0)
mean_expression_cluster_plot(obj,
                             'igf2bp1',
                             1)
 ### efna3b ####
degs_data$cluster[degs_data$gene == 'efna3b']
degs_data$class[degs_data$gene == 'efna3b']
biomart_basic[biomart_basic$entrezgene_accession %in%'efna3b',]


mean_expression_cluster_plot(obj,
                             'efna3b',
                             3)
mean_expression_cluster_plot(obj,
                             'efna3b',
                             18)

### cntn4 ####
degs_data$cluster[degs_data$gene == 'cntn4']
degs_data$class[degs_data$gene == 'cntn4']
biomart_basic[biomart_basic$entrezgene_accession %in%'cntn4',]
#15 is fake, move on

### chgb ####
degs_data$cluster[degs_data$gene == 'chgb']
degs_data$class[degs_data$gene == 'chgb']
biomart_basic[biomart_basic$entrezgene_accession %in%'chgb',]
#30 is fake

### cdh13 ####
degs_data$cluster[degs_data$gene == 'cdh13']
degs_data$class[degs_data$gene == 'cdh13']
biomart_basic[biomart_basic$entrezgene_accession %in%'cdh13',]

mean_expression_cluster_plot(obj,
                             'cdh13',
                             8)
mean_expression_cluster_plot(obj,
                             'cdh13',
                             29)

### brinp3a.1 ####
degs_data$cluster[degs_data$gene == 'brinp3a.1']
degs_data$class[degs_data$gene == 'brinp3a.1']
biomart_basic[biomart_basic$entrezgene_accession %in%'brinp3a.1',]

mean_expression_cluster_plot(obj,
                             'brinp3a.1',
                             2)
mean_expression_cluster_plot(obj,
                             'brinp3a.1',
                             6)

