#Negative binomial lower stringency updated
### here, I am going to do type III test first, THEN do FDR Adjustment, THEN
### do pairwise comparisons of those genes with tukey adjustment
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

obj <- readRDS('~/Desktop/optimal_clustering_rna_only.rds')


neg.bin.mult <- function(obj,
                         cluster, 
                         clustering = 'res0.8_50nn_40PC_45LSI',
                         n_cores = detectCores() - 1,
                         alpha = 0.05) {
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
    if(!exists('glmer_model')){
      glmer_model = NULL
    }
    
    if(!is.null(glmer_model)){
    

  av = car::Anova(glmer_model, type = 3)
    output_df <- data.frame(
      gene = valid_genes[i],
      av_p.value = av$`Pr(>Chisq)`[2],    
      warning = ifelse(length(glmer_model@optinfo$conv$lme4$code) != 0, substr(glmer_model@optinfo$conv$lme4$messages, 1, 50), NA),
      singular = ifelse(isSingular(glmer_model), TRUE, FALSE)
    )
    return(output_df)
    }
    
  }, mc.cores = n_cores
  )
  
  results_1 <- do.call(rbind, results)
  results_1 <- as.data.frame(results_1, stringsAsFactors = FALSE)
  results_1$av_q.value <- ifelse(
    is.na(results_1$warning),
    p.adjust(as.numeric(results_1$av_p.value),
             method = 'fdr',
             n = nrow(results_1)
             ),
    "NA")
  
  #define pairwise function to apply to a dataframe
  pairwise_function <- function(results_1){
    significant_genes <- results_1$gene[(as.numeric(results_1$av_q.value) < alpha) 
                                        & (!is.na(as.numeric(results_1$av_q.value)))]
          
           data = results_1
            data$f_m_p.value = NA
            data$d_m_p.value = NA
            data$d_f_p.value = NA
          
           #initialize columns in data frame
            
            for(j in significant_genes){
              
             message('Gene ', which(significant_genes ==j), ' of ', length(significant_genes))

              glmer_model <<- NULL
              
              i = which(valid_genes == j)
              
             #refit models
              dispersion <- dispersions.RAW[i]
              outcome <- df_counts_no_0_all_subjects[, i]
              
                  tryCatch({suppressMessages(
    glmer_model <- glmer(outcome ~ Status + (1 | subject),
                         offset = offset,
                         family = MASS::negative.binomial(theta = 1 / dispersion)))
     }, error = function(e) {
      return(NULL)
    })
    if(!exists('glmer_model')){
      glmer_model = NULL
    }
    
    if(!is.null(glmer_model)){
                  
                  pairs <- as.data.frame(pairs(emmeans(glmer_model, 'Status')))
                  
                  if(is.na(pairs$p.value[pairs$contrast == 'F - M'])){print(j)
                    break}

                  data$f_m_p.value[data$gene == j] = pairs$p.value[pairs$contrast == 'F - M']
                  data$d_m_p.value[data$gene == j] = pairs$p.value[pairs$contrast == 'D - M']
                  data$d_f_p.value[data$gene == j] = pairs$p.value[pairs$contrast == 'D - F']
    }
            }
            return(data)
          }
  
  full_results = pairwise_function(results_1)
  

  message('Complete')
  end_time <- Sys.time()  # End timing
  message(end_time - start_time)  # Print the time difference
  return(full_results)
}

for (i in 0:26) { ### already did 0
  print(i)
  output <- neg.bin.mult(obj,
                         cluster =i,
                         clustering = 'res0.8_50nn_40PC_45LSI')
  output$cluster = i
    write.csv(output, paste0('/Users/ggraham/Desktop/multiome_poa/DEG Outputs/05_31_2025 Neg Bin Anova First/cluster_',i,'.csv'))

}
# for some reason it stopped doing pairwise comparisons after 21?

### coalesce DEG list
final_data <- data.frame()
for (i in 0:26) { ### already did 0
  print(i)
    data = read.csv(paste0('/Users/ggraham/Desktop/multiome_poa/DEG Outputs/05_31_2025 Neg Bin Anova First/cluster_',i,'.csv'))
    subset_data = subset(data, av_q.value<0.05)
        if(nrow(subset_data)<1){next}
    final_data = rbind(final_data,subset_data )
}

#im still getting weird NAs in my pairwise comparisons so I rewrote the function
na_clusters <- unique(final_data$cluster[is.na(final_data$d_m_p.value)])
### even after I change the function I still get the problem here
for (i in na_clusters[2:length(na_clusters)]) {
  print(i)
  output <- neg.bin.mult(obj,
                         cluster =i,
                         clustering = 'res0.8_50nn_40PC_45LSI')
  output$cluster = i
    write.csv(output, paste0('/Users/ggraham/Desktop/multiome_poa/DEG Outputs/05_31_2025 Neg Bin Anova First/cluster_',i,'.csv'))

} ## and this doesnt fix it






