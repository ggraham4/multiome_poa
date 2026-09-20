#Negative binomial lower stringency 
{
  library(parallel)
  library(factoextra)
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


obj$cckb = ifelse(obj@meta.data$res0.8_50nn_40PC_45LSI == 6 & 
                    obj@assays$RNA$data['cckb',]>0, 
                  T,
                  F)

DimPlot(obj, group.by = 'cckb')


 cluster=T
 clustering = 'cckb'
 n_cores = detectCores() - 1
 
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
  Status <- factor(Status, levels = c("M", "D", "F"))  # pin level order so contrast
                                                        # naming doesn't depend on
                                                        # what ran earlier in the session
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
  
# Robust contrast extractor: works regardless of which order
# emmeans happened to name the contrast in (depends on current
# factor level order of Status, which can change between runs).
get_contrast <- function(pairs_res, level_a, level_b) {
  pdf <- as.data.frame(pairs_res)
  direct  <- pdf[pdf$contrast == paste(level_a, "-", level_b), ]
  reverse <- pdf[pdf$contrast == paste(level_b, "-", level_a), ]

  if (nrow(direct) == 1) {
    return(direct)
  } else if (nrow(reverse) == 1) {
    reverse$estimate <- -reverse$estimate   # A - B == -(B - A)
    if ("z.ratio" %in% names(reverse)) reverse$z.ratio <- -reverse$z.ratio
    return(reverse)
  } else {
    return(NULL)   # this pair doesn't exist -- a level was absent
  }
}

results <- mclapply(1:index, function(i) {
  dispersion <- dispersions.RAW[i]
  outcome <- df_counts_no_0_all_subjects[, i]

  out <- tryCatch({
    glmer_model <- suppressMessages(
      glmer(outcome ~ Status + (1 | subject),
            offset = offset,
            family = MASS::negative.binomial(theta = 1 / dispersion))
    )

    pairs_res <- pairs(emmeans(glmer_model, 'Status'), adjust = 'none')
    av <- car::Anova(glmer_model, type = 3)

    fm  <- get_contrast(pairs_res, "F", "M")
    dm  <- get_contrast(pairs_res, "D", "M")
    df_ <- get_contrast(pairs_res, "D", "F")

    data.frame(
      gene = valid_genes[i],
      f_m_estimate = if (!is.null(fm))  fm$estimate  else NA,
      f_m_p.value  = if (!is.null(fm))  fm$p.value   else NA,
      d_m_estimate = if (!is.null(dm))  dm$estimate  else NA,
      d_m_p.value  = if (!is.null(dm))  dm$p.value   else NA,
      d_f_estimate = if (!is.null(df_)) df_$estimate else NA,
      d_f_p.value  = if (!is.null(df_)) df_$p.value  else NA,
      av_p.value = av$`Pr(>Chisq)`[2],
      warning = ifelse(length(glmer_model@optinfo$conv$lme4$code) != 0,
                        substr(glmer_model@optinfo$conv$lme4$messages, 1, 50), NA),
      singular = isSingular(glmer_model)
    )
  }, error = function(e) {
    data.frame(gene = valid_genes[i], error_msg = conditionMessage(e))
  })

  out
}, mc.cores = n_cores)

  
results <- dplyr::bind_rows(results)
results <- as.data.frame(results, stringsAsFactors = FALSE)

  # genes that errored out will have an error_msg column populated and
  # NA f_m_p.value/d_m_p.value/d_f_p.value/av_p.value -- check these before
  # proceeding, e.g.: subset(results, !is.na(error_msg))

  results$f_m_p.value <- ifelse(results$f_m_p.value == 0, 1, results$f_m_p.value)
  results$d_m_p.value <- ifelse(results$d_m_p.value == 0, 1, results$d_m_p.value)
  results$d_f_p.value <- ifelse(results$d_f_p.value == 0, 1, results$d_f_p.value)
  
  results$av_q.value <- ifelse(test = is.na(results$warning),p.adjust(as.numeric(results$av_p.value), method = 'fdr',nrow(results)), "NA")
  

  message('Complete')
  end_time <- Sys.time()  # End timing
  message(end_time - start_time)  # Print the time difference


out = results
out$cluster ='cckb+ 6_poa_mixed'
write.csv(out, paste0('/Users/ggraham/Desktop/multiome_poa/DEG Outputs/2026_09_17 cckb 6_poa.csv'))


define_degs3 = define_degs_2 = function(data, alpha = 0.05){
  ref_data = read.csv("DEG Analyses/Expression DEGs/pairwise_patterns.csv")
  
  significant = subset(data, av_q.value < alpha)
  if(nrow(significant)<1){return(NULL)}
  significant$full_label = NA
  significant$first_word = NA
  significant$second_word = NA
  significant$short_label = NA

  newd <- data.frame()
  
  for(index in 1:nrow(significant)){
    data_of_interest = significant[index, ]
    
    if(data_of_interest$d_m_p.value < 0.05){
      ref_locs_d_m = which(ref_data$d_m != 'NS')
    } else{
      ref_locs_d_m = which(ref_data$d_m == 'NS')
    }
    
    if(data_of_interest$f_m_p.value < 0.05){
      ref_locs_f_m = which(ref_data$f_m != 'NS')
    } else{
      ref_locs_f_m = which(ref_data$f_m == 'NS')
    }
    
    if(data_of_interest$d_f_p.value < 0.05){
      ref_locs_d_f = which(ref_data$d_f != 'NS')
    } else{
      ref_locs_d_f = which(ref_data$d_f == 'NS')
    }
    
    locs_with_correct_signif_pattern = Reduce(intersect, list(ref_locs_d_f, ref_locs_f_m, ref_locs_d_m))
    
    if(data_of_interest$d_m_p.value < 0.05){
      if(data_of_interest$d_m_estimate < 0){
        dir_locs_d_m = which(ref_data$d_m == '<')
      } else{
        dir_locs_d_m = which(ref_data$d_m == '>')
      }
    } else{
      dir_locs_d_m = which(ref_data$d_m == 'NS')
    }
    
    if(data_of_interest$f_m_p.value < 0.05){
      if(data_of_interest$f_m_estimate < 0){
        dir_locs_f_m = which(ref_data$f_m == '<')
      } else{
        dir_locs_f_m = which(ref_data$f_m == '>')
      }
    } else{
      dir_locs_f_m = which(ref_data$f_m == 'NS')
    }
    
    if(data_of_interest$d_f_p.value < 0.05){
      if(data_of_interest$d_f_estimate < 0){
        dir_locs_d_f = which(ref_data$d_f == '<')
      } else{
        dir_locs_d_f = which(ref_data$d_f == '>')
      }
    } else{
      dir_locs_d_f = which(ref_data$d_f == 'NS')
    }
    
    dir_locs = Reduce(intersect, list(dir_locs_d_f, dir_locs_d_m, dir_locs_f_m))
    
    final_loc = intersect(dir_locs, locs_with_correct_signif_pattern)
    
    if(length(final_loc) > 0){
      full_label = ref_data[final_loc, ]$Full_label
      first_word = ref_data[final_loc, ]$First_word
      second_word = ref_data[final_loc, ]$Second_word
      short_label = ref_data[final_loc, ]$Short_label

      data_of_interest$full_label = full_label
      data_of_interest$first_word = first_word
      data_of_interest$second_word = second_word
      data_of_interest$short_label = short_label

      newd = rbind(newd, data_of_interest)
    }
  }
  return(newd)     
}

out_defined = define_degs3(out)
#gene_namer = readRDS('Functions/gene_namer.rds')
#out_defined$name = sapply(out_defined$gene, namer)

write.csv(out_defined, paste0('/Users/ggraham/Desktop/multiome_poa/DEG Outputs/defined 2026_09_17 cckb 6_poa.csv'))

out_defined = read.csv(paste0('/Users/ggraham/Desktop/multiome_poa/DEG Outputs/defined 2026_09_17 cckb 6_poa.csv'))
out = read.csv('/Users/ggraham/Desktop/multiome_poa/DEG Outputs/2026_09_17 cckb 6_poa.csv')