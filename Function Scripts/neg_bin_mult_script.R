neg_bin_mult <- function(obj,
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

saveRDS(neg_bin_mult, 'Functions/DEG_functions/neg_bin_mult.rds')
neg_bin_mult<- readRDS('Functions/DEG_functions/neg_bin_mult.rds')

