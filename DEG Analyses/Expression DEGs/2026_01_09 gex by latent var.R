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



obj <- readRDS('~/Desktop/optimal_clustering_rna_only.rds')

measures = read.csv('2025_12_26 all_data.csv')
measures$individual = measures$Fish
measures_sexstate =dplyr::select(measures,   c('individual', 'SexState'))

obj@meta.data = obj@meta.data%>%
  left_join(measures_sexstate, by = 'individual')

neg.bin.mult <- function(obj,
                         cluster, 
                         clustering = 'res0.8_50nn_40PC_45LSI',
                         n_cores = detectCores() - 1) {
  
  start_time <- Sys.time()  # Start timing
  
  subset_indices <- obj@meta.data[[clustering]] == cluster & 
                    (obj@meta.data$Status != 'NRM')
  
  message('Extracting Counts')
  counts <- obj@assays$RNA$counts[, subset_indices]
  combined_counts <- counts
  
  df_counts <- data.frame(t(combined_counts))
  colnames(df_counts) <- rownames(obj@assays$RNA)
  
  n_genes = ncol(df_counts)
  n_cells = nrow(df_counts)
  
  message("Making Counts Data Frame...")
  df_counts_meta <- data.frame(rownames(df_counts))
  df_counts_meta$id <- df_counts_meta$rownames.df_counts.
  df_counts_meta$rownames.df_counts. = NULL
  
  df_counts_meta$individual = obj$individual[subset_indices]

    df_counts_meta$SexState = as.numeric(obj$SexState[subset_indices]) 
  
  message("Removing Genes with 0 Counts...")
  df_counts_no_0 <- df_counts[, which(colSums(df_counts) != 0)]
  
  message("Making New Counts Data Frame Without 0s...")
  n_genes_no_0 = ncol(df_counts_no_0)
  
  
  df_counts_no_0 <- cbind(df_counts_no_0, df_counts_meta)
  df_counts_no_0_split_by_subject <- split(df_counts_no_0, f = df_counts_no_0$individual)
  
  message("Finding Good Genes for Subject...")
  
  for (l in 1:length(df_counts_no_0_split_by_subject)) {
    correct_gene_names <- colnames(df_counts_no_0)
    
    temp_subject_l <- data.frame(df_counts_no_0_split_by_subject[[l]]) 
    colnames(temp_subject_l) <- correct_gene_names
    
    temp_subject_l_counts <- temp_subject_l[, 1:n_genes_no_0]
    
    temp_subject_l_counts_no_0 <- temp_subject_l_counts 
    
    out <- data.frame(colnames(temp_subject_l_counts_no_0))
    assign(x = paste0("gene_list_subject_", l), value = get("out"))
  }
  
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
  
  SexState <- as.numeric(df_counts_no_0_all_subjects$SexState)
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
  
  coldata <- data.frame(SexState = SexState)
  
  # Calculate Dispersion (Using continuous design)
  fit <- glm_gp(as.matrix(count_matrix_final), 
                col_data = coldata, 
                size_factors = size_factors, 
                design = ~ SexState, 
                on_disk = FALSE)
                
  dispersions.RAW <- fit$overdispersion_shrinkage_list$ql_disp_estimate
  log.sizeFactors.RAW <- log(size_factors)
  
  offset <- log.sizeFactors.RAW
  index <- v
  
  results <- mclapply(1:index, function(i) {
    if(i %% 100 == 0) message('Calculating Gene ', paste0(i, ' of ', index, '...'))
    
    dispersion <- dispersions.RAW[i]
    outcome <- df_counts_no_0_all_subjects[, i]
    
    tryCatch({suppressMessages(
      glmer_model <- glmer(outcome ~ SexState + (1 | subject),
                           offset = offset,
                           family = MASS::negative.binomial(theta = 1 / dispersion),
                           control = glmerControl(optimizer = "bobyqa",
                                                  optCtrl = list(maxfun = 2e5))))
    }, error = function(e) {
      return(NULL)
    })
    
    if(!exists('glmer_model')){
      glmer_model = NULL
    }
    
    if(!is.null(glmer_model)){
      
      sum = (summary(glmer_model)$coefficients) %>% as.data.frame()
      av = car::Anova(glmer_model, type=3) %>% as.data.frame()
      
      output_df <- data.frame(
        gene = valid_genes[i],
        # sum[2,1] targets the Estimate for SexState (row 2, col 1)
        estimate = sum[2,1], 
        p.value = av$`Pr(>Chisq)`[2],    
        warning = ifelse(length(glmer_model@optinfo$conv$lme4$code) != 0, substr(glmer_model@optinfo$conv$lme4$messages, 1, 50), NA),
        singular = ifelse(isSingular(glmer_model), TRUE, FALSE)
      )
      return(output_df)
    }
    
  }, mc.cores = n_cores
  )
  
  results <- do.call(rbind, results)
  results <- as.data.frame(results, stringsAsFactors = FALSE)
  results$p.value <- ifelse(results$p.value == 0, 1, results$p.value)
  results$q.value <- ifelse(is.na(results$warning), p.adjust(as.numeric(results$p.value), method = 'fdr', n = nrow(results)), "NA")
  
  assign(paste0("results_", "cluster", cluster), results, envir = .GlobalEnv)
  
  message('Complete')
  end_time <- Sys.time()  # End timing
  message(end_time - start_time)  # Print the time difference
  return(results)
}

data = data.frame()
for (i in 26:0) {
  print(i)
  output <- neg.bin.mult(obj,
                         cluster =i,
                         clustering = 'res0.8_50nn_40PC_45LSI')
  output$cluster = i
  output <- rbind(output, data)
}

data <- data.frame()
for(i in 0:26){
  results <- get(paste0('results_cluster',i))
  results$q.value <- ifelse(is.na(results$warning), p.adjust(p = results$p.value, 'fdr', nrow(results)),'NA')
  results$cluster = i
  data = rbind(data, results)
    write.csv(results, paste0('/Users/ggraham/Desktop/multiome_poa/DEG Outputs/2026_01_09 Neg Bin latent/cluster_',i,'.csv'))
  
}

write.csv(data, paste0('/Users/ggraham/Desktop/multiome_poa/DEG Outputs/2026_01_09 Neg Bin latent_all_clusters.csv'))

