#Negative binomial lower stringency GJG ###

# Library9es
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
}

# read in RNA object
obj <- readRDS('~/Desktop/optimal_clustering_rna_only.rds')


# define negative binomial regression function in parallel #
neg.bin.mult <- function(obj, # a seurat object
                         cluster, # which cluster to analyze
                         clustering = 'res0.8_50nn_40PC_45LSI', # the clustering variable
                         n_cores = detectCores() - 1) {
  start_time <- Sys.time()  # Start timing
  
  message('Extracting Counts')
  # extract counts matrix for only males, females, and dominants (intermediates) in that cluster
  counts <- obj@assays$RNA$counts[, obj@meta.data[[clustering]] == cluster & (obj@meta.data$Status == "M" | obj@meta.data$Status == "F" | obj@meta.data$Status == "D")]
  combined_counts <- counts # rename the variable cause this is someone elses code I stole
  
  # convert to df and readd gene names as colnames
  df_counts <- data.frame(t(combined_counts))
  colnames(df_counts) <- rownames(obj@assays$RNA)
  
  # number of genes and cells for offset calculation 
  n_genes = ncol(df_counts)
  n_cells = nrow(df_counts)
  
  message("Making Counts Data Frame...")
  df_counts_meta <- data.frame(rownames(df_counts)) # make a dataframe of the genes
  df_counts_meta$id <- df_counts_meta$rownames.df_counts. # of the add into the df_counts dataframe as a column id
  df_counts_meta$rownames.df_counts. = NULL # get rid of that original column for some reason 
  
  # add individual ID and sex (Status) into the df
  df_counts_meta$individual = obj$individual[obj@meta.data[[clustering]] == cluster & (obj@meta.data$Status == "M" | obj@meta.data$Status == "F" | obj@meta.data$Status == "D")]
  df_counts_meta$Status = obj$Status[obj@meta.data[[clustering]] == cluster & (obj@meta.data$Status == "M" | obj@meta.data$Status == "F" | obj@meta.data$Status == "D")]
  
  # remove genes with 0 counts
  message("Removing Genes with 0 Counts...")
  df_counts_no_0 <- df_counts[, which(colSums(df_counts) != 0)]
  
  message("Making New Counts Data Frame Without 0s...")
  n_genes_no_0 = ncol(df_counts_no_0) # count non zero genes only
  
  # make a new df of the non zero genes and the meta data
  df_counts_no_0 <- cbind(df_counts_no_0, df_counts_meta)
  #split by individual
  df_counts_no_0_split_by_subject <- split(df_counts_no_0, f = df_counts_no_0$individual)
  

  # if you are reading these comments I sincerely apologize for this part
  # lynna's code kept messing up here so the only way I was able to get it to work
  # was to not filter each gene by if its present in each individual, 
  # that's why this file is called lower stringency 
  # basically that filtration just doesnt happen anymore
  
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
    out <- data.frame(colnames(temp_subject_l_counts_no_0))
    assign(x = paste0("gene_list_subject_", l), value = get("out"))
  }
  
  # GENERATE A LIST OF GENES FOR EACH SUBJECT (these should all be identical)
  good_gene_list <- gene_list_subject_1$colnames.temp_subject_l_counts_no_0.
  # intersect gene lists for each subject, again, should all be the same
  for (m in 2:length(df_counts_no_0_split_by_subject)) {
    temp_good_gene_list_m <- data.frame(value = get(paste0("gene_list_subject_", m)))
    temp_good_gene_list_m <- temp_good_gene_list_m$colnames.temp_subject_l_counts_no_0.
    good_gene_list <- intersect(good_gene_list, temp_good_gene_list_m)
  }
  
  #the number of good genes
  p <- length(good_gene_list)
  
  message('Making Gene Data Frame for Each Subject...')
  # genes that are non-zero and present in all individuals
  valid_genes <- good_gene_list[good_gene_list %in% colnames(df_counts_no_0)]
  # number of valid genes
  v <- length(valid_genes)
  
  # subset to valid genes
  df_counts_no_0_all_subjects <- df_counts_no_0[, valid_genes]
  # make a matrix of that df
  count_matrix_final <- as.matrix(df_counts_no_0_all_subjects)
  
  # and then transpose it such that genes are now rows and cells are columns this is what size factors expects
  count_matrix_final <- as.data.frame(t(count_matrix_final))
  
  #add meta data back into
  df_counts_no_0_all_subjects <- cbind(df_counts_no_0_all_subjects, df_counts_meta)
  
  # ensure status and subeject are factors
  Status <- as.factor(df_counts_no_0_all_subjects$Status)
  subject <- as.factor(df_counts_no_0_all_subjects$individual)
  
  message("Estimating Dispersion Using Gamma-Poisson...")
  # cluster size is the number of cells present in the cluster
  cluster_size <- ncol(count_matrix_final)
  
  # what calcualteSumFactors is doing is scaling and normalizing 
  # the data by deconvolution,
  #> "Briefly, a pool of cells is selected and the expression profiles for
  #>  those cells are summed together. The pooled expression profile is 
  #>  normalized against an average reference pseudo-cell, constructed by
  #>   averaging the counts across all cells. This defines a size factor for
  #>    the pool as the median ratio between the count sums and the average 
  #>    across all genes.
  #>    
  #>    https://rdrr.io/bioc/scran/man/computeSumFactors.html
  
  size_factors <- calculateSumFactors(count_matrix_final,
                                      clusters = NULL,
                                      ref.clust = NULL,
                                      max.cluster.size = cluster_size,
                                      positive = TRUE,
                                      scaling = NULL,
                                      min.mean = NULL,
                                      subset.row = NULL)
  
  # make a data frame of status
  coldata <- data.frame(Status)
  
  #here, dispersions are being calculated using a gamma poisson glm
  fit <- glm_gp(as.matrix(count_matrix_final), # the matrix
                col_data = coldata,  # status 
                size_factors = size_factors, # size factors corrects for different counts based on sequencing depth
                design = ~ Status, # for each gene, correct for status
                on_disk = FALSE)
  # extract dispersions
  dispersions.RAW <- fit$overdispersion_shrinkage_list$ql_disp_estimate
  # log transform
  log.sizeFactors.RAW <- log(size_factors)
  
  # some more renaming for some reason
  offset <- log.sizeFactors.RAW
  index <- v
  
  results <- mclapply(1:index, function(i) { 
    # in parallel, for each gene we are running a mixed negative binomial regression 
    # with theta being 1/ dispersion, and a random effect of individual
    # and an for each cell from log.sizeFactors
    
    message('Calculating Gene', paste0(i, ' of ', index, '...'))
    # extract dispersion and gene counts
    dispersion <- dispersions.RAW[i]
    outcome <- df_counts_no_0_all_subjects[, i]
    
    
    # try catch cause these models throw a lot of errors
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
    
      # once the model is run, compute pairwise comparisons
      # using emmeans with no p.value adjustment
    pairs <- pairs(emmeans(glmer_model, 'Status'), adjust = 'none')
    
    # calculate type 3 tests using car:: Anova, this is what will get
    # the fdr adjustment not the pairwise comparisons 
    
  av = car::Anova(glmer_model, type = 3)
  
  # add to df
    output_df <- data.frame(
      gene = valid_genes[i],
      f_m_estimate = as.data.frame(pairs[pairs@grid$contrast == 'F - M'])[, 2],
      f_m_p.value = as.data.frame(pairs[pairs@grid$contrast == 'F - M'])[, 6],
      d_m_estimate = as.data.frame(pairs[pairs@grid$contrast == 'D - M'])[, 2],
      d_m_p.value = as.data.frame(pairs[pairs@grid$contrast == 'D - M'])[, 6],
      d_f_estimate = as.data.frame(pairs[pairs@grid$contrast == 'D - F'])[, 2],
      d_f_p.value = as.data.frame(pairs[pairs@grid$contrast == 'D - F'])[, 6],
      av_p.value = av$`Pr(>Chisq)`[2],    
      
      # extract any warnings from the model
      warning = ifelse(length(glmer_model@optinfo$conv$lme4$code) != 0, substr(glmer_model@optinfo$conv$lme4$messages, 1, 50), NA),
      #record if the model was singular
      singular = ifelse(isSingular(glmer_model), TRUE, FALSE)
    )
    return(output_df)
    }
    
  }, mc.cores = n_cores
  )
  
  # bind all results
  results <- do.call(rbind, results)
  results <- as.data.frame(results, stringsAsFactors = FALSE)
  
  # if any p values = 0, I assume the model failed and assign them 1 
  results$f_m_p.value <- ifelse(results$f_m_p.value == 0, 1, results$f_m_p.value)
  results$d_m_p.value <- ifelse(results$d_m_p.value == 0, 1, results$d_m_p.value)
  results$d_f_p.value <- ifelse(results$d_f_p.value == 0, 1, results$d_f_p.value)
  
  # fdr adjust anova p values by how many genes were analyzed
  results$av_q.value <- ifelse(is.na(results$warning), p.adjust(as.numeric(results$av_p.value), method = 'fdr'), n = nrow(results), "NA")
  
  #assign df to global environment
  assign(paste0("results_", "cluster", cluster), results, envir = .GlobalEnv)
  
  message('Complete')
  end_time <- Sys.time()  # End timing
  message(end_time - start_time)  # Print the time difference
  return(results)
}

# bind 
data = data.frame()
for (i in 0:26) {
  print(i)
  output <- neg.bin.mult(obj,
                         cluster =i,
                         clustering = 'res0.8_50nn_40PC_45LSI')
  output$cluster = i
  output <- rbind(output, data)
}

# output 
data <- data.frame()
for(i in 0:26){
  results <- get(paste0('results_cluster',i))
    results$f_m_q.value <- ifelse(is.na(results$warning), p.adjust(as.numeric(results$f_m_p.value), method = "fdr", n = nrow(results)), "NA")
  results$d_m_q.value <- ifelse(is.na(results$warning), p.adjust(as.numeric(results$d_m_p.value), method = "fdr", n = nrow(results)), "NA")
  results$d_f_q.value <- ifelse(is.na(results$warning), p.adjust(as.numeric(results$d_f_p.value), method = "fdr", n = nrow(results)), "NA")
  results$av_q.value <- ifelse(is.na(results$warning), p.adjust(as.numeric(results$av_p.value), method = 'fdr', n = nrow(results)), "NA")
  results$av_q.value <- ifelse(is.na(results$warning), p.adjust(p = results$av_p.value, 'fdr', nrow(results)),'NA')
  results$cluster = i
  data = rbind(data, results)
    write.csv(results, paste0('/Users/ggraham/Desktop/multiome_poa/DEG Outputs/05_11_2025 Neg Bin w Doms New_clustering/cluster_',i,'.csv'))
  
}

  #write.csv(data, paste0('/Users/ggraham/Desktop/multiome_poa/DEG Outputs/05_11_2025 Neg Bin w Doms New_clustering_all_clusters.csv'))

