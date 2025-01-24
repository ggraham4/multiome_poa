prop_deg_function <- function(object = obj, cluster = 19, clustering = 'harmony.wnn_res0.4_clusters'){
  start <- Sys.time()
  library(lme4)
  library(dplyr)
  library(parallel)
  library(car)
  
  options(dplyr.summarise.inform = FALSE)
  
  #extract counts matrix 
  counts <- t(as.matrix(object@assays$RNA$data[, object@meta.data[[clustering]] == cluster & (object@meta.data$Status == "M" | object@meta.data$Status == "F" | object@meta.data$Status == "D")]))
  
  binary_counts <- (counts>0)+0
  
  n_cells <- nrow(binary_counts)
  n_genes <- ncol(binary_counts)
  
  #convert to matrix
  counts_matrix <- matrix(data = binary_counts, 
                          nrow =n_cells,
                          ncol= n_genes)
  
  rownames(counts_matrix) <- rownames(binary_counts)
  colnames(counts_matrix) <- colnames(binary_counts)
  
  ##Remove 0 genes 
  filtered_cols_matrix <- counts_matrix[,which(colSums(counts_matrix) != 0)]
  
  #make meta data column
  meta_data <- data.frame(
    cells = rownames(object@meta.data[object@meta.data[[clustering]] == cluster & (object@meta.data$Status == "M" | object@meta.data$Status == "F" | object@meta.data$Status == "D"),]),
    Status = object@meta.data$Status[object@meta.data[[clustering]] == cluster & (object@meta.data$Status == "M" | object@meta.data$Status == "F" | object@meta.data$Status == "D")],
    individual =object@meta.data$individual[object@meta.data[[clustering]] == cluster & (object@meta.data$Status == "M" | object@meta.data$Status == "F" | object@meta.data$Status == "D")]
  )
  
  
  n_cells_individual <- meta_data%>%
    group_by(individual)%>%
    summarize(n_cells = n())
  
  genes <- colnames(filtered_cols_matrix)
  genes <- genes[!is.na(genes)]
  
  
  deg_function <- function(gene){
    message(paste0('gene ', which(genes == gene), ' of ', n_genes))
    
    gene_expression <- filtered_cols_matrix[, gene, drop = FALSE]
    
    meta_data$gene <- gene_expression
    
    joined_data_restrictions <- meta_data%>%
      group_by(Status)%>%
      summarize(cells_expressing = sum(gene))
    
    if (any(joined_data_restrictions$cells_expressing < 1)) return(NULL)
    
    data_for_matrix <- meta_data %>%
      group_by(individual, Status) %>%
      summarise(
        cells_expressing = sum(gene),
        cells_total = n()
      )
    
    model_matrix <- with(data_for_matrix, cbind(cells_expressing, cells_total - cells_expressing))
    
    
    if (length(unique(model_matrix[, 1])) == 1) return(NULL)
    
    
    glmer_model <- tryCatch({suppressMessages(glmer(model_matrix~Status + (1|individual), family = binomial('logit'), data = data_for_matrix))
    }, error = function(e) {
    return(NULL)
  })
    if (is.null(glmer_model)) return(NULL)
    
    pairwise_comps <- as.data.frame(pairs(emmeans(glmer_model, 'Status'), adjust ='none'))
    
    data_for_output <- data_for_matrix %>%
      group_by(Status) %>%
      summarise(
        prop_expressing = sum(cells_expressing) / sum(cells_total),
        cells_expressing = sum(cells_expressing),
        cells_total = sum(cells_total)
      ) %>%
      pivot_wider(names_from = 'Status', values_from = c(cells_expressing,cells_total, prop_expressing))
    
    data_for_output$anova_p.value <- car::Anova(glmer_model, type = 'III')[2,3]
    data_for_output$d_f_p.value <- pairwise_comps$p.value[pairwise_comps$contrast=='D - F']
    data_for_output$d_m_p.value <- pairwise_comps$p.value[pairwise_comps$contrast=='D - M']
    data_for_output$f_m_p.value <- pairwise_comps$p.value[pairwise_comps$contrast=='F - M']
    data_for_output$singular = isSingular(glmer_model)
    data_for_output$warning = ifelse(length(glmer_model@optinfo$conv$lme4$code) != 0, substr(glmer_model@optinfo$conv$lme4$messages, 1, 50), NA)
    data_for_output$gene <- gene
    data_for_output$cluster <- cluster
    
    return(data_for_output)
    
  }
  
  library(parallelsugar)
  deg_output <- parallelsugar::mclapply(X=genes, FUN=deg_function
                       ,mc.cores= detectCores()-1
  )
  
  deg_output2 <- do.call(rbind, deg_output)
  
  deg_output2$anova_q.value <- p.adjust(deg_output2$anova_p.value, 'fdr', nrow(deg_output2))
  
  deg_output2$d_f_q.value <- p.adjust(deg_output2$d_f_p.value, 'fdr', nrow(deg_output2))
  deg_output2$d_m_q.value <- p.adjust(deg_output2$d_m_p.value, 'fdr', nrow(deg_output2))
  deg_output2$f_m_q.value <- p.adjust(deg_output2$f_m_p.value, 'fdr', nrow(deg_output2))
  
  end <- Sys.time()
  print(end-start)
  return(deg_output2)
  
}

saveRDS(prop_deg_function, 'Functions/DEG_functions/prop_deg_function.rds')
prop_deg_function.rds<- readRDS('Functions/DEG_functions/prop_deg_function.rds')
