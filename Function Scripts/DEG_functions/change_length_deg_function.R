change_length_function <- function(object = obj, 
                              cluster = 19,
                              clustering = 'harmony.wnn_res0.4_clusters'){
  start <- Sys.time()
  library(lme4)
  library(dplyr)
  #library(parallelsugar)
  
  options(dplyr.summarise.inform = FALSE)
  
  #extract counts matrix 
  counts <- t(as.matrix(object@assays$RNA$data[, object@meta.data[[clustering]] == cluster & object@meta.data$Status == 'D']))
  
  filtered_cols_matrix <- counts[,which(colSums(counts) != 0)]
  
  meta_data <- data.frame(
    cells = rownames(object@meta.data[object@meta.data[[clustering]] == cluster & object@meta.data$Status == 'D',]),
    Status = object@meta.data$Status[object@meta.data[[clustering]] == cluster & object@meta.data$Status == 'D'],
    individual =object@meta.data$individual[object@meta.data[[clustering]] == cluster & object@meta.data$Status == 'D'],
    change_length = object@meta.data$Change_Length[object@meta.data[[clustering]] == cluster &object@meta.data$Status == 'D']
  )
  
  genes <- colnames(filtered_cols_matrix)
  genes <- genes[!is.na(genes)]
  
  n_genes <- length(genes)
  
  deg_function <- function(gene){
    
    message(paste0('gene ', which(genes == gene), ' of ', n_genes))
    
    gene_expression <- filtered_cols_matrix[, gene]
    meta_data$gene <- gene_expression
    
    joined_data_restrictions <- meta_data%>%
      group_by(Status)%>%
      summarize(expression = sum(gene))
    
    if (any(joined_data_restrictions$expression < 1)) return(NULL)
    
    joined_data_for_model <- meta_data%>%
      group_by(Status, individual)%>%
      summarize(expression = mean(gene),
                se_expression = sd(gene)/sqrt(n()),
                change_length = mean(change_length)
      )%>%na.omit()
    
    if (length(unique(joined_data_for_model$expression)) == 1) return(NULL)
    
    ###length model
    length_model <- lm(expression~change_length, data=joined_data_for_model)
    
    #summary
    length_summary <- summary(length_model)
    length_coefs <- as.data.frame(length_summary$coefficients)
    
    #values 
    length_estimate <- length_coefs[2,1]
    length_summary_p.value <- length_coefs[2,4]

    
    data_for_output <- data.frame(cluster = cluster,
                                  gene = gene,
                                  length_estimate=length_estimate,
                                  length_summary_p.value=length_summary_p.value
    )
    return(data_for_output)
  }
  
  #deg_output <- parallelsugar::mclapply(X=genes, FUN=deg_function,mc.cores= detectCores()-1)
  deg_output <- lapply(X=genes, FUN=deg_function)
  
  deg_output2 <- do.call(rbind, deg_output)
  
  deg_output2$length_summary_q.value <- p.adjust(deg_output2$length_summary_p.value, 'fdr', nrow(deg_output2))

  
  
  return(deg_output2)
  
}