mean_expression_cluster_data <- function(object, gene, cluster, clustering = 'harmony.wnn_res0.4_clusters'){
      options(dplyr.summarise.inform = FALSE)

    counts <- t(object@assays$RNA$data[,object@meta.data[[clustering]] == cluster])
  Counts_of_interest <- as.data.frame(counts[,gene])
    Counts_of_interest$expression <- Counts_of_interest[,1]
  Counts_of_interest$individual <- object@meta.data$individual[object@meta.data[[clustering]] == cluster]
    Counts_of_interest$Status <- object@meta.data$Status[object@meta.data[[clustering]] == cluster]

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
  return(results)
}

saveRDS(mean_expression_cluster_data, 'Functions/mean_expression_cluster_data.rds')
mean_expression_cluster_data<- readRDS('Functions/mean_expression_cluster_data.rds')
