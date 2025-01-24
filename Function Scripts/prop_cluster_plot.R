  prop_cluster_plot <- function(object,
                         gene,
                         cluster,
                         clustering = 'harmony.wnn_res0.4_clusters'
){
  options(dplyr.summarise.inform = FALSE)


  counts <- suppressWarnings(t(as.matrix(object@assays$RNA$data[, object@meta.data[[clustering]] == cluster ]))[,gene])
  
  binary_counts <- (counts>0)+0
  
  n_cells <- nrow(binary_counts)
  n_genes <- ncol(binary_counts)
  
  

  #make meta data column
  full_data <- data.frame(
    cells = rownames(object@meta.data[object@meta.data[[clustering]] == cluster,]),
    Status = object@meta.data$Status[object@meta.data[[clustering]] == cluster],
    individual =object@meta.data$individual[object@meta.data[[clustering]] == cluster]
  )
    full_data$gene_expression <- binary_counts
    
    full_data <- full_data%>%
      group_by(individual, Status)%>%
    summarize(prop = sum(gene_expression)/n(),
              se = sqrt(prop*(1-prop)/n())
              )
  
    full_data$factor <- ifelse(full_data$Status == "NRM", 0, NA)
  full_data$factor <- ifelse(full_data$Status == "M", 1, full_data$factor)
  full_data$factor <- ifelse(full_data$Status == "D", 2, full_data$factor)
  full_data$factor <- ifelse(full_data$Status == "E", 3, full_data$factor)
  full_data$factor <- ifelse(full_data$Status == "NF", 4, full_data$factor)
  full_data$factor <- ifelse(full_data$Status == "F", 5, full_data$factor)

  full_data$Status <- factor(full_data$Status, levels <- c('NRM','M','D','E','NF',"F"))
  
  
  plot <- ggplot(full_data, aes(x = fct_reorder(individual, factor), y = prop, color = Status))+
    geom_boxplot(aes(group = Status, fill = Status), alpha = 0.25, outlier.shape = NA)+
    geom_point()+
    geom_pointrange(aes(x = individual, y = prop, ymin = prop - se, ymax = prop+se))+
    theme_classic()+
    labs(x  ='FishID', y = ' Proportion +/- SE ', title = paste0(gene,'_cluster_',cluster))+
    theme(axis.text.x = element_text(angle = -45, vjust = 1, hjust=0))
  return(plot)
  
  }
  
saveRDS(prop_cluster_plot, 'Functions/prop_cluster_plot.rds')
prop_cluster_plot<- readRDS( 'Functions/prop_cluster_plot.rds')
