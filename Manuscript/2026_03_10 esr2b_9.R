library(Seurat)
library(ggplot2)
library(ggsignif)

degs =read.csv('/Users/ggraham/Desktop/multiome_poa/DEG Outputs/FINAL degs classified w singular.csv')
degs_esr2b = subset(degs, gene =='esr2b')
obj  = readRDS("~/Desktop/optimal_clustering_rna_only.rds")
obj$Status = factor(obj$Status, levels = c('NRM','M',"D",'E','NF','F'))

obj@meta.data$Phase = as.character(obj@meta.data$Status)
obj@meta.data$Phase= ifelse(obj@meta.data$Phase == 'D', 'I', obj@meta.data$Phase)
obj@meta.data$Phase= ifelse(obj@meta.data$Phase == 'E', 'LI', obj@meta.data$Phase)

obj@meta.data$Phase = factor(obj@meta.data$Phase, levels = c('NRM',
                                                             'M',
                                                             'I',
                                                             'LI',
                                                             'NF',
                                                             'F'))


p_value_annotation <- function(p) {
  if (p < 0.001) {
    return("***")
  } else if (p < 0.01) {
    return("**")
  } else if (p < 0.05) {
    return("*")
  } else {
    return(paste0("p = ", round(p, 3)))
  }
}

plotter_function_final = function(obj,data, deg){
  
  deg_data = subset(data, gene == deg & cluster ==9)
  d_f_p.value = deg_data$d_f_p.value
  d_m_p.value = deg_data$d_m_p.value
  f_m_p.value = deg_data$f_m_p.value
  
  d_f_annotation = p_value_annotation(d_f_p.value)
  d_m_annotation = p_value_annotation(d_m_p.value)
  f_m_annotation = p_value_annotation(f_m_p.value)
  
  meta.data = subset(obj@meta.data, final_clusters ==9)
  meta.data$expression = obj@assays$RNA$data[deg,obj@meta.data$final_clusters == 9]
  
  meta_grouped = meta.data%>%
    group_by(Phase, individual)%>%
    summarize(mean_expression = mean(expression))%>%
    subset(Phase != 'NRM')
  
  meta_grouped$Phase = factor(meta_grouped$Phase, levels = c(
                                                             'M',
                                                             'I',
                                                             'LI',
                                                             'NF',
                                                             'F'))

  
    r = ggplot(meta_grouped, aes(x = Phase, y = mean_expression))+
    geom_boxplot(data = subset(meta_grouped, Phase != "NF"),
                 aes(x = Phase, y =mean_expression ,
                     fill = Phase), outlier.shape = NA)+
    geom_point(size =1)+
    theme_classic()+
        scale_x_discrete(drop = FALSE)+  
    labs(x = 'Phase', y = 'Normalized Expression')+
    geom_signif(xmin = c(1),
                xmax = c(1.9),
                y_position = c(max(meta_grouped$mean_expression*1.1)),
                annotation =c(d_m_annotation),
                color = "black",
                tip_length = c(0,0), textsize = 3)+
          geom_signif(xmin = c(1),
                xmax = c(5),
                y_position = c(max(meta_grouped$mean_expression*1.3)),
                annotation =c(f_m_annotation),
                color = "black",
                tip_length = c(0,0), textsize = 3)+
    geom_signif(xmin = c(2.1),
                xmax = c(5),
                y_position = c(max(meta_grouped$mean_expression*1.1)),
                annotation =c(d_f_annotation),
                color = "black",
                tip_length = c(0,0), textsize = 3)+
      theme(legend.position = 'none',)
    r
  return(r)
}



esr2b_9 = plotter_function_final(obj, degs_esr2b, 'esr2b')+
  labs( y = 'mean expression')
esr2b_9

  ggsave(plot = esr2b_9,
       file = 'esr2b.svg',
       device = "svg",
       units = "in",
       width = 2.5,
       height = 2.5,
       path = "Manuscript/Plots/")

