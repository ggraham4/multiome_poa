library(Seurat)
library(ggplot2)
library(ggsignif)

degs =read.csv('/Users/ggraham/Desktop/multiome_poa/DEG Outputs/FINAL degs classified w singular.csv')
degs_6 = subset(degs, cluster =='6')
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

colors = c('#1965B0', '#4EB265', '#F7F056', '#DC050C')



plotter_function_final = function(obj, deg, cluster_id,
                                  data     = NULL,
                                  d_m_p    = NULL,
                                  d_f_p    = NULL,
                                  f_m_p    = NULL,
                                  cluster_col = 'res0.8_50nn_40PC_45LSI') {

  p_value_annotation <- function(p) {
    if (p < 0.001) {
      return(list(label = '***', size = 7, adjust = 0.6))
    } else if (p < 0.01) {
      return(list(label = '**',  size = 7, adjust = 0.6))
    } else if (p < 0.05) {
      return(list(label = '*',   size = 7, adjust = 0.6))
    } else {
      return(list(label = paste0('p = ', round(p, 3)), size = 2, adjust = 0))
    }
  }

  colors = c('#1965B0', '#4EB265', '#F7F056', '#DC050C')

  # --- p-values: DEG table takes priority, else use supplied values ---
  if (!is.null(data)) {
    deg_data = subset(data, gene == deg & cluster == cluster_id)
    d_m_p = deg_data$d_m_p.value
    d_f_p = deg_data$d_f_p.value
    f_m_p = deg_data$f_m_p.value
  }

  get_fmt = function(p) {
    if (is.null(p) || is.na(p)) return(list(label = '', size = 3, adjust = 0))
    p_value_annotation(p)
  }

  d_m_fmt = get_fmt(d_m_p)
  d_f_fmt = get_fmt(d_f_p)
  f_m_fmt = get_fmt(f_m_p)

  # --- subset metadata & pull expression ---
  cell_mask = obj@meta.data[[cluster_col]] == cluster_id
  meta_data = obj@meta.data[cell_mask, ]
  meta_data$expression = obj@assays$RNA$data[deg, cell_mask]

  # --- per-individual means, drop NRM ---
  meta_grouped = meta_data %>%
    group_by(Phase, individual) %>%
    summarize(mean_expression = mean(expression), .groups = 'drop') %>%
    filter(Phase != 'NRM')

  meta_grouped$Phase = factor(meta_grouped$Phase,
                               levels = c('M', 'I', 'LI', 'NF', 'F'))

  # --- y-axis positions ---
  y_max      = max(meta_grouped$mean_expression, na.rm = TRUE)
  y_min      = min(meta_grouped$mean_expression, na.rm = TRUE) * 0.95
  y_lower    = y_max * 1.10
  y_upper    = y_max * 1.35
  y_axis_max = y_max * 1.50

  p = ggplot(meta_grouped, aes(x = Phase, y = mean_expression)) +
    geom_boxplot(data = subset(meta_grouped, Phase != 'NF'),
                 aes(x = Phase, y = mean_expression, fill = Phase),
                 outlier.shape = NA, width = 0.6) +
    geom_point(size = 1,
               position = position_jitter(width = 0.1, seed = 42)) +
    theme_classic() +
    scale_x_discrete(drop = FALSE) +
    scale_y_continuous(limits = c(y_min, y_axis_max)) +
    scale_fill_manual(values = colors) +
    labs(x     = 'Phase',
         y     = 'Normalized Expression',
         title = paste0(deg, '  |  cluster ', cluster_id)) +
    theme(legend.position = 'none')

  # --- brackets only when a p-value is available ---
  if (d_m_fmt$label != '') {
    p = p + geom_signif(xmin = 1, xmax = 1.9,
                        y_position = y_lower, annotation = d_m_fmt$label,
                        color = 'black', tip_length = c(0, 0),
                        textsize = d_m_fmt$size, vjust = d_m_fmt$adjust)
  }
  if (d_f_fmt$label != '') {
    p = p + geom_signif(xmin = 2.1, xmax = 5,
                        y_position = y_lower, annotation = d_f_fmt$label,
                        color = 'black', tip_length = c(0, 0),
                        textsize = d_f_fmt$size, vjust = d_f_fmt$adjust)
  }
  if (f_m_fmt$label != '') {
    p = p + geom_signif(xmin = 1, xmax = 5,
                        y_position = y_upper, annotation = f_m_fmt$label,
                        color = 'black', tip_length = c(0, 0),
                        textsize = f_m_fmt$size, vjust = f_m_fmt$adjust)
  }

  return(p)
}
saveRDS(plotter_function_final,'Functions/plotter_function_final.rds')

gnrh1_6 = plotter_function_final(obj, 'LOC111571064', 6)
 ggsave(plot = gnrh1_6,
       file = paste0('gnrh1_not_sig.svg'),
       device = "svg",
       units = "in",
       width = 2,
       height = 2,
       path = paste0("Manuscript/Plots/Manuscript v1.2/6_NPO_degs"))
 
 cckb_6 = plotter_function_final(obj, 'cckb', 6)
 ggsave(plot = cckb_6,
       file = paste0('cckb_not_sig.svg'),
       device = "svg",
       units = "in",
       width = 2,
       height = 2,
       path = paste0("Manuscript/Plots/Manuscript v1.2/6_NPO_degs"))
 
drd3_6 = plotter_function_final(obj, degs, 'drd3', cluster_id = 6)
drd3_6

for(i in degs_6$gene){
  p = plotter_function_final(obj, degs, i, cluster_id = 6)
  
    ggsave(plot = p,
       file = paste0(i,'.svg'),
       device = "svg",
       units = "in",
       width = 2,
       height = 2,
       path = paste0("Manuscript/Plots/Manuscript v1.2/6_NPO_degs"))
}

obj$cyto = CytoTRACE(obj@assays$RNA$data%>%as.matrix())$CytoTRACE

FeaturePlot(obj, 'cyto', reduction ='harmony_wnn.umap')


marks_22 = FindMarkers(obj, 22)
