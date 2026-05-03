library(Seurat)
library(patchwork)
library(tidyverse)
library(ggplot2)
library(Polychrome)
library(emmeans)
library(ggsignif)
clown_go = readRDS("Functions/clown_go2")
library(clusterProfiler)

obj = readRDS("~/Desktop/optimal_clustering_rna_only.rds")

status_to_phase = c(
  "D"  = 'I',
  'M'  = 'M',
  'F'  = 'F',
  'NF' = 'NF',
  'NRM'= 'NRM',
  'E'  = 'LI'
)

sub_1 = FindSubCluster(obj, 1, graph.name = 'harmony.wsnn')
Idents(sub_1) <- 'sub.cluster'
sub_1 = subset(sub_1, final_clusters == 1)
sub_1$Status = factor(sub_1$Status, levels = c('NRM', 'M', 'D', 'E', 'NF', 'F'))

sub_1@meta.data$Phase = status_to_phase[as.character(sub_1@meta.data$Status)]
sub_1@meta.data$Phase = factor(sub_1@meta.data$Phase,
                                levels = c('NRM', 'M', 'I', 'LI', 'NF', 'F'))

colors = c('#1965B0', '#4EB265', '#F7F056', '#DC050C')

p_value_annotation <- function(p) {
  if (p < 0.001) {
    return(list(label = "***", size = 7, adjust = 0.6))
  } else if (p < 0.01) {
    return(list(label = "**", size = 7, adjust = 0.6))
  } else if (p < 0.05) {
    return(list(label = "*", size = 7, adjust = 0.6))
  } else {
    return(list(label = paste0("p = ", round(p, 3)), size = 2, adjust = 0))
  }
}

degs     = read.csv('/Users/ggraham/Desktop/multiome_poa/DEG Outputs/FINAL degs classified w singular.csv')
deg_data = subset(degs, gene == 'LOC111577263' & cluster == 1)

d_m_fmt = p_value_annotation(deg_data$d_m_p.value)
d_f_fmt = p_value_annotation(deg_data$d_f_p.value)
f_m_fmt = p_value_annotation(deg_data$f_m_p.value)

cell_mask            = sub_1@meta.data$final_clusters == 1
meta_data            = sub_1@meta.data[cell_mask, ]
meta_data$expression = sub_1@assays$RNA$data['LOC111577263', cell_mask]

meta_grouped = meta_data %>%
  group_by(Phase, individual) %>%
  summarize(mean_expression = mean(expression), .groups = 'drop') %>%
  filter(Phase != 'NRM')

meta_grouped$Phase = factor(meta_grouped$Phase,
                             levels = c('M', 'I', 'LI', 'NF', 'F'))

y_max      = max(meta_grouped$mean_expression, na.rm = TRUE)
y_min      = min(meta_grouped$mean_expression, na.rm = TRUE) * 0.95
y_lower    = y_max * 1.10
y_upper    = y_max * 1.35
y_axis_max = y_max * 1.50

arom_plot = ggplot(meta_grouped, aes(x = Phase, y = mean_expression)) +
  geom_boxplot(data = subset(meta_grouped, Phase != 'NF'),
               aes(x = Phase, y = mean_expression, fill = Phase),
               outlier.shape = NA, width = 0.6) +
  geom_point(size = 1, position = position_jitter(width = 0.1, seed = 42)) +
  theme_classic() +
  scale_x_discrete(drop = FALSE) +
  scale_y_continuous(limits = c(y_min, y_axis_max)) +
  scale_fill_manual(values = colors) +
  labs(x = 'Phase', y = 'Normalized Expression', title = 'LOC111577263 | cluster 1') +
  theme(legend.position = 'none') +
  geom_signif(xmin = 1, xmax = 1.9,
              y_position = y_lower, annotation = d_m_fmt$label,
              color = 'black', tip_length = c(0, 0),
              textsize = d_m_fmt$size, vjust = d_m_fmt$adjust) +
  geom_signif(xmin = 2.1, xmax = 5,
              y_position = y_lower, annotation = d_f_fmt$label,
              color = 'black', tip_length = c(0, 0),
              textsize = d_f_fmt$size, vjust = d_f_fmt$adjust) +
  geom_signif(xmin = 1, xmax = 5,
              y_position = y_upper, annotation = f_m_fmt$label,
              color = 'black', tip_length = c(0, 0),
              textsize = f_m_fmt$size, vjust = f_m_fmt$adjust)

arom_plot

  
    ggsave(plot = arom_plot,
       file = paste0('aromatase','.svg'),
       device = "svg",
       units = "in",
       width = 2,
       height = 2,
       path = paste0("Manuscript/Plots/Manuscript v1.2/1_RGC_degs"))
