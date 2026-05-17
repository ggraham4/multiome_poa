library(Seurat)
library(patchwork)
library(tidyverse)
library(ggplot2)
library(Polychrome)
library(emmeans)
library(ggsignif)
library(clusterProfiler)

obj = readRDS("~/Desktop/optimal_clustering_rna_only.rds")

sub_1 = FindSubCluster(obj, 1, graph.name = 'harmony.wsnn')
Idents(sub_1) <- 'sub.cluster'
sub_1 = subset(sub_1, final_clusters == 1)
sub_1$Status = factor(sub_1$Status, levels = c('NRM', 'M', 'D', 'E', 'NF', 'F'))

dimplot = DimPlot(sub_1, label = F) +
  theme_void() +
  theme(legend.position = 'none')
dimplot
DimPlot(sub_1, label = T)

# ---- shared helpers ----
colors = c('#1965B0', '#4EB265', '#F7F056','#7BAFDE', '#DC050C')

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

# ---- cell count data ----
cells_ind = sub_1@meta.data %>%
  group_by(individual) %>%
  summarize(n_cells = n())

cells_sub_ind = sub_1@meta.data %>%
  group_by(individual, Status, sub.cluster) %>%
  summarize(n_cells_in = n(), .groups = 'drop')

cells_total = cells_ind %>%
  right_join(cells_sub_ind, by = 'individual')

cells_total$Phase = as.character(cells_total$Status)
cells_total$Phase = ifelse(cells_total$Phase == 'D',  'I',  cells_total$Phase)
cells_total$Phase = ifelse(cells_total$Phase == 'E',  'LI', cells_total$Phase)
cells_total$Phase = factor(cells_total$Phase, levels = c('M', 'I', 'LI', 'NF', 'F'))

# ---- GLM for p-values ----
out_df = data.frame()
for (i in 0:5) {
  cells   = subset(cells_total, sub.cluster == paste0('1_', i) & Status != 'NRM')
  mat     = cbind(cells$n_cells_in, cells$n_cells - cells$n_cells_in)
  mod     = glm(mat ~ Status, family = binomial(), data = cells)
  av_     = car::Anova(mod, 3)
  pair    = pairs(emmeans(mod, 'Status'), adjust = 'none') %>% as.data.frame()

  newd = data.frame(
    cluster   = i,
    av_p      = av_$`Pr(>Chisq)`[1],
    d_m_p     = pair$p.value[pair$contrast == 'M - D'],
    f_m_p     = pair$p.value[pair$contrast == 'M - F'],
    d_f_p     = pair$p.value[pair$contrast == 'D - F']
  )
  out_df = rbind(out_df, newd)
}
out_df$signif = ifelse(out_df$av_p < 0.05, '*', NA)

# ---- plotter ----
prop_plotter_sub1 = function(cells_total, subcluster_id, d_m_p, f_m_p, d_f_p,
                              y_label = '% of 1_RGC') {

  dat = subset(cells_total, sub.cluster == subcluster_id & Status != 'NRM')

  summary_dat = dat %>%
    group_by(Phase) %>%
    summarise(
      mean_prop = mean(n_cells_in / n_cells, na.rm = TRUE),
      se        = sd(n_cells_in / n_cells,   na.rm = TRUE) / sqrt(n()),
      .groups   = 'drop'
    )

  d_m_fmt = p_value_annotation(d_m_p)
  f_m_fmt = p_value_annotation(f_m_p)
  d_f_fmt = p_value_annotation(d_f_p)

  y_max      = max(dat$n_cells_in / dat$n_cells, na.rm = TRUE)
  y_min      = 0
  y_lower    = y_max * 1.10
  y_upper    = y_max * 1.35
  y_axis_max = y_max * 1.50

  ggplot(dat, aes(x = Phase, y = n_cells_in / n_cells)) +
    geom_bar(data        = summary_dat,
             aes(x = Phase, y = mean_prop, fill = Phase),
             stat        = 'identity',
             width       = 0.6,
             inherit.aes = FALSE) +
    geom_errorbar(data        = summary_dat,
                  aes(x = Phase, ymin = mean_prop - se, ymax = mean_prop + se),
                  width       = 0.25,
                  linewidth   = 0.7,
                  inherit.aes = FALSE) +
    geom_point(size     = 1,
               position = position_jitter(width = 0.1, seed = 42)) +
    theme_classic() +
    scale_x_discrete(drop = FALSE) +
    scale_y_continuous(labels = scales::percent,
                       limits = c(y_min, y_axis_max),
                       expand = c(0, 0)) +
    scale_fill_manual(values = colors) +
    labs(x     = 'Phase',
         y     = y_label,
         title = paste0('Cluster ', subcluster_id)) +
    theme(legend.position = 'none') +
    geom_signif(xmin       = 1, xmax = 1.9,
                y_position = y_lower,
                annotation = d_m_fmt$label,
                color      = 'black', tip_length = c(0, 0),
                textsize   = d_m_fmt$size, vjust = d_m_fmt$adjust) +
    geom_signif(xmin       = 2.1, xmax = 5,
                y_position = y_lower,
                annotation = d_f_fmt$label,
                color      = 'black', tip_length = c(0, 0),
                textsize   = d_f_fmt$size, vjust = d_f_fmt$adjust) +
    geom_signif(xmin       = 1, xmax = 5,
                y_position = y_upper,
                annotation = f_m_fmt$label,
                color      = 'black', tip_length = c(0, 0),
                textsize   = f_m_fmt$size, vjust = f_m_fmt$adjust)
}

# ---- plots ----
prop_1_0 = prop_plotter_sub1(cells_total, '1_0',
              d_m_p = out_df$d_m_p[out_df$cluster == 0],
              f_m_p = out_df$f_m_p[out_df$cluster == 0],
              d_f_p = out_df$d_f_p[out_df$cluster == 0])

prop_1_1 = prop_plotter_sub1(cells_total, '1_1',
              d_m_p = out_df$d_m_p[out_df$cluster == 1],
              f_m_p = out_df$f_m_p[out_df$cluster == 1],
              d_f_p = out_df$d_f_p[out_df$cluster == 1])

prop_1_2 = prop_plotter_sub1(cells_total, '1_2',
              d_m_p = out_df$d_m_p[out_df$cluster == 2],
              f_m_p = out_df$f_m_p[out_df$cluster == 2],
              d_f_p = out_df$d_f_p[out_df$cluster == 2])

prop_1_3 = prop_plotter_sub1(cells_total, '1_3',
              d_m_p = out_df$d_m_p[out_df$cluster == 3],
              f_m_p = out_df$f_m_p[out_df$cluster == 3],
              d_f_p = out_df$d_f_p[out_df$cluster == 3])

# ---- save ----
for (i in 0:3) {
  p = get(paste0('prop_1_', i))
  ggsave(plot   = p,
         file   = paste0('rgc_prop_1_', i, '.svg'),
         device = 'svg',
         units  = 'in',
         width  = 2,
         height = 2,
         path   = 'Manuscript/Plots/Manuscript v1.2/RGC_sub_props')
}
