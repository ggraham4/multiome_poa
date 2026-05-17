library(Seurat)
library(tidyverse)
library(ggplot2)
library(ggsignif)
library(lme4)
library(car)
library(emmeans)
library(CytoTRACE)

# ---- data prep ----
obj   = readRDS("~/Desktop/optimal_clustering_rna_only.rds")
sub_6 = FindSubCluster(obj, 6, graph.name = 'harmony.wsnn')
Idents(sub_6) <- 'sub.cluster'
sub_6 = subset(sub_6, final_clusters == 6)
sub_6$Status = factor(sub_6$Status, levels = c('NRM', 'M', 'D', 'E', 'NF', 'F'))

sub_6@meta.data$Phase = as.character(sub_6@meta.data$Status)
sub_6@meta.data$Phase = ifelse(sub_6@meta.data$Phase == 'D',  'I',  sub_6@meta.data$Phase)
sub_6@meta.data$Phase = ifelse(sub_6@meta.data$Phase == 'E',  'LI', sub_6@meta.data$Phase)
sub_6@meta.data$Phase = factor(sub_6@meta.data$Phase, levels = c('NRM', 'M', 'I', 'LI', 'NF', 'F'))

cyto      = CytoTRACE(sub_6@assays$RNA$data %>% as.matrix())
sub_6$cyto = cyto$CytoTRACE

degs_plasticity = c(
  'LOC111588913',
  'cntn4',
  'LOC111567620',
  'pcdh10b',
  'sdc2',
  'LOC111585095',
  'bcan'
)

plasticity_positive = colSums(sub_6@assays$RNA$data[degs_plasticity, ]) > 0
sub_6$gene_pos = plasticity_positive

# ---- helpers ----
colors = c('#1965B0', '#4EB265', '#F7F056', '#DC050C')

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

# ---- model ----
temp = sub_6@meta.data
temp = subset(temp, Status != 'NRM')

model  = lmer(cyto ~ nCount_RNA + (Status * gene_pos) + (1 | individual), data = temp)
av_    = car::Anova(model, 3)

type3_p_status   = av_$`Pr(>Chisq)`[rownames(av_) == 'Status']
type3_p_gene_pos = av_$`Pr(>Chisq)`[rownames(av_) == 'gene_pos']
if (type3_p_status   > 0.05) message('Type III not significant for Status')
if (type3_p_gene_pos > 0.05) message('Type III not significant for gene_pos')

pair_ = pairs(emmeans(model, c('Status'), by = 'gene_pos'), adjust = 'none') %>%
  as.data.frame()

genepos_label = p_value_annotation(type3_p_gene_pos)$label

get_p = function(g1, g2, gpos) {
  row = pair_ %>% filter(
    gene_pos == gpos,
    (contrast == paste(g1, '-', g2)) | (contrast == paste(g2, '-', g1))
  )
  if (nrow(row) == 0) return(NA)
  row$p.value[1]
}

# ---- per-individual summary ----
ind_summary = temp %>%
  group_by(individual, Status, gene_pos) %>%
  summarize(mean_cyto = mean(cyto, na.rm = TRUE), .groups = 'drop')

ind_summary$Phase = as.character(ind_summary$Status)
ind_summary$Phase = ifelse(ind_summary$Phase == 'D',  'I',  ind_summary$Phase)
ind_summary$Phase = ifelse(ind_summary$Phase == 'E',  'LI', ind_summary$Phase)
ind_summary$Phase = factor(ind_summary$Phase, levels = c('M', 'I', 'LI', 'NF', 'F'))

# ---- shared y axis ----
y_min_global = min(ind_summary$mean_cyto, na.rm = TRUE) * 0.95
y_max_global = max(ind_summary$mean_cyto, na.rm = TRUE)
y_lower      = y_max_global * 1.10
y_upper      = y_max_global * 1.35
y_axis_max   = y_max_global * 1.50

# ---- panel builder ----
make_panel = function(gpos, title_suffix) {
  dat     = subset(ind_summary, gene_pos == gpos)
  dat_box = subset(dat, Phase != 'NF')

  p_mi = get_p('M', 'D', gpos)
  p_mf = get_p('M', 'F', gpos)
  p_if = get_p('D', 'F', gpos)

  ggplot(dat, aes(x = Phase, y = mean_cyto)) +
    geom_boxplot(data = dat_box,
                 aes(x = Phase, y = mean_cyto, fill = Phase),
                 outlier.shape = NA, width = 0.6) +
    geom_point(size = 1,
               position = position_jitter(width = 0.1, seed = 42)) +
    theme_classic() +
    scale_x_discrete(drop = FALSE) +
    scale_y_continuous(limits = c(y_min_global, y_axis_max)) +
    scale_fill_manual(values = colors) +
    labs(x        = 'Phase',
         y        = 'Mean CytoTRACE',
         title    = paste0('Plasticity module ', title_suffix),
         subtitle = paste0('pos vs neg: ', genepos_label)) +
    theme(legend.position = 'none') +
    geom_signif(xmin = 1, xmax = 1.9,
                y_position = y_lower,
                annotation = p_value_annotation(p_mi)$label,
                color = 'black', tip_length = c(0, 0),
                textsize = p_value_annotation(p_mi)$size,
                vjust    = p_value_annotation(p_mi)$adjust) +
    geom_signif(xmin = 2.1, xmax = 5,
                y_position = y_lower,
                annotation = p_value_annotation(p_if)$label,
                color = 'black', tip_length = c(0, 0),
                textsize = p_value_annotation(p_if)$size,
                vjust    = p_value_annotation(p_if)$adjust) +
    geom_signif(xmin = 1, xmax = 5,
                y_position = y_upper,
                annotation = p_value_annotation(p_mf)$label,
                color = 'black', tip_length = c(0, 0),
                textsize = p_value_annotation(p_mf)$size,
                vjust    = p_value_annotation(p_mf)$adjust)
}

# ---- combine ----
p_pos = make_panel(TRUE,  'positive')
p_neg = make_panel(FALSE, 'negative') + labs(subtitle = element_blank())

cyto_ecm = (p_pos | p_neg) + plot_layout(axes = 'collect')
cyto_ecm


ggsave(plot = cyto_ecm,
       file = paste0('cyto_ecm_presence.svg'),
       device = "svg",
       units = "in",
       width =4,
       height = 2.5,
       path = paste0("Manuscript/Plots/Manuscript v1.2/"))
