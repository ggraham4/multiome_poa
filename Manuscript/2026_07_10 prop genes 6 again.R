library(Seurat)
library(tidyverse)
library(ggplot2)
library(emmeans)
library(ggsignif)

obj  = readRDS("C:/Users/Gabe/Desktop/nemo.orig_harmony.integration_all_testd_clusters.rds")
sub_6 = subset(obj,res0.8_50nn_40PC_45LSI==6)



plot_proportions_sub_6 = function(sub_6, gene) {
  
  colors = c('#1965B0', '#4EB265', '#F7F056', '#7BAFDE', '#DC050C')

make_summary = function(dat, numerator_col, denominator_col) {
  dat %>%
    group_by(Phase) %>%
    summarise(
      mean_prop = mean(.data[[numerator_col]] / .data[[denominator_col]], na.rm = TRUE),
      se = sd(.data[[numerator_col]] / .data[[denominator_col]], na.rm = TRUE) / sqrt(n()),
      .groups = 'drop'
    )
}

format_pval = function(p) {
  if (p >= 0.05) {
    list(label = paste0('p = ', round(p, 3)), size = 2, adjust = 0)
  } else if (p < 0.001) {
    list(label = '***', size = 7, adjust = .6)
  } else if (p < 0.01) {
    list(label = '**', size = 7, adjust = .6)
  } else {
    list(label = '*', size = 7,adjust = .6)
  }
}
  
  sub_6$gene = ifelse(sub_6@assays$RNA$data[gene,] > 0, T, F)
  
  gex = sub_6@meta.data %>%
    group_by(individual, Status) %>%
    summarize(n_gene = sum(gene == T))
  
  tot = sub_6@meta.data %>%
    group_by(individual, Status) %>%
    summarize(total_cells = n())
  
  tog = gex %>%
    right_join(tot, by = c('individual'))
  
  tog_6_gene = subset(tog, Status.x != c('NRM'))
  
  ### Stats ###
  gene_glm = glm(cbind(tog_6_gene$n_gene, tog_6_gene$total_cells - tog_6_gene$n_gene) ~ Status.x,
                 data = tog_6_gene,
                 family = binomial('logit'))
  
  gene_pairs = as.data.frame(pairs(emmeans(gene_glm, 'Status.x'), adjust = 'none'))
  
  get_p = function(g1, g2) {
    row = gene_pairs %>% filter(
      (contrast == paste(g1, '-', g2)) | (contrast == paste(g2, '-', g1))
    )
    if (nrow(row) == 0) return(NA)
    row$p.value[1]
  }
  
  p_M_I = get_p('M', 'D')
  p_M_F = get_p('M', 'F')
  p_I_F = get_p('D', 'F')
  
  tog_6_gene$Phase = as.character(tog_6_gene$Status.x)
  tog_6_gene$Phase = ifelse(tog_6_gene$Phase == 'D', 'I', tog_6_gene$Phase)
  tog_6_gene$Phase = ifelse(tog_6_gene$Phase == 'E', 'LI', tog_6_gene$Phase)
  tog_6_gene$Phase = factor(tog_6_gene$Phase, levels = c('M','I','LI','NF','F'))
  
  summary_gene = make_summary(tog_6_gene, 'n_gene', 'total_cells')
  
  # Two tiers: lower for MI and IF, upper for MF
  #y_min      = min((tog_6_gene$n_gene / tog_6_gene$total_cells), na.rm = TRUE) * 0.95
  y_max = max((tog_6_gene$n_gene / tog_6_gene$total_cells), na.rm = TRUE)
  y_lower = y_max * 1.10
  y_upper = y_max * 1.25
  y_axis_max = y_max * 1.40
  
  
  p=gene_prop = ggplot(tog_6_gene, aes(x = Phase, y = n_gene / total_cells)) +
    geom_bar(data = summary_gene,
             aes(x = Phase, y = mean_prop, fill = Phase),
             stat = 'identity', inherit.aes = FALSE) +
    geom_errorbar(data = summary_gene,
                  aes(x = Phase, ymin = mean_prop - se, ymax = mean_prop + se),
                  width = 0.4, inherit.aes = FALSE) +
    geom_point(size = 1) +
    theme_classic() +
    scale_x_discrete(drop = FALSE) +
    scale_y_continuous(labels = scales::percent,
                       limits = c(0, y_axis_max)) +
    labs(x = 'Phase', y = paste0('Proportion ', gene, '+')) +
    theme(legend.position = 'none') +
    scale_fill_manual(values = colors) +
    # Lower tier: M vs I
    geom_signif(xmin = 1, xmax = 1.9,
                y_position = y_lower,
                annotation = format_pval(p_M_I)$label,
                color = "black", 
                tip_length = c(0.00, 0.00), 
                textsize = format_pval(p_M_I)$size,
                 vjust = format_pval(p_M_I)$adjust) +
    # Lower tier: I vs F
    geom_signif(xmin = 2.1, xmax = 5,
                y_position = y_lower,
                annotation = format_pval(p_I_F)$label,
                color = "black", 
                tip_length = c(0.0, 0.0), 
                textsize = format_pval(p_I_F)$size,
                 vjust = format_pval(p_I_F)$adjust) +
    # Upper tier: M vs F
    geom_signif(xmin = 1, xmax = 5,
                y_position = y_upper,
                annotation = format_pval(p_M_F)$label,
                color = "black",
                tip_length = c(0.0, 0.0),
                textsize = format_pval(p_M_F)$size,
                vjust = format_pval(p_M_F)$adjust)
 # return(anova(gene_glm, test= 'Chisq'))
  return(p)
}
 nmbr= plot_proportions_sub_6(sub_6, 'nmbr')

#ggsave(plot = nmbr,
#       file = paste0('prop_nmbr_6.svg'),
#       device = "svg",
#       units = "in",
#       width =2,
#       height = 2,
#       path = paste0("Manuscript/Plots/Manuscript v1.2.1/6_props"))

 pgr= plot_proportions_sub_6(sub_6, 'pgr')

#ggsave(plot = pgr,
#       file = paste0('prop_pgr_6.svg'),
#       device = "svg",
#       units = "in",
#       width =2,
#       height = 2,
#       path = paste0("Manuscript/Plots/Manuscript v1.2.1/6_props"))

 pgrlike= plot_proportions_sub_6(sub_6, 'LOC111568069')

#ggsave(plot = pgrlike,
#       file = paste0('prop_pgrlike_6.svg'),
#       device = "svg",
#       units = "in",
#       width =2,
#       height = 2,
#       path = paste0("Manuscript/Plots/Manuscript v1.2.1/6_props"))

npy7r= plot_proportions_sub_6(sub_6, 'npy7r')

#ggsave(plot = npy7r,
#       file = paste0('prop_npy7r_6.svg'),
#       device = "svg",
#       units = "in",
#       width =2,
#       height = 2,
#       path = paste0("Manuscript/Plots/Manuscript v1.2.1/6_props"))

cckb= plot_proportions_sub_6(sub_6, 'cckb')

ggsave(plot = cckb,
       file = paste0('prop_cckb_6.svg'),
       device = "svg",
       units = "in",
       width =2,
       height = 2,
       path = paste0("Manuscript/Plots/Manuscript v1.3/6_props"))

tacr3a = plot_proportions_sub_6(sub_6, 'tacr3a')

ggsave(plot = tacr3a,
       file = paste0('prop_tacr3a_6.svg'),
       device = "svg",
       units = "in",
       width =2,
       height = 2,
       path = paste0("Manuscript/Plots/Manuscript v1.3/6_props"))

drd3 = plot_proportions_sub_6(sub_6, 'drd3')

ggsave(plot = drd3,
       file = paste0('prop_drd3_6.svg'),
       device = "svg",
       units = "in",
       width =2,
       height = 2,
       path = paste0("Manuscript/Plots/Manuscript v1.3/6_props"))



