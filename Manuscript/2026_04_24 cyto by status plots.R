library(CytoTRACE)
library(Seurat)
library(tidyverse)
library(lme4)
library(car)
refactor =readRDS('Functions/refactor_status.rds')
cyto_by_status = function(obj_sub, cluster_name) {
  
  colors = c('#1965B0', '#4EB265', '#F7F056', '#7BAFDE', '#DC050C')
  
  format_pval = function(p) {
    if (p >= 0.05) {
      list(label = paste0('p = ', round(p, 3)), size = 2, adjust = 0)
    } else if (p < 0.001) {
      list(label = '***', size = 7, adjust = .6)
    } else if (p < 0.01) {
      list(label = '**', size = 7, adjust = .6)
    } else {
      list(label = '*', size = 7, adjust = .6)
    }
  }
  
  mod = lmer(cyto ~ Status + (1|individual), 
             data = subset(obj_sub@meta.data, Status != 'NRM'))
  
  av_ = car::Anova(mod, 3)
  p_status = av_$`Pr(>Chisq)`[rownames(av_) == 'Status']
  status_label = format_pval(p_status)$label
  
  pair_ = pairs(emmeans(mod, 'Status'), adjust = 'none') %>% as.data.frame()
  
  get_p = function(g1, g2) {
    row = pair_ %>% filter(
      (contrast == paste(g1, '-', g2)) | (contrast == paste(g2, '-', g1))
    )
    if (nrow(row) == 0) return(NA)
    row$p.value[1]
  }
  
  p_M_I = get_p('M', 'D')
  p_M_F = get_p('M', 'F')
  p_I_F = get_p('D', 'F')
  
  ind_summary = obj_sub@meta.data %>%
    group_by(individual, Status) %>%
    summarize(mean_cyto = mean(cyto), .groups = 'drop') %>%
    refactor()
  
  ind_summary$Phase = as.character(ind_summary$Status)
  ind_summary$Phase = ifelse(ind_summary$Phase == 'D',  'I',  ind_summary$Phase)
  ind_summary$Phase = ifelse(ind_summary$Phase == 'E',  'LI', ind_summary$Phase)
  ind_summary$Phase = factor(ind_summary$Phase, levels = c('M','I','LI','NF','F'))
  
  ind_summary = subset(ind_summary, Status != 'NRM')
  
  dat_box = subset(ind_summary, Phase != 'NF')
  
  y_min      = min(ind_summary$mean_cyto, na.rm = TRUE) * 0.95
  y_max      = max(ind_summary$mean_cyto, na.rm = TRUE)
  y_lower    = y_max * 1.10
  y_upper    = y_max * 1.25
  y_axis_max = y_max * 1.40
  
  ggplot(ind_summary, aes(x = Phase, y = mean_cyto)) +
    geom_boxplot(data = dat_box,
                 aes(x = Phase, y = mean_cyto, fill = Phase),
                 outlier.shape = NA, width = 0.6) +
    geom_point(size = 1, position = position_jitter(width = 0.1, seed = 42)) +
    theme_classic() +
    scale_x_discrete(drop = FALSE) +
    scale_y_continuous(limits = c(y_min, y_axis_max)) +
    scale_fill_manual(values = colors) +
    labs(x = 'Phase', y = 'Mean CytoTRACE',
         title = paste0('Cluster ', cluster_name, ' | Status: ', status_label)) +
    theme(legend.position = 'none') +
    geom_signif(xmin = 1, xmax = 1.9,
                y_position = y_lower,
                annotation = format_pval(p_M_I)$label,
                color = "black", tip_length = c(0, 0),
                textsize = format_pval(p_M_I)$size,
                vjust = format_pval(p_M_I)$adjust) +
    geom_signif(xmin = 2.1, xmax = 5,
                y_position = y_lower,
                annotation = format_pval(p_I_F)$label,
                color = "black", tip_length = c(0, 0),
                textsize = format_pval(p_I_F)$size,
                vjust = format_pval(p_I_F)$adjust) +
    geom_signif(xmin = 1, xmax = 5,
                y_position = y_upper,
                annotation = format_pval(p_M_F)$label,
                color = "black", tip_length = c(0, 0),
                textsize = format_pval(p_M_F)$size,
                vjust = format_pval(p_M_F)$adjust)
}

#saveRDS(cyto_by_status, 'Functions/cyto_by_status.rds')


# read in data
obj  = readRDS("~/Desktop/optimal_clustering_rna_only.rds")
sub_6 = subset(obj,res0.8_50nn_40PC_45LSI==6)
sub_6$cyto = CytoTRACE(sub_6@assays$RNA$data%>%as.matrix())$CytoTRACE

sub_0 = subset(obj,res0.8_50nn_40PC_45LSI==0)
sub_0$cyto = CytoTRACE(sub_0@assays$RNA$data%>%as.matrix())$CytoTRACE

sub_9 = subset(obj,res0.8_50nn_40PC_45LSI==9)
sub_9$cyto = CytoTRACE(sub_9@assays$RNA$data%>%as.matrix())$CytoTRACE

sub_23 = subset(obj,res0.8_50nn_40PC_45LSI==23)
sub_23$cyto = CytoTRACE(sub_23@assays$RNA$data%>%as.matrix())$CytoTRACE

cyto_by_status_23 = cyto_by_status(sub_23, '23')
# no diff to be expected 

cyto_by_status_6 = cyto_by_status(sub_6, '6')

cyto_by_status_0 =cyto_by_status(sub_0, '0')

cyto_by_status_9=cyto_by_status(sub_9, '9')

ggsave(plot = cyto_by_status_6,
       file = 'cyto_by_status_6.svg',
       device = "svg",
       units = "in",
       width = 2,
       height = 2,
       path = paste0("Manuscript/Plots/Manuscript v1.2/cyto_by_status"))


ggsave(plot = cyto_by_status_0,
       file = 'cyto_by_status_0.svg',
       device = "svg",
       units = "in",
       width = 2,
       height = 2,
       path = paste0("Manuscript/Plots/Manuscript v1.2/cyto_by_status"))


ggsave(plot = cyto_by_status_9,
       file = 'cyto_by_status_9.svg',
       device = "svg",
       units = "in",
       width = 2,
       height = 2,
       path = paste0("Manuscript/Plots/Manuscript v1.2/cyto_by_status"))



