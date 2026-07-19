#libs
library(ggplot2)
library(ggsignif)
library(tidyverse)
library(emmeans)
library(Seurat)
library(CytoTRACE)

# read in data, subset to 6_NPO
obj  = readRDS("~/Desktop/optimal_clustering_rna_only.rds")
sub_6 = subset(obj,res0.8_50nn_40PC_45LSI==6)
degs = read.csv('/Users/ggraham/Desktop/multiome_poa/DEG Outputs/FINAL degs classified w singular.csv')
degs_6 = subset(degs, cluster == 6)


# define genes to examine
genes_interest = c(
"LOC111568069"
)

sub_6$cyto = CytoTRACE(sub_6@assays$RNA$data%>%as.matrix())$CytoTRACE
# loop throgh genes interest and test for cyto differences
cdp =  function(sub_6, gene){
  
  if(!require(lme4)){library(lme4)}
  if(!require(ggplot2)){library(ggplot2)}
  if(!require(car)){library(car)}
  if(!require(ggsignif)){library(ggsignif)}
  if(!require(tidyverse)){library(tidyverse)}
  if(!require(emmeans)){library(emmeans)}
  
  colors = c('#1965B0', '#4EB265', '#F7F056', '#DC050C')
  
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
  
  if(!'cyto' %in% colnames(sub_6@meta.data)){
    stop('"cyto" column not in sub_6@meta.data')
  }
  
  temp = sub_6@meta.data
  sub_6_gex = sub_6@assays$RNA$data[gene,]
  gene_positive_cells = sub_6_gex > 0
  temp$gene_pos = gene_positive_cells
  temp = subset(temp, Status != 'NRM')
  
  # --- key change: restrict to gene+ cells only ---
  temp = subset(temp, gene_pos == TRUE)
  
  if(nrow(temp) == 0){
    stop('No gene-positive cells remain after subsetting (Status != "NRM" & gene_pos == TRUE)')
  }
  
  # model: Status only, still regressing out nCount_RNA, random intercept for individual
  model = lmer(cyto ~ nCount_RNA + Status + (1|individual),
               data = temp)
  
  av_ = car::Anova(model, 3)
  type_3_p.value_phase = av_$`Pr(>Chisq)`[rownames(av_) == 'Status']
  if(type_3_p.value_phase > 0.05) message('Type III test is not significant for Phase')
  
  pair_ = pairs(emmeans(model, 'Status'), adjust = 'none') %>% as.data.frame()
  
  get_p = function(g1, g2) {
    row = pair_ %>% filter(
      (contrast == paste(g1, '-', g2)) | (contrast == paste(g2, '-', g1))
    )
    if (nrow(row) == 0) return(NA)
    row$p.value[1]
  }
  
  ind_summary = temp %>%
    group_by(individual, Status) %>%
    summarize(mean_cyto = mean(cyto, na.rm = TRUE), .groups = 'drop')
  
  ind_summary$Phase = as.character(ind_summary$Status)
  ind_summary$Phase = ifelse(ind_summary$Phase == 'D',  'I',  ind_summary$Phase)
  ind_summary$Phase = ifelse(ind_summary$Phase == 'E',  'LI', ind_summary$Phase)
  ind_summary$Phase = factor(ind_summary$Phase, levels = c('M','I','LI','NF','F'))
  
  y_min_global = min(ind_summary$mean_cyto, na.rm = TRUE)*0.95
  y_max_global = max(ind_summary$mean_cyto, na.rm = TRUE)
  y_lower    = y_max_global * 1.10
  y_upper    = y_max_global * 1.35
  y_axis_max = y_max_global * 1.50
  
  dat = ind_summary
  dat_box = subset(dat, Phase != 'NF')
  
  p_mi = get_p('M', 'D')
  p_mf = get_p('M', 'F')
  p_if = get_p('D', 'F')
  
  p = ggplot(dat, aes(x = Phase, y = mean_cyto)) +
    geom_boxplot(data = dat_box,
                 aes(x = Phase, y = mean_cyto, fill = Phase),
                 outlier.shape = NA, width = 0.6) +
    geom_point(size = 1,
               position = position_jitter(width = 0.1, seed = 42)) +
    theme_classic() +
    scale_x_discrete(drop = FALSE) +
    scale_y_continuous(limits = c(y_min_global, y_axis_max)) +
    scale_fill_manual(values = colors) +
    labs(x = 'Phase',
         y = 'Mean CytoTRACE') +
    theme(legend.position = 'none') +
    geom_signif(xmin = 1, xmax = 1.9,
                y_position = y_lower,
                annotation = format_pval(p_mi)$label,
                color = "black", tip_length = c(0, 0),
                textsize = format_pval(p_mi)$size,
                vjust = format_pval(p_mi)$adjust) +
    geom_signif(xmin = 2.1, xmax = 5,
                y_position = y_lower,
                annotation = format_pval(p_if)$label,
                color = "black", tip_length = c(0, 0),
                textsize = format_pval(p_if)$size,
                vjust = format_pval(p_if)$adjust) +
    geom_signif(xmin = 1, xmax = 5,
                y_position = y_upper,
                annotation = format_pval(p_mf)$label,
                color = "black", tip_length = c(0, 0),
                textsize = format_pval(p_mf)$size,
                vjust = format_pval(p_mf)$adjust)
  
  p
}

for(i in genes_interest){
  c = cdp(sub_6, i)
  
  ggsave(plot = c,
         file = paste0('arlike', '_cyto_6.svg'),
         device = "svg",
         units = "in",
         width = 2,
         height = 2,
         path = paste0("Manuscript/Plots/Manuscript v1.3/cyto_genes_6"))
  
}

