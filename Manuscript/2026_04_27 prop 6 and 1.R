library(Seurat)
library(patchwork)
library(tidyverse)
library(ggplot2)
library(Polychrome)
library(emmeans)
library(ggsignif)
  clown_go = readRDS("Functions/clown_go2")  
library(clusterProfiler)




obj  = readRDS("~/Desktop/optimal_clustering_rna_only.rds")
colors = c('#1965B0', '#4EB265', '#F7F056', '#7BAFDE', '#DC050C')

status_to_phase = c(
  "D"='I',
  'M' = 'M',
  'F' ='F',
  'NF' ='NF',
  'NRM' = 'NRM',
  'E' = 'LI'
)

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

proportion_plotter = function(data_frame, cluster,
                              cluster_col = 'res0.8_50nn_40PC_45LSI',
                              d_m_p = NULL, f_m_p = NULL, d_f_p = NULL) {
  
  colors = c('#1965B0', '#4EB265', '#F7F056', '#7BAFDE', '#DC050C')
  
  dag = subset(data_frame, data_frame[[cluster_col]] == cluster & Status != 'NRM')
  dag$Phase = status_to_phase[dag$Status]
  dag$Phase = factor(dag$Phase, levels = c('M', 'I', 'LI', 'NF', 'F'))
  dag$prop  = dag$ncells / dag$total_cells
  
  dag_summary = dag %>%
    group_by(Phase) %>%
    summarise(
      mean_prop = mean(prop, na.rm = TRUE),
      se        = sd(prop,   na.rm = TRUE) / sqrt(n()),
      .groups   = 'drop'
    )
  
  # --- y scaling (same multiplier logic as expression/cyto functions) ---
  y_max      = max(dag$prop, na.rm = TRUE)
  y_min      = 0
  y_lower    = y_max * 1.10
  y_upper    = y_max * 1.35
  y_axis_max = y_max * 1.50
  
  # --- format p-values if supplied ---
  d_m_fmt = if (!is.null(d_m_p)) p_value_annotation(d_m_p) else list(label = '', size = 3, adjust = 0)
  f_m_fmt = if (!is.null(f_m_p)) p_value_annotation(f_m_p) else list(label = '', size = 3, adjust = 0)
  d_f_fmt = if (!is.null(d_f_p)) p_value_annotation(d_f_p) else list(label = '', size = 3, adjust = 0)
  
  p = ggplot(dag, aes(x = Phase, y = prop)) +
    geom_bar(data    = dag_summary,
             aes(x = Phase, y = mean_prop, fill = Phase),
             stat    = 'identity',
             width   = 0.6,
             inherit.aes = FALSE) +
    geom_errorbar(data = dag_summary,
                  aes(x = Phase, ymin = mean_prop - se, ymax = mean_prop + se),
                  width     = 0.4,
                  linewidth = 0.5,
                  inherit.aes = FALSE) +
    geom_point(size     = 1,
               position = position_jitter(width = 0.1, seed = 42)) +
    scale_x_discrete(drop = FALSE) +
    scale_y_continuous(labels    = scales::percent,
                       limits    = c(y_min, y_axis_max),
                       expand    = c(0, 0)) +
    scale_fill_manual(values = colors) +
    theme_classic() +
    labs(x     = 'Phase',
         y     = '% of Cells',
         title = paste0('Cluster ', cluster)) +
    theme(legend.position = 'none')
  
  # --- add brackets only when p-values are supplied ---
  if (!is.null(d_m_p)) {
    p = p + geom_signif(xmin = 1, xmax = 1.9,
                        y_position = y_lower,
                        annotation = d_m_fmt$label,
                        color      = 'black', tip_length = c(0, 0),
                        textsize   = d_m_fmt$size, vjust = d_m_fmt$adjust)
  }
  if (!is.null(d_f_p)) {
    p = p + geom_signif(xmin = 2.1, xmax = 5,
                        y_position = y_lower,
                        annotation = d_f_fmt$label,
                        color      = 'black', tip_length = c(0, 0),
                        textsize   = d_f_fmt$size, vjust = d_f_fmt$adjust)
  }
  if (!is.null(f_m_p)) {
    p = p + geom_signif(xmin = 1, xmax = 5,
                        y_position = y_upper,
                        annotation = f_m_fmt$label,
                        color      = 'black', tip_length = c(0, 0),
                        textsize   = f_m_fmt$size, vjust = f_m_fmt$adjust)
  }
  
  return(p)
}

neurons_only <- subset(obj, 
                     #oligos
                     res0.8_50nn_40PC_45LSI!=2&
                     #microglia
                    res0.8_50nn_40PC_45LSI!=11&
                    #opcs
                    res0.8_50nn_40PC_45LSI!=13&
                    #dividing glia
                    res0.8_50nn_40PC_45LSI!=26&
                    #leuko
                    res0.8_50nn_40PC_45LSI!=20&
                    #ependymal
                    res0.8_50nn_40PC_45LSI!=15
                    &  res0.8_50nn_40PC_45LSI!=1
                    )

total_cells_neuron = neurons_only@meta.data%>%
  group_by(individual, Status)%>%
  summarize(total_cells = n())

n_cells_neuron=neurons_only@meta.data%>%
  group_by(individual, res0.8_50nn_40PC_45LSI)%>%
  summarize(ncells = n())

joint_neuron = total_cells_neuron%>%
  left_join(n_cells_neuron, by = 'individual')

joint_neuron$Status = as.character(joint_neuron$Status)

dat_neuron = data.frame()
for(i in unique(joint_neuron$res0.8_50nn_40PC_45LSI)){
  sub= subset(joint_neuron, res0.8_50nn_40PC_45LSI ==i & Status %in% c('M','D','F'))
  mat = cbind(sub$ncells, sub$total_cells-sub$ncells)
  
  mod = glm(mat ~ Status, data = sub, family = 'binomial')
  p = anova(mod, 'III')
  p.value = p$`Pr(>Chi)`[2]
  
  pair = pairs(emmeans(mod, 'Status'), adjust = 'none')%>%as.data.frame()
  
  newd = data.frame(cluster = i,
                    p = p.value,
                    d_m_p.value = pair$p.value[pair$contrast == 'D - M'],
                    f_m_p.value = pair$p.value[pair$contrast == 'F - M'],
                    d_f_p.value = pair$p.value[pair$contrast == 'D - F'],
                    d_m_estimate = pair$estimate[pair$contrast == 'D - M'],
                    f_m_estimate = pair$estimate[pair$contrast == 'F - M'],
                    d_f_estimate = pair$estimate[pair$contrast == 'D - F']
)
  dat_neuron = rbind(newd, dat_neuron)

  
}

dat_neuron$av_q.value = p.adjust(dat_neuron$p, 'fdr',nrow(dat_neuron))
dat_neuron$sig = ifelse(dat_neuron$av_q.value < 0.05, '*', NA)

total_cells = obj@meta.data%>%
  group_by(individual, Status)%>%
  summarize(total_cells = n())

n_cells=obj@meta.data%>%
  group_by(individual, res0.8_50nn_40PC_45LSI)%>%
  summarize(ncells = n())

joint = total_cells%>%
  left_join(n_cells, by = 'individual')

joint$Status = as.character(joint$Status)

dat = data.frame()
for(i in 0:26){
  sub= subset(joint, res0.8_50nn_40PC_45LSI ==i & Status %in% c('M','D','F'))
  mat = cbind(sub$ncells, sub$total_cells-sub$ncells)
  
  mod = glm(mat ~ Status, data = sub, family = 'binomial')
  p = anova(mod, 'III')
  p.value = p$`Pr(>Chi)`[2]
  
  pair = pairs(emmeans(mod, 'Status'), adjust = 'none')%>%as.data.frame()
  
  newd = data.frame(cluster = i,
                    p = p.value,
                    d_m_p.value = pair$p.value[pair$contrast == 'D - M'],
                    f_m_p.value = pair$p.value[pair$contrast == 'F - M'],
                    d_f_p.value = pair$p.value[pair$contrast == 'D - F'],
                    d_m_estimate = pair$estimate[pair$contrast == 'D - M'],
                    f_m_estimate = pair$estimate[pair$contrast == 'F - M'],
                    d_f_estimate = pair$estimate[pair$contrast == 'D - F']
)
  dat = rbind(newd, dat)

  
}

dat$av_q.value = p.adjust(dat$p, 'fdr', 27)
dat$sig = ifelse(dat$av_q.value < 0.05, '*', NA)
dat$singular = F

clust_1 = proportion_plotter(joint, 1,
            d_m_p = dat$d_m_p.value[dat$cluster == 1],
            f_m_p = dat$f_m_p.value[dat$cluster == 1],
            d_f_p = dat$d_f_p.value[dat$cluster == 1]) +
  labs(y = '% of Cells')

clust_1

ggsave(plot = clust_1,
       file = "clust_1_prop.svg",
       device = "svg",
       units = "in",
       width = 2,
       height = 2,
       path = '/Users/ggraham/Desktop/multiome_poa/Manuscript/Plots/Manuscript v1.2/prop_cells')



clust_6 = proportion_plotter(joint_neuron, 6,
            d_m_p = dat_neuron$d_m_p.value[dat_neuron$cluster == 6],
            f_m_p = dat_neuron$f_m_p.value[dat_neuron$cluster == 6],
            d_f_p = dat_neuron$d_f_p.value[dat_neuron$cluster == 6]) +
  labs(y = '% of Neurons')

clust_6


ggsave(plot = clust_6,
       file = "clust_6_prop.svg",
       device = "svg",
       units = "in",
       width = 2,
       height = 2,
       path = '/Users/ggraham/Desktop/multiome_poa/Manuscript/Plots/Manuscript v1.2/prop_neurons')


clust_23 = proportion_plotter(joint_neuron, 23,
            d_m_p = dat_neuron$d_m_p.value[dat_neuron$cluster == 23],
            f_m_p = dat_neuron$f_m_p.value[dat_neuron$cluster == 23],
            d_f_p = dat_neuron$d_f_p.value[dat_neuron$cluster == 23]) +
  labs(y = '% of Neurons')

clust_23
