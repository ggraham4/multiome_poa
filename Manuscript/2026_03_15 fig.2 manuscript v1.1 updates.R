library(Seurat)
library(patchwork)
library(tidyverse)
library(ggplot2)
library(Polychrome)
library(emmeans)
library(ggsignif)
  clown_go = readRDS("Functions/clown_go2")  
library(clusterProfiler)



colors = c("red", "#006400", "blue",'#000000', 'purple','gray','brown','orange')

obj  = readRDS("~/Desktop/optimal_clustering_rna_only.rds")
status_to_phase = c(
  "D"='I',
  'M' = 'M',
  'F' ='F',
  'NF' ='NF',
  'NRM' = 'NRM',
  'E' = 'LI'
)


proportion_plotter = function(data_frame, cluster){
  dag = subset(data_frame, res0.8_50nn_40PC_45LSI==cluster&
                 Status != 'NRM')
  
  dag$Phase = status_to_phase[dag$Status]
  dag$Phase = factor(dag$Phase, levels = c(
    'M',
    'I',
    'LI',
    'NF',
    "F"
  ))
  
  dag$prop = dag$ncells / dag$total_cells
  
  dag_summary = dag %>%
    filter(Phase != "NF") %>%
    group_by(Phase) %>%
    summarise(
      mean_prop = mean(prop, na.rm = TRUE),
      se = sd(prop, na.rm = TRUE) / sqrt(n()),
      .groups = 'drop'
    )
  
  r = ggplot(dag, aes(x = Phase, y = prop))+
    geom_bar(data = dag_summary,
             aes(x = Phase, y = mean_prop, fill = Phase),
             stat = 'identity', inherit.aes = FALSE)+
    geom_errorbar(data = dag_summary,
                  aes(x = Phase, ymin = mean_prop - se, ymax = mean_prop + se),
                  width = 0.4, inherit.aes = FALSE, linewidth = 1)+
    geom_point(size = 1)+
    theme_classic()+
    scale_x_discrete(drop = FALSE)+
    scale_y_continuous(labels = scales::percent)+
    labs(x = 'Phase', y = 'Proportion of ')
  
  return(r)
  
}
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

textsize = 3

olig = proportion_plotter(joint, 2)+
  labs(y = '% of Cells')+
  geom_signif(xmin = c(1),xmax = c(1.9),y_position = c(0.2), annotation =c('***'), color = "black",tip_length = c(0,0), textsize = textsize)+
  geom_signif(xmin = c(1),xmax = c(5),y_position = c(0.25), annotation =c('***'), color = "black",tip_length = c(0,0), textsize = textsize)+
  geom_signif(xmin = c(2.1),xmax = c(5),y_position = c(0.2), annotation =c('***'), color = "black",tip_length = c(0,0), textsize = textsize)+
  ylim(0, 0.32)+
      scale_y_continuous(labels = scales::percent)+
  theme(legend.position = 'none')

olig

ggsave(plot = olig,
       file = "olig_prop.svg",
       device = "svg",
       units = "in",
       width = 2,
       height = 2,
       path = "Manuscript/Plots/Manuscript v1.1/Fig.2")

rgc = proportion_plotter(joint, 1)+
  labs(y = '% of Cells')+
  geom_signif(xmin = c(1),xmax = c(1.9),y_position = c(0.20), annotation =c('**'), color = "black",tip_length = c(0,0), textsize = textsize)+
  geom_signif(xmin = c(1),xmax = c(5),y_position = c(0.27), annotation =c('p = 0.907'), color = "black",tip_length = c(0,0), textsize = textsize)+
  geom_signif(xmin = c(2.1),xmax = c(5),y_position = c(0.20), annotation =c('***'), color = "black",tip_length = c(0,0), textsize = textsize)+
  ylim(0, 0.34)+
      scale_y_continuous(labels = scales::percent)+
  theme(legend.position = 'none')
rgc

ggsave(plot = rgc,
       file = "rgc_prop.svg",
       device = "svg",
       units = "in",
       width = 2,
       height = 2,
       path = "Manuscript/Plots/Manuscript v1.1/Fig.2")

mg = proportion_plotter(joint, 11)+
  labs(y = '% of Cells')+
  geom_signif(xmin = c(1),xmax = c(1.9),y_position = c(0.07), annotation =c('*'), color = "black",tip_length = c(0,0), textsize = textsize)+
  geom_signif(xmin = c(1),xmax = c(5),y_position = c(0.08), annotation =c('**'), color = "black",tip_length = c(0,0), textsize = textsize)+
  geom_signif(xmin = c(2.1),xmax = c(5),y_position = c(0.07), annotation =c('p = 0.347'), color = "black",tip_length = c(0,0), textsize = textsize)+
  ylim(0, 0.32)+
      scale_y_continuous(labels = scales::percent)+
  theme(legend.position = 'none')
mg
ggsave(plot = mg,
       file = "mg_prop.svg",
       device = "svg",
       units = "in",
       width = 2,
       height = 2,
       path = "Manuscript/Plots/Manuscript v1.1/Fig.2")


#### looking at neurons only #####
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


clust_0 = proportion_plotter(joint_neuron, 0)+
  labs(y = '% of Neurons')+
  geom_signif(xmin = c(1),xmax = c(1.9),y_position = c(0.35), annotation =c('*'), color = "black",tip_length = c(0,0), textsize = textsize)+
  geom_signif(xmin = c(1),xmax = c(5),y_position = c(0.4), annotation =c('*'), color = "black",tip_length = c(0,0), textsize = textsize)+
  geom_signif(xmin = c(2.1),xmax = c(5),y_position = c(0.35), annotation =c('***'), color = "black",tip_length = c(0,0), textsize = textsize)+
  ylim(0.2, 0.4)+
      scale_y_continuous(labels = scales::percent)+
  theme(legend.position = 'none'
        )
clust_0

ggsave(plot = clust_0,
       file = "clust_0_prop.svg",
       device = "svg",
       units = "in",
       width = 2,
       height = 2,
       path = "Manuscript/Plots/Manuscript v1.1/Fig.2")



clust_6 = proportion_plotter(joint_neuron, 6)+
  labs(y = '% of Neurons')+
  geom_signif(xmin = c(1),xmax = c(1.9),y_position = c(0.08), annotation =c('p = 0.510'), color = "black",tip_length = c(0,0), textsize = textsize)+
  geom_signif(xmin = c(1),xmax = c(5),y_position = c(0.11), annotation =c('**'), color = "black",tip_length = c(0,0), textsize = textsize)+
  geom_signif(xmin = c(2.1),xmax = c(5),y_position = c(0.08), annotation =c('*'), color = "black",tip_length = c(0,0), textsize = textsize)+
  ylim(0.2, 0.4)+
      scale_y_continuous(labels = scales::percent)+
  theme(legend.position = 'none'
        )
clust_6

ggsave(plot = clust_6,
       file = "clust6_prop.svg",
       device = "svg",
       units = "in",
       width = 2,
       height = 2,
       path = "Manuscript/Plots/Manuscript v1.1/Fig.2")


clust_9 = proportion_plotter(joint_neuron, 22)+
  labs(y = '% of Neurons')+
  geom_signif(xmin = c(1),xmax = c(1.9),y_position = c(0.033), annotation =c('p = 0.439'), color = "black",tip_length = c(0,0), textsize = textsize)+
  geom_signif(xmin = c(1),xmax = c(5),y_position = c(0.04), annotation =c('**'), color = "black",tip_length = c(0,0), textsize = textsize)+
  geom_signif(xmin = c(2.1),xmax = c(5),y_position = c(0.033), annotation =c('**'), color = "black",tip_length = c(0,0), textsize = textsize)+
  ylim(0.2, 0.4)+
      scale_y_continuous(labels = scales::percent)+
  theme(legend.position = 'none'
        )
clust_9

ggsave(plot = clust_9,
       file = "clust_9_prop.svg",
       device = "svg",
       units = "in",
       width = 2,
       height = 2,
       path = "Manuscript/Plots/Manuscript v1.1/Fig.2")



