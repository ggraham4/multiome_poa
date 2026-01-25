library(Seurat)
library(emmeans)
library(tidyverse)
define_degs_2 = readRDS('Functions/define_degs2')

obj <- readRDS('~/Desktop/optimal_clustering_rna_only.rds')


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
dat2 = define_degs_2(dat, singular_test = F)
concise = dat2[,c('cluster', 'full_label')]


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
dat_neuron$singular = F
dat2_neuron = define_degs_2(dat_neuron, singular_test = F)
concise_neuron = dat2_neuron[,c('cluster', 'full_label')]

joint_neuron$Status = factor(joint_neuron$Status, levels = c('NRM',
                                                             'M',
                                                             'D',
                                                             'E',
                                                             'NF',
                                                             'F'))
ggplot(subset(joint_neuron, res0.8_50nn_40PC_45LSI==6), aes(x = Status, y = ncells/total_cells))+
  geom_boxplot()+
  geom_point()+
  labs(title = 'Cluster 6')

ggplot(subset(joint_neuron, res0.8_50nn_40PC_45LSI==4), aes(x = Status, y = ncells/total_cells))+
  geom_boxplot()+
  geom_point()


ggplot(subset(joint_neuron, res0.8_50nn_40PC_45LSI==0), aes(x = Status, y = ncells/total_cells))+
  geom_boxplot()+
  geom_point()+
    labs(title = 'Cluster 0')

ggplot(subset(joint_neuron, res0.8_50nn_40PC_45LSI==18), aes(x = Status, y = ncells/total_cells))+
  geom_boxplot()+
  geom_point()+
    labs(title = 'Cluster 18')


#### I think this is too liberal, so I am going to look at proportions using linear ####
dat_neuron2 = data.frame()
for(i in unique(joint_neuron$res0.8_50nn_40PC_45LSI)){
  sub= subset(joint_neuron, res0.8_50nn_40PC_45LSI ==i & Status %in% c('M','D','F'))
  mat = (sub$ncells/sub$total_cells)
  
  mod = lm(mat ~ Status, data = sub)
  p = anova(mod, test = 'Chisq')
  p.value = p$`Pr(>F)`[1]
  
  pair = pairs(emmeans(mod, 'Status'), adjust = 'none')%>%as.data.frame()
  
  newd = data.frame(cluster = i,
                    p = p.value,
                    m_d_p.value = pair$p.value[pair$contrast == 'M - D'],
                    m_f_p.value = pair$p.value[pair$contrast == 'M - F'],
                    d_f_p.value = pair$p.value[pair$contrast == 'D - F'],
                    m_d_estimate = pair$estimate[pair$contrast == 'M - D'],
                    m_f_estimate = pair$estimate[pair$contrast == 'M - F'],
                    d_f_estimate = pair$estimate[pair$contrast == 'D - F']
)
  dat_neuron2 = rbind(newd, dat_neuron2)
}
dat_neuron2$av_q.value = p.adjust(dat_neuron2$p, 'fdr',nrow(dat_neuron2))
dat_neuron2$sig = ifelse(dat_neuron2$av_q.value < 0.05, '*', NA)
dat_neuron2$singular = F
dat2_neuron2 = define_degs_2(dat_neuron2, singular_test = F)
concise_neuron = dat2_neuron2[,c('cluster', 'full_label')]
# ok no significant differences, too conservative
