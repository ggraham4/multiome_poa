{
  library(ggsignif)
  library(patchwork)
  library(ggplot2)
  library(Seurat)
  library(dplyr)
  library(tidyverse)
  library(CytoTRACE)
  library(BiocManager)
  library(lme4)
  library(enrichR)
  library(WGCNA)
  `%notin%` = Negate(`%in%`)
  library(car)
  library(emmeans)
  library(glmGamPoi)
  library(scran)
  library(parallel)
    library(parallel)
  library(reshape2)
library(ggplot2)
  library(Seurat)
  library(tidyr)
  library(lme4)
  library(dplyr)
  library(MASS)
  library(Signac)
  library('glmGamPoi')
  library(scran)
  library(emmeans)
  library(openxlsx)
  library(ggplot2)
  library(stringr)
  library(forcats)
  library(clusterProfiler)
  library(biomaRt)
  library(hdWGCNA)
  
}
obj = readRDS("~/Desktop/optimal_clustering_rna_only.rds")


measures = read.csv('Measures/all_data.csv')
measures$Time_Day_2= as.numeric(measures$Time_Day_2)
measures$Behaviors_Day_2= as.numeric(measures$Behaviors_Day_2)
measures$Log_11KT = as.numeric(measures$Log_11KT)

status_to_phase = list('D'='I',
                       'E' = 'LI',
                       'EP' = 'LIP',
                       'S' = 'IP',
                       'M'='M',
                       'F'='F',
                       'NF'='NF',
                       'NM'='NM')
measures$Phase = unlist(status_to_phase[measures$Status])
measures$Phase = factor(measures$Phase, levels =c('F','M','I','IP','LI','LIP','NF','NM'))
measures$Dominance = ifelse(measures$Phase %in%c('F','I','LI','NF'), 'Dominant', 'Subordinate')
measures$Behaviors_Day_2 = as.numeric(measures$Behaviors_Day_2)



obj@meta.data$Status = factor(obj@meta.data$Status, levels = c('NRM',
                                                               'M',
                                                               'D',
                                                               'E',
                                                               'NF',
                                                               'F'))
neur = obj@meta.data%>%
  group_by(individual, Status)%>%
  summarize(n_cells = n())

joint=measures%>%right_join(neur, by = join_by('Fish'=='individual'))%>%
  subset(!is.na(n_cells)&
           Phase!='NRM')
  
joint$Phase = factor(joint$Phase, levels = c('M','I','LI','NF','F'))

ggplot(joint, aes(x = Phase, y = residuals_len))+
  geom_boxplot()+
  geom_point()+
  labs(y = 'Neurons Corrected for Length')
  
neurons_size = lm(n_cells~mass_final_cm, data = joint)
joint$residuals_mass = resid(neurons_size)

  
ggplot(joint, aes(x = Phase, y = residuals_mass))+
  geom_boxplot()+
  geom_point()+
  labs(y = 'Neurons Corrected for Mass')

neurons_len = lm(n_cells~length_final_cm, data = joint)
joint$residuals_len = resid(neurons_len)

  
ggplot(joint, aes(x = Phase, y = residuals_len))+
  geom_boxplot()+
  geom_point()+
  labs(y = 'Neurons Corrected for Length')


# proportion of cells that are neurons
obj@meta.data$neuron = ifelse(!obj@meta.data$final_clusters %in% c(1,2,11,13,20,15,26 ), 'neuron', 'not')

obj@meta.data$rgc = ifelse(obj@meta.data$final_clusters %in% c(1), 'rgc', 'other')
obj@meta.data$mic = ifelse(obj@meta.data$final_clusters %in% c(11), 'mic', 'other')
obj@meta.data$olig = ifelse(obj@meta.data$final_clusters %in% c(2), 'olig', 'other')
obj@meta.data$div = ifelse(obj@meta.data$final_clusters %in% c(26), 'div', 'other')
obj@meta.data$epe = ifelse(obj@meta.data$final_clusters %in% c(15), 'epe', 'other')

obj@meta.data$clust_0 = ifelse(obj@meta.data$final_clusters %in% c('0'), 'zero', 'other')
obj@meta.data$clust_3 = ifelse(obj@meta.data$final_clusters %in% c('3'), 'zero', 'other')
obj@meta.data$clust_6 = ifelse(obj@meta.data$final_clusters %in% c('6'), 'zero', 'other')

prop_neuron =obj@meta.data%>%
  group_by(individual, Status)%>%
  subset(individual != 'GH')%>%
  summarize(n_neuron = sum(neuron =='neuron'),
            prop_neuron = sum(neuron=='neuron')/n())

prop_rgc =obj@meta.data%>%
  group_by(individual, Status)%>%
  subset(individual != 'GH')%>%
  summarize(n_rgc = sum(rgc =='rgc'),
            prop_rgc = sum(rgc=='rgc')/n())

prop_epe =obj@meta.data%>%
  group_by(individual, Status)%>%
  subset(individual != 'GH')%>%
  summarize(n_epe = sum(epe =='epe'),
            prop_epe = sum(epe=='epe')/n())

prop_0 =obj@meta.data%>%
  group_by(individual, Status)%>%
  subset(individual != 'GH')%>%
  summarize(clust_0_ = sum(clust_0 =='zero'),
            prop_clust_0 = sum(clust_0=='zero')/n())
prop_3 =obj@meta.data%>%
  group_by(individual, Status)%>%
  subset(individual != 'GH')%>%
  summarize(clust_3_ = sum(clust_3 =='zero'),
            prop_clust_3= sum(clust_3=='zero')/n())

prop_6 =obj@meta.data%>%
  group_by(individual, Status)%>%
  subset(individual != 'GH')%>%
  summarize(clust_6_ = sum(clust_6 =='zero'),
            prop_clust_6= sum(clust_6=='zero')/n())


prop_div =obj@meta.data%>%
  group_by(individual, Status)%>%
  subset(individual != 'GH')%>%
  summarize(n_div = sum(div =='div'),
            prop_div = sum(div=='div')/n())


prop_olig =obj@meta.data%>%
  group_by(individual, Status)%>%
  subset(individual != 'GH')%>%
  summarize(n_olig = sum(olig =='olig'),
            prop_olig = sum(olig=='olig')/n())


prop_mic =obj@meta.data%>%
  group_by(individual, Status)%>%
  subset(individual != 'GH')%>%
  summarize(n_mic = sum(mic =='mic'),
            prop_mic = sum(mic=='mic')/n())

joint$prop_neuron = prop_neuron$prop_neuron
joint$prop_rgc= prop_rgc$prop_rgc
joint$prop_mic= prop_mic$prop_mic
joint$prop_olig= prop_olig$prop_olig
joint$prop_epe= prop_epe$prop_epe
joint$prop_clust_0= prop_0$prop_clust_0

joint$prop_clust_3= prop_3$prop_clust_3
joint$prop_clust_6= prop_6$prop_clust_6

joint$rgc_neuron = joint$prop_rgc + joint$prop_neuron

joint$div = prop_div$prop_div

ggplot(joint, aes(x = Phase, y = prop_neuron))+
  geom_boxplot()+
  geom_point()+
  labs(y = 'Neurons')

ggplot(joint, aes(x = Phase, y = prop_clust_0))+
  geom_boxplot()+
  geom_point()+
  labs(y = 'Cluster 0')

ggplot(joint, aes(x = Phase, y = prop_clust_3))+
  geom_boxplot()+
  geom_point()+
  labs(y = 'Cluster 3')

ggplot(joint, aes(x = Phase, y = prop_clust_6))+
  geom_boxplot()+
  geom_point()+
  labs(y = 'Cluster 6')



ggplot(joint, aes(x = Phase, y = prop_olig))+
  geom_boxplot()+
  geom_point()+
  labs(y = 'Oligos')

ggplot(joint, aes(x = Phase, y = prop_epe))+
  geom_boxplot()+
  geom_point()+
  labs(y = 'Ependymal')

ggplot(joint, aes(x = Phase, y = div))+
  geom_boxplot()+
  geom_point()+
  labs(y = 'Dividing Glia')

ggplot(joint, aes(x = Phase, y = prop_olig+prop_neuron+prop_rgc+prop_mic+div))+
  geom_boxplot()+
  geom_point()+
  labs(y = 'Oligos + Neuron + Radial Glia + Microglia')


ggplot(joint, aes(x = Phase, y = prop_rgc))+
  geom_boxplot()+
  geom_point()+
  labs(y = 'Radial Glia')

ggplot(joint, aes(x = Phase, y = rgc_neuron))+
  geom_boxplot()+
  geom_point()+
  labs(y = 'Radial Glia + Neuron')

ggplot(joint, aes(x = Phase, y = rgc_neuron))+
  geom_boxplot()+
  geom_point()+
  labs(y = 'Olig + Neuron')

ggplot(joint, aes(x = Phase, y = prop_mic))+
  geom_boxplot()+
  geom_point()+
  labs(y = 'Microglia')


indiv_cells = obj@meta.data%>%
  group_by(individual, Status, final_clusters)%>%
  summarize(n_cells = n())

total_cells = obj@meta.data%>%
  group_by(individual, Status)%>%
  summarize(total_cells = n())

indiv_cells_join = indiv_cells%>%
  right_join(total_cells, by = 'individual')


indiv_cells_join$Status.x = factor(indiv_cells_join$Status.x, levels = c(
                                                               'M',
                                                               'D',
                                                               'E',
                                                               'NF',
                                                               'F'))
ggplot(subset(indiv_cells_join, Status.x %in% c("M",'D','E','F')), aes(x = final_clusters, y = n_cells/total_cells, fill = Status.x))+
  stat_summary(geom = 'bar', fun = 'mean', position = 'dodge')



