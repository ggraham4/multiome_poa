library(Seurat)
library(tidyverse)
library(ggplot2)
library(tidyverse)
`%notin%` = Negate(`%in%`)
library(emmeans)

# read in data
obj <- readRDS("~/Desktop/nemo.orig_harmony.integration_all_testd_clusters.rds")
Idents(obj) <- "res0.8_50nn_40PC_45LSI"
#find subclusters
obj_subcluster = FindSubCluster(obj,1, graph.name='harmony.wsnn')

# remove glial cells, we are only interested in neurogenesis
obj_subcluster = subset(obj_subcluster, 
                                sub.cluster %notin% c(
                                  2, # olig
                                  13, # opc
                                  20, # leuko
                                  26, # DG
                                  11, # mg
                                  15, # eg
                                  '1_4', # tanycyte
                                  '1_5' # astrocyte like
                                ) &
                                  Status != 'NRM')


# summarize cells per subcluster
cells_subcluster = obj_subcluster@meta.data%>%
  group_by(individual, Status, sub.cluster)%>%
  summarize(n_cells = n())

# summarize total cells
cells_total = obj_subcluster@meta.data%>%
  group_by(individual)%>%
  summarize(total_cells = n())

# join and get proportions
joined = cells_subcluster%>%
  right_join(cells_total, by = 'individual')%>%
  mutate(proportion = n_cells/total_cells)%>%
  subset(individual != 'GH')

# define all cells not part of this lineage as "other"
joined$cluster_label = ifelse(joined$sub.cluster %in% c('1_1',
                                                        '1_3',
                                                        '1_2',
                                                        '1_0',
                                                        9),
                              joined$sub.cluster,
                              'other_cells'
                              )


# create a proportion matrix
proportion_matrix_data = joined%>%
  group_by(individual)%>%
  dplyr::select(individual,Status, proportion, cluster_label)%>%
  pivot_wider(values_from = 'proportion', names_from  = 'cluster_label', values_fn = sum)

proportion_matrix = as.matrix(proportion_matrix_data[,-c(1:2)])

# create a generic function to populate the transition matrix with transition
# probabilities

transition_probabilities = function(matrix=proportion_matrix){
  
  p.1_1.1_3 = matrix[,'1_1']/matrix[,'1_3']
  p.1_3.1_2 =matrix[,'1_3']/matrix[,'1_2']
  p.1_2.1_0 =matrix[,'1_2']/matrix[,'1_0']   
  p.1_0.9=matrix[,'1_0']/matrix[,'9']  
  p.9.other =matrix[,'9']/matrix[,'other_cells']  
  
   pStay.1_1 = matrix[,'1_1'] / (matrix[,'1_1'] + matrix[,'1_3'])
  pStay.1_3 = matrix[,'1_3'] / (matrix[,'1_3'] + matrix[,'1_2'])
  pStay.1_2 = matrix[,'1_2'] / (matrix[,'1_2'] + matrix[,'1_0'])
  pStay.1_0 = matrix[,'1_0'] / (matrix[,'1_0'] + matrix[,'9'])
  pStay.9   = matrix[,'9']   / (matrix[,'9']   + matrix[,'other_cells'])
  
  data.frame(
    p.1_1.1_3, p.1_3.1_2, p.1_2.1_0, p.1_0.9, p.9.other,
    pStay.1_1, pStay.1_3, pStay.1_2, pStay.1_0, pStay.9
  )
  
  
}
# create a transition matrix
transition_matrix = transition_probabilities()
transition_matrix$individual = proportion_matrix_data$individual
transition_matrix$Status = proportion_matrix_data$Status

transitions_long = transition_matrix%>%
  pivot_longer(cols = starts_with('p'), values_to ='probability', names_to = 'probability_type')
transitions_long$Status = factor(transitions_long$Status, levels = c('M',
                                                                     'D',
                                                                     'E',
                                                                     'NF',
                                                                     'F'))

transitions_long_transitions = transitions_long[!str_detect(pattern = 'pStay', string=transitions_long$probability_type),]

transitions_long_transitions$probability_type = factor(transitions_long_transitions$probability_type, 
                                                       levels = c('p.1_1.1_3',
                                                                  'p.1_3.1_2',
                                                                  'p.1_2.1_0',
                                                                  'p.1_0.9',
                                                                  'p.9.other'))
ggplot(transitions_long_transitions, aes(x = Status, y = probability , color = Status))+
  geom_boxplot()+
  geom_point()+
  facet_wrap(~probability_type, scales = 'free')


transitions_stay = transitions_long[str_detect(pattern = 'pStay', string=transitions_long$probability_type),]

transitions_stay$probability_type = factor(transitions_stay$probability_type, 
                                                       levels = c('pStay.1_1',
                                                                  'pStay.1_3',
                                                                  'pStay.1_2',
                                                                  'pStay.1_0',
                                                                  'pStay.9'))
ggplot(transitions_stay, aes(x = Status, y = probability , color = Status))+
  geom_boxplot()+
  geom_point()+
  facet_wrap(~probability_type, scales = 'free')

manova_matrix = cbind(transition_matrix[,1:(ncol(transition_matrix)-2)])%>%as.matrix()
manova = lm(manova_matrix~transition_matrix$Status)
summary.aov(manova)


model_stay_1_3 = aov(transition_matrix$pStay.1_3~transition_matrix$Status)
summary(model_stay_1_3)