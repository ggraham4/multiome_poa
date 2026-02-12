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

p = DimPlot(obj)+
  theme_void()+
  theme(legend.position = "none")
p

#ggsave(plot = p,
 #      file = "UMAP.tiff",
  #     device = "tiff",
   #    units = "in",
    #   width = 4,
     #  height = 4,
      # path = "Manuscript/Plots/Fig.2")

DimPlot(obj, label = T)

### DEGs
degs =read.csv('/Users/ggraham/Desktop/multiome_poa/DEG Outputs/FINAL degs classified w singular.csv')

plot_degs = degs%>%
  group_by(cluster, short_label)%>%
  summarize(n = n())

degs = ggplot(plot_degs, aes(x = (as.numeric(cluster)), y = n, fill = short_label))+
  geom_bar(stat= 'identity')+
  theme_minimal()+
  scale_fill_manual(values = colors)+
  labs(x = 'Cluster', y = 'DEGs', fill = 'Classification')+
  theme(axis.text.x = element_text(angle = 90))+
  scale_x_continuous(breaks = c(0:26))
degs

#ggsave(plot = degs,
 #      file = "degs.svg",
  #     device = "svg",
   #    units = "in",
    #   width = 8.5,
     #  height = 2.5,
      # path = "Manuscript/Plots/Fig.2")

## proportions

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
  
  r = ggplot(dag, aes(x = Phase, y = ncells/total_cells))+
    geom_boxplot(data = subset(dag, Phase != "NF"),
                 aes(x = Phase, y =ncells/total_cells ,
                     fill = Phase), outlier.shape = NA)+
    geom_point(size =0.5)+
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


#ggsave(plot = olig,
 #      file = "olig_prop.svg",
  #     device = "svg",
   #    units = "in",
    #   width = 1.5,
     #  height = 1.5,
      # path = "Manuscript/Plots/Fig.2")

rgc = proportion_plotter(joint, 1)+
  labs(y = '% of Cells')+
  geom_signif(xmin = c(1),xmax = c(1.9),y_position = c(0.20), annotation =c('**'), color = "black",tip_length = c(0,0), textsize = textsize)+
  geom_signif(xmin = c(1),xmax = c(5),y_position = c(0.27), annotation =c('p = 0.907'), color = "black",tip_length = c(0,0), textsize = textsize)+
  geom_signif(xmin = c(2.1),xmax = c(5),y_position = c(0.20), annotation =c('***'), color = "black",tip_length = c(0,0), textsize = textsize)+
  ylim(0, 0.34)+
      scale_y_continuous(labels = scales::percent)+
  theme(legend.position = 'none', axis.title.y = element_blank())
rgc

#ggsave(plot = rgc,
 #      file = "rgc_prop.svg",
  #     device = "svg",
   #    units = "in",
    #   width = 1.35,
     #  height = 1.5,
      # path = "Manuscript/Plots/Fig.2")

mg = proportion_plotter(joint, 11)+
  labs(y = '% of Cells')+
  geom_signif(xmin = c(1),xmax = c(1.9),y_position = c(0.07), annotation =c('*'), color = "black",tip_length = c(0,0), textsize = textsize)+
  geom_signif(xmin = c(1),xmax = c(5),y_position = c(0.08), annotation =c('**'), color = "black",tip_length = c(0,0), textsize = textsize)+
  geom_signif(xmin = c(2.1),xmax = c(5),y_position = c(0.07), annotation =c('p = 0.347'), color = "black",tip_length = c(0,0), textsize = textsize)+
  ylim(0, 0.32)+
      scale_y_continuous(labels = scales::percent)+
  theme(legend.position = 'none', axis.title.y = element_blank())
mg
#ggsave(plot = mg,
 #      file = "mg_prop.svg",
  #     device = "svg",
   #    units = "in",
    #   width = 1.35,
     #  height = 1.5,
      # path = "Manuscript/Plots/Fig.2")


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
  geom_signif(xmin = c(1),xmax = c(5),y_position = c(0.37), annotation =c('*'), color = "black",tip_length = c(0,0), textsize = textsize)+
  geom_signif(xmin = c(2.1),xmax = c(5),y_position = c(0.35), annotation =c('***'), color = "black",tip_length = c(0,0), textsize = textsize)+
  ylim(0.2, 0.4)+
      scale_y_continuous(labels = scales::percent)+
  theme(legend.position = 'none',# axis.title.y = element_blank()
        )
clust_0

#ggsave(plot = clust_0,
 #      file = "clust_0_prop.svg",
  #     device = "svg",
   #    units = "in",
    #   width = 1.5,
     #  height = 1.5,
      # path = "Manuscript/Plots/Fig.2")


clust_6 = proportion_plotter(joint_neuron, 6)+
  labs(y = '% of Neurons')+
  geom_signif(xmin = c(1),xmax = c(1.9),y_position = c(0.08), annotation =c('p = 0.510'), color = "black",tip_length = c(0,0), textsize = textsize)+
  geom_signif(xmin = c(1),xmax = c(5),y_position = c(0.11), annotation =c('**'), color = "black",tip_length = c(0,0), textsize = textsize)+
  geom_signif(xmin = c(2.1),xmax = c(5),y_position = c(0.08), annotation =c('*'), color = "black",tip_length = c(0,0), textsize = textsize)+
  ylim(0.2, 0.4)+
      scale_y_continuous(labels = scales::percent)+
  theme(legend.position = 'none', axis.title.y = element_blank()
        )
clust_6

#ggsave(plot = clust_6,
 #      file = "clust6_prop.svg",
  #     device = "svg",
   #    units = "in",
    #   width = 1.35,
     #  height = 1.5,
      # path = "Manuscript/Plots/Fig.2")

clust_24 = proportion_plotter(joint_neuron, 22)+
  labs(y = '% of Neurons')+
  geom_signif(xmin = c(1),xmax = c(1.9),y_position = c(0.033), annotation =c('p = 0.229'), color = "black",tip_length = c(0,0), textsize = textsize)+
  geom_signif(xmin = c(1),xmax = c(5),y_position = c(0.04), annotation =c('**'), color = "black",tip_length = c(0,0), textsize = textsize)+
  geom_signif(xmin = c(2.1),xmax = c(5),y_position = c(0.033), annotation =c('***'), color = "black",tip_length = c(0,0), textsize = textsize)+
  ylim(0.2, 0.4)+
      scale_y_continuous(labels = scales::percent)+
  theme(legend.position = 'none', axis.title.y = element_blank()
        )
clust_24

#ggsave(plot = clust_24,
 #      file = "clust_24_prop.svg",
  #     device = "svg",
   #    units = "in",
    #   width = 1.35,
     #  height = 1.5,
      # path = "Manuscript/Plots/Fig.2")



# NT type 
nt_type = DotPlot(neurons_only, c('LOC111588076',
                        'gad2',
                        'LOC111584103',
                        'slc17a6b'),
        dot.min = .1,
        scale=F,
        dot.scale = 3)+
  coord_flip()+
  scale_x_discrete(labels =c('gad1b',
                             'gad2',
                             'slc17a6a',
                             'sc17a6b'))+
  labs(y = 'Cluster')+
  theme(axis.title.y = element_blank())

ggsave(plot = nt_type,
       file = "nt_type.svg",
       device = "svg",
       units = "in",
       width = 5.5,
       height = 1.5,
       path = "Manuscript/Plots/Fig.2")


# DEG enrichment? 
degs =read.csv('/Users/ggraham/Desktop/multiome_poa/DEG Outputs/FINAL degs classified w singular.csv')

clown_go = readRDS("Functions/clown_go2")  
clown_go(degs$gene[degs$second_word=='Upregulated'], p= 1, fdr = 1)%>%dotplot()
clown_go(degs$gene[degs$second_word=='Downregulated'], p= 1, fdr = 1)%>%dotplot()
clown_go(degs$gene[degs$second_word=='Downregulated'])%>%dotplot()
clown_go(degs$gene[degs$second_word=='Downregulated'])
clown_go(degs$gene, p= 1, fdr = 1)%>%dotplot()

clown_go(degs$gene[degs$cluster==6])%>%dotplot()


clown_go(degs$gene[degs$second_word=='Upregulated'], p=1)%>%dotplot()

clown_go(degs$gene[degs$first_word=='Early'], p= 1, fdr = 1)%>%dotplot()
clown_go(degs$gene[degs$first_word=='Late'], p= 1, fdr = 1)%>%dotplot()
clown_go(degs$gene[degs$first_word=='Late'])%>%dotplot()

clown_go(degs$gene[degs$first_word=='Transiently'], p= 1, fdr = 1)%>%dotplot()
clown_go(degs$gene[degs$first_word=='Progressively'], p= 1, fdr = 1)%>%dotplot()

clown_go(degs$gene[degs$first_word=='Transiently'])%>%dotplot()
clown_go(degs$gene[degs$first_word=='Progressively'])%>%dotplot()

 table(degs$type[degs$second_word=='Upregulated'])%>%sort()
 table(degs$type[degs$second_word=='Downregulated'])%>%sort()
sum(degs$second_word=='Downregulated')


# cluster 6
all_degs_6_go = clown_go(degs$gene[degs$cluster==6])%>%dotplot()+
  labs(title = '6_POA_2 All DEGs')+
  theme(axis.text.y = element_text(size = 4),
        text = element_text(size =4),
        axis.text.x = element_text(size = 4),
        axis.title.x = element_text(size=4))

degs_6_go_df=clown_go(degs$gene[degs$cluster==6])
d = data.frame(degs_6_go_df@result$Count,
degs_6_go_df@result$Description)
  
#ggsave(plot = all_degs_6_go,
 #      file = "all_degs_6_go.svg",
#       device = "svg",
  #     units = "in",
   #    width = 4,
    #   height = 2,
     #  path = "Manuscript/Plots/Fig.2")


clown_go(degs$gene[degs$cluster==6 & degs$second_word=='Upregulated'])%>%dotplot()
#clown_go(degs$gene[degs$cluster==6 & degs$second_word=='Downregulated'])%>%dotplot()

up_degs_6_go = clown_go(degs$gene[degs$cluster==6& degs$second_word=='Upregulated'])%>%dotplot()+
  labs(title = '6_POA_2 Upregulated DEGs')+
  theme(axis.text.y = element_text(size = 4),
        text = element_text(size =4),
        axis.text.x = element_text(size = 4),
        axis.title.x = element_text(size=4))

#ggsave(plot = up_degs_6_go,
#       file = "up_degs_6_go.svg",
#       device = "svg",
#       units = "in",
#       width = 4,
#       height = 2,
#       path = "Manuscript/Plots/Fig.2")

clown_go(degs$gene[degs$cluster==1 ])%>%dotplot()
clown_go(degs$gene[degs$cluster==1 & degs$second_word=='Downregulated'])%>%dotplot()
clown_go(degs$gene[degs$cluster==1 & degs$second_word=='Upregulated'])%>%dotplot()

clown_go(degs$gene[degs$cluster==1 ],p = 1, fdr =1)%>%dotplot()

##### Classifier #####
mf_data = read.csv('/Users/ggraham/Desktop/multiome_poa/Sex Classifier and Linearity/progress_classifier/logistic_mf_validation_07_07_2025.csv')
data = read_csv("/Users/ggraham/Desktop/multiome_poa/Sex Classifier and Linearity/progress_classifier/logistic_individual_probabilities_07_05_2025.csv")
mean_dominants = data%>%
  subset(status == 'D')%>%
  group_by(cluster)%>%
  summarize(mean_dom =mean(prediction),
            n_degs = mean(n_degs))

data_with_dom_mean = data%>%
  dplyr::select(-c(n_degs))%>%
  right_join(mean_dominants, by = 'cluster')

t_data = data.frame()
for(i in unique(mf_data$cluster)){
  
  data =mf_data%>%subset(cluster == i & status == 'm')

male_p = t.test(data$proba, mu =0)$p.value
data2 =mf_data%>%subset(cluster == i & status == 'f')

female_p=t.test(data2$proba, mu = 1)$p.value 
newd = data.frame(cluster =i, 
                  male_p = male_p,
                  female_p = female_p)
t_data = rbind(t_data, newd)

}
t_data$different = ifelse(t_data$male_p<0.1| t_data$female_p<0.1, '*', NA) ## this is not a very good approach I think

cluster_performance = mf_data %>%
  group_by(cluster, status) %>%
  summarize(
    mean_prob = mean(proba),
    sd_prob = sd(proba),
    n = n(),
    .groups = 'drop'
  ) %>%
  pivot_wider(
    names_from = status, 
    values_from = c(mean_prob, sd_prob, n)
  )

# Keep clusters based on practical thresholds
confident_clusters = cluster_performance %>%
  filter(
    mean_prob_m < 0.10 &      # Males classified as male
    mean_prob_f > 0.9        # Females classified as female
  ) %>%
  pull(cluster)

ggplot(mf_data, aes(x = as.factor(cluster), y = proba, color = status))+
  stat_summary(fun = 'mean', geom = 'crossbar')+
  stat_summary(fun.data = "mean_se", geom = "errorbar", width = 0.2) + # why does the t test get this wrong
  geom_point(color = 'black', aes(shape = status))

ggplot(subset(mf_data, cluster %in% confident_clusters), aes(x = as.factor(cluster), y = proba, color = status))+
  stat_summary(fun = 'mean', geom = 'crossbar')+
  stat_summary(fun.data = "mean_se", geom = "errorbar", width = 0.2) + # why does the t test get this wrong
  geom_point(color = 'black', aes(shape = status))

ggplot(subset(data_with_dom_mean, cluster %in% confident_clusters), aes(x = fct_reorder(as.factor(cluster), mean_dom), y = prediction, color = status, shape = status))+
  stat_summary(geom= 'point', fun='mean', size=3)+
  theme_minimal()+
  labs(x = 'Cluster', y = 'Prediction', title = 'Individuals', shape = 'Status', color = 'Status')

plott =subset(data_with_dom_mean, cluster %in% confident_clusters)%>%
  group_by(status, cluster)%>%
  summarize(mean_prediction = mean(prediction),
            se_prediction = sd(prediction)/sqrt(n()))%>%
  subset(status == 'D')

plott$neurotransmitter = ifelse(plott$cluster %in% c(24,0, 6, 9), 'Mixed',NA)
plott$neurotransmitter = ifelse(plott$cluster %in% c(1,2, 11, 13, 15, 20, 26), 'Glial',plott$neurotransmitter)
plott$neurotransmitter = ifelse(plott$cluster %in% c(19, 3, 14, 25), 'GABAergic',plott$neurotransmitter)
plott$neurotransmitter = ifelse(is.na(plott$neurotransmitter), 'Glutamatergic', plott$neurotransmitter)


library(forcats)
classifier = ggplot(plott, aes(x = fct_reorder(as.factor(cluster), mean_prediction, .desc= T), y = mean_prediction))+
  geom_errorbar(aes(x =fct_reorder(as.factor(cluster), mean_prediction,.desc= T), y =mean_prediction ,ymin =mean_prediction -se_prediction,ymax = mean_prediction+se_prediction ), width = 0.4)+
  geom_point(aes( shape = neurotransmitter), size = 1.5)+
    geom_point(aes( color = neurotransmitter, shape = neurotransmitter), size = 1)+
  labs(y = 'Prediction', x= 'Cluster')+
  theme_minimal()+
  theme(legend.position = 'none')+
  scale_color_manual(values = c('purple', 'black', '#5bb450', 'maroon'))+
  scale_shape_manual(values =c(15, 16, 17, 18))

classifier

ggsave(plot = classifier,
       file = "classifier.svg",
       device = "svg",
       units = "in",
       width = 2.5,
       height = 1.5,
       path = "Manuscript/Plots/Fig.2")

### linearity ####
mecd = readRDS("Functions/mean_expression_cluster_data.rds")

### read in DEGs
deg_data = data.frame()
base_path = "DEG Outputs/05_11_2025 Neg Bin w Doms New_clustering/"
base_dir = dir("DEG Outputs/05_11_2025 Neg Bin w Doms New_clustering")
for(i in base_dir){
  data = read.csv(paste0(base_path, i))
  
  subset_data = subset(data, av_q.value < 0.1) #### EDIT MADE HERE, SINGULAR IN ERROR
  
  data_to_append = subset_data%>%
    dplyr::select(gene, cluster)
  deg_data = rbind(deg_data, data_to_append)
  
}

#summarize deg data
expression_data_all_clusters = list()
for(i in unique(deg_data$cluster)){
  print(i)
  deg_data_subset = subset(deg_data, cluster ==i)
  expression_data = data.frame()
    for(g in unique(deg_data_subset$gene)){
      gex = mecd(object = obj,gene = g,cluster = i, clustering='res0.8_50nn_40PC_45LSI')
      
      gex$gene = g
      expression_data = rbind(expression_data, gex)
    }
  
  expression_data_pivoted = expression_data%>%
    dplyr::select(individual, Status, mean, gene)%>%
    pivot_wider(names_from = gene, values_from=mean)
  
  expression_data_all_clusters[[paste0('cluster_',i)]] = expression_data_pivoted
  
}

## run PCA
pca_data = data.frame()
for(i in unique(deg_data$cluster)){
  print(i)
  
  cluster_data = expression_data_all_clusters[[paste0('cluster_',i)]]
  
  pca_matrix = as.matrix(cluster_data[,3:ncol(cluster_data)])
  
  if(ncol(pca_matrix)<2){next}
  pca = prcomp(pca_matrix, scale = T)
  
  pc_1 = pca$x[,1]
  pc_2 = pca$x[,2]
  
  loading_data = data.frame('individual' = cluster_data$individual,
                            'status' = cluster_data$Status,
                            'PC1' = pc_1,
                            'PC2' = pc_2,
                            'cluster' = i)
  pca_data = rbind(pca_data, loading_data)
  
}

#fviz_pca_ind(pca, habillage = cluster_data$Status)

## do linearity calculation for EACH INDIVIDUAL
linearity_individual_data = data.frame()

for(i in unique(pca_data$cluster)){
  subset_data = subset(pca_data, cluster == i)
  
  # Calculate male and female centroids (these remain as endpoints)
  male_data = subset(subset_data, status =='M')
  female_data = subset(subset_data, status =='F')
  
  male_centroid = c(mean(male_data$PC1), mean(male_data$PC2))
  female_centroid = c(mean(female_data$PC1), mean(female_data$PC2))
  
  # Process each  dominant animal
  dominant_data = subset(subset_data, status =='D')
  if(nrow(dominant_data) > 0){
    for(j in 1:nrow(dominant_data)){
      individual_point = c(dominant_data$PC1[j], dominant_data$PC2[j])
      
      # Calculate distances
      distance_matrix = matrix(c(male_centroid,
                                female_centroid,
                                individual_point),
                              3, 2, byrow = T)
      
      m_f_distance = stats::dist(distance_matrix[c(1,2),], method = 'euclidean')[1]
      m_ind_distance = stats::dist(distance_matrix[c(1,3),], method = 'euclidean')[1]
      f_ind_distance = stats::dist(distance_matrix[c(2,3),], method = 'euclidean')[1]
      
      linearity_index = m_f_distance/(m_ind_distance + f_ind_distance)
      
      newd = data.frame(cluster = i,
                       individual = dominant_data$individual[j],
                       status = 'D',
                       m_f_distance = m_f_distance,
                       m_ind_distance = m_ind_distance,
                       f_ind_distance = f_ind_distance,
                       linearity_index = linearity_index,
                       progress = f_ind_distance/m_f_distance)
      
      linearity_individual_data = rbind(linearity_individual_data, newd)
    }
  }
  
  # Process each  expanded animal
  expanded_data = subset(subset_data, status =='E')
  if(nrow(expanded_data) > 0){
    for(j in 1:nrow(expanded_data)){
      individual_point = c(expanded_data$PC1[j], expanded_data$PC2[j])
      
      # Calculate distances
      distance_matrix = matrix(c(male_centroid,
                                female_centroid,
                                individual_point),
                              3, 2, byrow = T)
      
      m_f_distance = stats::dist(distance_matrix[c(1,2),], method = 'euclidean')[1]
      m_ind_distance = stats::dist(distance_matrix[c(1,3),], method = 'euclidean')[1]
      f_ind_distance = stats::dist(distance_matrix[c(2,3),], method = 'euclidean')[1]
      
      linearity_index = m_f_distance/(m_ind_distance + f_ind_distance)
      
      newd = data.frame(cluster = i,
                       individual = expanded_data$individual[j],
                       status = 'E',
                       m_f_distance = m_f_distance,
                       m_ind_distance = m_ind_distance,
                       f_ind_distance = f_ind_distance,
                       linearity_index = linearity_index,
                       progress = f_ind_distance/m_f_distance)
      
      linearity_individual_data = rbind(linearity_individual_data, newd)
    }
  }
  
  # Process each nf
  nf_data = subset(subset_data, status =='NF')
  if(nrow(nf_data) > 0){
    for(j in 1:nrow(nf_data)){
      individual_point = c(nf_data$PC1[j], nf_data$PC2[j])
      
      # Calculate distances
      distance_matrix = matrix(c(male_centroid,
                                female_centroid,
                                individual_point),
                              3, 2, byrow = T)
      
      m_f_distance = stats::dist(distance_matrix[c(1,2),], method = 'euclidean')[1]
      m_ind_distance = stats::dist(distance_matrix[c(1,3),], method = 'euclidean')[1]
      f_ind_distance = stats::dist(distance_matrix[c(2,3),], method = 'euclidean')[1]
      
      linearity_index = m_f_distance/(m_ind_distance + f_ind_distance)
      
      newd = data.frame(cluster = i,
                       individual = nf_data$individual[j],
                       status = 'NF',
                       m_f_distance = m_f_distance,
                       m_ind_distance = m_ind_distance,
                       f_ind_distance = f_ind_distance,
                       linearity_index = linearity_index,
                       progress = f_ind_distance/m_f_distance)
      
      linearity_individual_data = rbind(linearity_individual_data, newd)
    }
  }
}

summary_data = linearity_individual_data %>%
  group_by(cluster, status) %>%
  summarise(
    n = n(),
    mean_linearity = mean(linearity_index),
    se_linearity = sd(linearity_index)/sqrt(n()),
    mean_progress = mean(progress),
    se_progress = sd(progress)/sqrt(n()),
    .groups = 'drop'
  )

summary_data$neurotransmitter = ifelse(summary_data$cluster %in% c(24,0, 6, 9), 'Mixed',NA)
summary_data$neurotransmitter = ifelse(summary_data$cluster %in% c(1,2, 11, 13, 15, 20, 26), 'Glial',summary_data$neurotransmitter)
summary_data$neurotransmitter = ifelse(summary_data$cluster %in% c(19, 3, 14, 25), 'GABAergic',summary_data$neurotransmitter)
summary_data$neurotransmitter = ifelse(is.na(summary_data$neurotransmitter), 'Glutamatergic', summary_data$neurotransmitter)

linearity =ggplot(subset(summary_data, status == 'D'), aes(x = fct_reorder(as.factor(cluster), .x=mean_linearity,.desc = T), y = mean_linearity)) +
  geom_errorbar(aes(ymin = mean_linearity - se_linearity, 
                    ymax = mean_linearity + se_linearity),
                width = 0.4) +
      geom_point( aes(shape = neurotransmitter), size = 1.5) +
    geom_point(aes(color =neurotransmitter, shape = neurotransmitter), size = 1) +
  ylim(0, 1) +
  labs(x = "Cluster", y = "Linearity Index") +
  theme_minimal()+
  theme(legend.position = 'none')+
  scale_color_manual(values = c('purple', 'black', '#5bb450', 'maroon'))+
  scale_shape_manual(values =c(15, 16, 17, 18))
linearity

ggsave(plot = linearity,
       file = "linearity.svg",
       device = "svg",
       units = "in",
       width = 3.5,
       height = 1.5,
       path = "Manuscript/Plots/Fig.2")


####linearity and classifier correlation #####
lin_clas = summary_data%>%
  subset(status=='D')%>%
  right_join(plott,by = 'cluster')

corr = ggplot(lin_clas, aes(x = mean_linearity, y = mean_prediction))+
  geom_errorbar(aes(x = mean_linearity, y = mean_prediction, 
                    ymin = mean_prediction-se_prediction,
                    ymax = mean_prediction+se_prediction), width = 0.0)+
    geom_errorbarh(aes(x = mean_linearity, y = mean_prediction, 
                    xmin = mean_linearity-se_linearity,
                    xmax = mean_linearity+se_linearity), height = 0.0)+
  theme_minimal()+
      geom_point(aes(shape = neurotransmitter.x), size =1.5)+
    geom_point(aes(color = neurotransmitter.x, shape = neurotransmitter.x), size =1)+
    scale_color_manual(values = c('purple', 'black', '#5bb450', 'maroon'))+
  scale_shape_manual(values =c(15, 16, 17, 18))+
  labs(x = 'Mean +/- SE Linearity', y = 'Mean +/- SE Prediction')+
  geom_text(aes(label = cluster, hjust =1,
                vjust = 1))+
  theme(legend.position ='none')
corr

ggsave(plot = corr,
       file = "linearity_classifier_corr.svg",
       device = "svg",
       units = "in",
       width = 2.5,
       height = 2,
      path = "Manuscript/Plots/Fig.2")


### degs of interest, I need to write a function for this bs I think

degs_6 = subset(degs, cluster == 6)

p_value_annotation <- function(p) {
  if (p < 0.001) {
    return("***")
  } else if (p < 0.01) {
    return("**")
  } else if (p < 0.05) {
    return("*")
  } else {
    return(paste0("p = ", round(p, 3)))
  }
}
degs =read.csv('/Users/ggraham/Desktop/multiome_poa/DEG Outputs/FINAL degs classified w singular.csv')

deg = 'pgr'
data = degs

obj@meta.data$Phase = obj@meta.data$Status
obj@meta.data$Phase= ifelse(obj@meta.data$Phase == 'D', 'I', obj@meta.data$Phase)
obj@meta.data$Phase= ifelse(obj@meta.data$Phase == 'E', 'LI', obj@meta.data$Phase)

obj@meta.data$Phase = factor(obj@meta.data$Phase, levels = c('NRM',
                                                             'M',
                                                             'I',
                                                             'LI',
                                                             'NF',
                                                             'F'))

plotter_function = function(data, deg){
  
  deg_data = subset(data, gene == deg)
  d_f_p.value = deg_data$d_f_p.value
  d_m_p.value = deg_data$d_m_p.value
  f_m_p.value = deg_data$f_m_p.value
  
  d_f_annotation = p_value_annotation(d_f_p.value)
  d_m_annotation = p_value_annotation(d_m_p.value)
  f_m_annotation = p_value_annotation(f_m_p.value)
  
  meta.data = subset(obj@meta.data, final_clusters ==6)
  meta.data$expression = obj@assays$RNA$data[deg,obj@meta.data$final_clusters == 6]
  
  meta_grouped = meta.data%>%
    group_by(Phase, individual)%>%
    summarize(mean_expression = mean(expression))%>%
    subset(Phase != 'NRM')
  
  meta_grouped$Phase = factor(meta_grouped$Phase, levels = c(
                                                             'M',
                                                             'I',
                                                             'LI',
                                                             'NF',
                                                             'F'))

  
    r = ggplot(meta_grouped, aes(x = Phase, y = mean_expression))+
    geom_boxplot(data = subset(meta_grouped, Phase != "NF"),
                 aes(x = Phase, y =mean_expression ,
                     fill = Phase), outlier.shape = NA)+
    geom_point(size =0.5)+
    theme_classic()+
        scale_x_discrete(drop = FALSE)+  
    labs(x = 'Phase', y = 'Normalized Expression')+
    geom_signif(xmin = c(1),
                xmax = c(1.9),
                y_position = c(max(meta_grouped$mean_expression*1.1)),
                annotation =c(d_m_annotation),
                color = "black",
                tip_length = c(0,0), textsize = 3)+
          geom_signif(xmin = c(1),
                xmax = c(5),
                y_position = c(max(meta_grouped$mean_expression*1.3)),
                annotation =c(f_m_annotation),
                color = "black",
                tip_length = c(0,0), textsize = 3)+
    geom_signif(xmin = c(2.1),
                xmax = c(5),
                y_position = c(max(meta_grouped$mean_expression*1.1)),
                annotation =c(d_f_annotation),
                color = "black",
                tip_length = c(0,0), textsize = 3)+
      theme(legend.position = 'none', axis.title.y = element_blank())
    r
  return(r)
}

for(i in data$gene[data$cluster=='6']){
  plot = plotter_function(data, i)
  
  #ggsave(plot = plot,
  #     file = paste0('deg_', i, '.svg'),
   #    device = "svg",
    #   units = "in",
     #  width = 1.35,
      # height = 1.5,
       #path = "Manuscript/Plots/Fig.3/degs_6")

  
}







