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

ggsave(plot = p,
       file = "UMAP.tiff",
       device = "tiff",
       units = "in",
       width = 4,
       height = 4,
       path = "Manuscript/Plots/Fig.2")

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

ggsave(plot = degs,
       file = "degs.svg",
       device = "svg",
       units = "in",
       width = 8.5,
       height = 2.5,
       path = "Manuscript/Plots/Fig.2")

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


ggsave(plot = olig,
       file = "olig_prop.svg",
       device = "svg",
       units = "in",
       width = 1.5,
       height = 1.5,
       path = "Manuscript/Plots/Fig.2")

rgc = proportion_plotter(joint, 1)+
  labs(y = '% of Cells')+
  geom_signif(xmin = c(1),xmax = c(1.9),y_position = c(0.20), annotation =c('**'), color = "black",tip_length = c(0,0), textsize = textsize)+
  geom_signif(xmin = c(1),xmax = c(5),y_position = c(0.27), annotation =c('p = 0.907'), color = "black",tip_length = c(0,0), textsize = textsize)+
  geom_signif(xmin = c(2.1),xmax = c(5),y_position = c(0.20), annotation =c('***'), color = "black",tip_length = c(0,0), textsize = textsize)+
  ylim(0, 0.34)+
      scale_y_continuous(labels = scales::percent)+
  theme(legend.position = 'none', axis.title.y = element_blank())
rgc

ggsave(plot = rgc,
       file = "rgc_prop.svg",
       device = "svg",
       units = "in",
       width = 1.35,
       height = 1.5,
       path = "Manuscript/Plots/Fig.2")

mg = proportion_plotter(joint, 11)+
  labs(y = '% of Cells')+
  geom_signif(xmin = c(1),xmax = c(1.9),y_position = c(0.07), annotation =c('*'), color = "black",tip_length = c(0,0), textsize = textsize)+
  geom_signif(xmin = c(1),xmax = c(5),y_position = c(0.08), annotation =c('**'), color = "black",tip_length = c(0,0), textsize = textsize)+
  geom_signif(xmin = c(2.1),xmax = c(5),y_position = c(0.07), annotation =c('p = 0.347'), color = "black",tip_length = c(0,0), textsize = textsize)+
  ylim(0, 0.32)+
      scale_y_continuous(labels = scales::percent)+
  theme(legend.position = 'none', axis.title.y = element_blank())
mg
ggsave(plot = mg,
       file = "mg_prop.svg",
       device = "svg",
       units = "in",
       width = 1.35,
       height = 1.5,
       path = "Manuscript/Plots/Fig.2")


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

ggsave(plot = clust_0,
       file = "clust_0_prop.svg",
       device = "svg",
       units = "in",
       width = 1.5,
       height = 1.5,
       path = "Manuscript/Plots/Fig.2")


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

ggsave(plot = clust_6,
       file = "clust6_prop.svg",
       device = "svg",
       units = "in",
       width = 1.35,
       height = 1.5,
       path = "Manuscript/Plots/Fig.2")

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

ggsave(plot = clust_24,
       file = "clust_24_prop.svg",
       device = "svg",
       units = "in",
       width = 1.35,
       height = 1.5,
       path = "Manuscript/Plots/Fig.2")



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
  
ggsave(plot = all_degs_6_go,
       file = "all_degs_6_go.svg",
       device = "svg",
       units = "in",
       width = 4,
       height = 2,
       path = "Manuscript/Plots/Fig.2")


clown_go(degs$gene[degs$cluster==6 & degs$second_word=='Upregulated'])%>%dotplot()
#clown_go(degs$gene[degs$cluster==6 & degs$second_word=='Downregulated'])%>%dotplot()

up_degs_6_go = clown_go(degs$gene[degs$cluster==6& degs$second_word=='Upregulated'])%>%dotplot()+
  labs(title = '6_POA_2 Upregulated DEGs')+
  theme(axis.text.y = element_text(size = 4),
        text = element_text(size =4),
        axis.text.x = element_text(size = 4),
        axis.title.x = element_text(size=4))

ggsave(plot = up_degs_6_go,
       file = "up_degs_6_go.svg",
       device = "svg",
       units = "in",
       width = 4,
       height = 2,
       path = "Manuscript/Plots/Fig.2")

clown_go(degs$gene[degs$cluster==1 ])%>%dotplot()
clown_go(degs$gene[degs$cluster==1 & degs$second_word=='Downregulated'])%>%dotplot()
clown_go(degs$gene[degs$cluster==1 & degs$second_word=='Upregulated'])%>%dotplot()

clown_go(degs$gene[degs$cluster==1 ],p = 1, fdr =1)%>%dotplot()

        