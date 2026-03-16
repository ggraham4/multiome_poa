library(Seurat)
library(patchwork)
library(tidyverse)
library(ggplot2)
library(Polychrome)
library(emmeans)
library(ggsignif)
library(clusterProfiler)

obj  = readRDS("~/Desktop/optimal_clustering_rna_only.rds")

sub_1 = FindSubCluster(obj,1, graph.name='harmony.wsnn')
Idents(sub_1) <- 'sub.cluster'
sub_1 = subset(sub_1, final_clusters ==1)
sub_1$Status = factor(sub_1$Status, levels = c('NRM','M',"D",'E','NF','F'))

dimplot = DimPlot(sub_1, label = F)+
  theme_void()+
  theme(legend.position = "none")
dimplot
DimPlot(sub_1, label = T)


### proportions ####
cells_ind = sub_1@meta.data%>%
  group_by(individual)%>%
  summarize(n_cells = n())

cells_sub_ind = sub_1@meta.data%>%
  group_by(individual, Status, sub.cluster)%>%
  summarize(n_cells_in = n())

cells_total = cells_ind%>%
  right_join(cells_sub_ind, by = 'individual')
cells_total$prop = cells_total$n_cells_in/cells_total$n_cells

ggplot(subset(cells_total, Status %in% c('M','D',"F")), aes(x = sub.cluster, y = prop, color = Status))+
  geom_boxplot(alpha = 0)+
  geom_point(position = position_jitterdodge(1))

ggplot(subset(cells_total, sub.cluster =='1_0'), aes(x = sub.cluster, y = prop, color = Status))+
  geom_boxplot(alpha = 0)+
  geom_point(position = position_jitterdodge(1))

ggplot(subset(cells_total, sub.cluster =='1_1'), aes(x = sub.cluster, y = prop, color = Status))+
  geom_boxplot(alpha = 0)+
  geom_point(position = position_jitterdodge(1))

out_df = data.frame()
for(i in 0:5){
cells = subset(cells_total, sub.cluster == paste0('1_',i) & Status !='NRM')
matrix_1 = cbind(cells$n_cells_in,cells$n_cells-cells$n_cells_in)

matrix_1 = glm(matrix_1~Status, family = binomial(), data = cells)
av_ = car::Anova(matrix_1, 3)
pairs = pairs(emmeans(matrix_1, 'Status'), adjust ='none')%>%as.data.frame()

newd = data.frame(cluster = i,
                  av_p = av_$`Pr(>Chisq)`[1],
                 d_m_p =  pairs$p.value[pairs$contrast == 'M - D'],
                 f_m_p =  pairs$p.value[pairs$contrast == 'M - F'],
                 d_f = pairs$p.value[pairs$contrast == 'D - F'])
out_df = rbind(out_df, newd)

}
out_df$signif = ifelse(out_df$av_p<0.05, '*', NA)

cells_total$Phase = as.character(cells_total$Status)
cells_total$Phase = ifelse(cells_total$Status=='D', 'I', cells_total$Phase)
cells_total$Phase = ifelse(cells_total$Phase=='E', 'LI', cells_total$Phase)
cells_total$Phase = factor(cells_total$Phase, levels = c('M', 'I', "LI", 'NF','F'))

# helper to make summary for bars
make_summary = function(dat) {
  dat %>%
    filter(Phase != "NF") %>%
    group_by(Phase) %>%
    summarise(
      mean_prop = mean(n_cells_in / n_cells, na.rm = TRUE),
      se = sd(n_cells_in / n_cells, na.rm = TRUE) / sqrt(n()),
      .groups = 'drop'
    )
}

textsize = 3

# 1_0
cells_total_1_0 = subset(cells_total, sub.cluster == paste0('1_',0) & Status !='NRM')
summary_1_0 = make_summary(cells_total_1_0)

prop_1_0 = ggplot(cells_total_1_0, aes(x = Phase, y = n_cells_in/n_cells))+
  geom_bar(data = summary_1_0,
           aes(x = Phase, y = mean_prop, fill = Phase),
           stat = 'identity', inherit.aes = FALSE)+
  geom_errorbar(data = summary_1_0,
                aes(x = Phase, ymin = mean_prop - se, ymax = mean_prop + se),
                width = 0.2, inherit.aes = FALSE)+
  geom_point(size = 0.5)+
  theme_classic()+
  scale_x_discrete(drop = FALSE)+
  scale_y_continuous(labels = scales::percent)+
  labs(x = 'Phase', y = '% of 1_RGC')+
  geom_signif(xmin = c(1),xmax = c(1.9),y_position = c(0.55), annotation =c('p=0.110'), color = "black",tip_length = c(0,0), textsize = textsize)+
  geom_signif(xmin = c(1),xmax = c(5),y_position = c(0.6), annotation =c('**'), color = "black",tip_length = c(0,0), textsize = textsize)+
  geom_signif(xmin = c(2.1),xmax = c(5),y_position = c(0.55), annotation =c('***'), color = "black",tip_length = c(0,0), textsize = textsize)+
  ylim(0, 0.32)+
  scale_y_continuous(labels = scales::percent)+
  theme(legend.position = 'none')

# 1_1
cells_total_1_1 = subset(cells_total, sub.cluster == paste0('1_',1) & Status !='NRM')
summary_1_1 = make_summary(cells_total_1_1)

prop_1_1 = ggplot(cells_total_1_1, aes(x = Phase, y = n_cells_in/n_cells))+
  geom_bar(data = summary_1_1,
           aes(x = Phase, y = mean_prop, fill = Phase),
           stat = 'identity', inherit.aes = FALSE)+
  geom_errorbar(data = summary_1_1,
                aes(x = Phase, ymin = mean_prop - se, ymax = mean_prop + se),
                width = 0.2, inherit.aes = FALSE)+
  geom_point(size = 0.5)+
  theme_classic()+
  scale_x_discrete(drop = FALSE)+
  scale_y_continuous(labels = scales::percent)+
  labs(x = 'Phase', y = '% of 1_RGC')+
  geom_signif(xmin = c(1),xmax = c(1.9),y_position = c(0.4), annotation =c('p=0.110'), color = "black",tip_length = c(0,0), textsize = textsize)+
  geom_signif(xmin = c(1),xmax = c(5),y_position = c(0.45), annotation =c('**'), color = "black",tip_length = c(0,0), textsize = textsize)+
  geom_signif(xmin = c(2.1),xmax = c(5),y_position = c(0.4), annotation =c('***'), color = "black",tip_length = c(0,0), textsize = textsize)+
  ylim(0, 0.32)+
  scale_y_continuous(labels = scales::percent)+
  theme(legend.position = 'none')

# 1_2
cells_total_1_2 = subset(cells_total, sub.cluster == paste0('1_',2) & Status !='NRM')
summary_1_2 = make_summary(cells_total_1_2)

prop_1_2 = ggplot(cells_total_1_2, aes(x = Phase, y = n_cells_in/n_cells))+
  geom_bar(data = summary_1_2,
           aes(x = Phase, y = mean_prop, fill = Phase),
           stat = 'identity', inherit.aes = FALSE)+
  geom_errorbar(data = summary_1_2,
                aes(x = Phase, ymin = mean_prop - se, ymax = mean_prop + se),
                width = 0.2, inherit.aes = FALSE)+
  geom_point(size = 0.5)+
  theme_classic()+
  scale_x_discrete(drop = FALSE)+
  scale_y_continuous(labels = scales::percent)+
  labs(x = 'Phase', y = '% of 1_RGC')+
  geom_signif(xmin = c(1),xmax = c(1.9),y_position = c(0.32), annotation =c('p=0.662'), color = "black",tip_length = c(0,0), textsize = textsize)+
  geom_signif(xmin = c(1),xmax = c(5),y_position = c(0.38), annotation =c('***'), color = "black",tip_length = c(0,0), textsize = textsize)+
  geom_signif(xmin = c(2.1),xmax = c(5),y_position = c(0.32), annotation =c('**'), color = "black",tip_length = c(0,0), textsize = textsize)+
  scale_y_continuous(labels = scales::percent)+
  theme(legend.position = 'none')

# 1_3
cells_total_1_3 = subset(cells_total, sub.cluster == paste0('1_',3) & Status !='NRM')
summary_1_3 = make_summary(cells_total_1_3)

prop_1_3 = ggplot(cells_total_1_3, aes(x = Phase, y = n_cells_in/n_cells))+
  geom_bar(data = summary_1_3,
           aes(x = Phase, y = mean_prop, fill = Phase),
           stat = 'identity', inherit.aes = FALSE)+
  geom_errorbar(data = summary_1_3,
                aes(x = Phase, ymin = mean_prop - se, ymax = mean_prop + se),
                width = 0.2, inherit.aes = FALSE)+
  geom_point(size = 0.5)+
  theme_classic()+
  scale_x_discrete(drop = FALSE)+
  scale_y_continuous(labels = scales::percent)+
  labs(x = 'Phase', y = '% of 1_RGC')+
  geom_signif(xmin = c(1),xmax = c(1.9),y_position = c(0.25), annotation =c('**'), color = "black",tip_length = c(0,0), textsize = textsize)+
  geom_signif(xmin = c(1),xmax = c(5),y_position = c(0.27), annotation =c('p=0.524'), color = "black",tip_length = c(0,0), textsize = textsize)+
  geom_signif(xmin = c(2.1),xmax = c(5),y_position = c(0.25), annotation =c('*'), color = "black",tip_length = c(0,0), textsize = textsize)+
  scale_y_continuous(labels = scales::percent)+
  theme(legend.position = 'none')


for(i in 0:3){
  p = get(paste0('prop_1_',i))
  
  ggsave(plot = p,
       file = paste0('rgc_prop_1_',i,'.svg'),
       device = "svg",
       units = "in",
       width = 2,
       height = 2,
       path = "Manuscript/Plots/Manuscript v1.1/Fig.3")

  # no pairwise differences in 1_1 thats why i didnt include it
  
}
