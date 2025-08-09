#write.csv(cyto_data, 'X:/Fish Lab/Experiments/sex change single nuc POA/Seurat Outputs/011125 all neurons cyto.csv')
cyto_data <- read.csv( 'X:/Fish Lab/Experiments/sex change single nuc POA/Seurat Outputs/011125 all neurons cyto.csv')

cyto_data$label <- ifelse(cyto_data$cluster==27, '27', NA)

ggplot(cyto_data,aes(x =as.factor(cluster), y = cyto, fill = label))+
  geom_boxplot()

cyto_data$status <- factor(cyto_data$status, levels = c('NRM','M','D','E','NF','F'))
ggplot(subset(cyto_data, cluster ==27),aes(x =status, y = cyto, fill = status))+
  geom_violin()+
  geom_boxplot(fill = 'white',  width = 0.1, size =1)+
  theme_minimal()

cyto_data_27 <- subset(cyto_data, cluster ==27  & status %in% c('M','D','F'))
model_27 <- lmer(cyto~status+(1|individual), data =cyto_data_27)
car::Anova(model_27, type = 'III')
pairs(emmeans(model_27, 'status'), adjust = 'none')
#no difference

cluster_data <-obj@meta.data%>%
  group_by(Status, individual, harmony.wnn_res0.4_clusters)%>%
  summarize(cells_in_cluster = n())

total_data <-obj@meta.data%>%
  group_by(Status, individual)%>%
  summarize(cells_total = n())

cluster_data_joined <- cluster_data%>%
  right_join(total_data, by ='individual')
  
cluster_data_joined$cells_not_in_cluster <-cluster_data_joined$cells_total- cluster_data_joined$cells_in_cluster
gen_data <- data.frame()
for(i in 0:31){
  print(i)
  cluster_data_joined2 <- subset(cluster_data_joined, harmony.wnn_res0.4_clusters==i)
  
ind_matrix <- matrix(NA, nrow(cluster_data_joined2), 2)
ind_matrix[,1] <- cluster_data_joined2$cells_in_cluster
ind_matrix[,2] <- cluster_data_joined2$cells_not_in_cluster

model <- glmer(ind_matrix~Status.y + (1|individual), data = cluster_data_joined2, family = binomial('logit'))

av <- as.data.frame(car::Anova(model, type = 'III'))

pairs <- as.data.frame(pairs(emmeans(model, 'Status.y'), adjust ='none'))
d_f_p.value = pairs$p.value[pairs$contrast == 'D - F']
d_m_p.value = pairs$p.value[pairs$contrast == 'D - M']
f_m_p.value = pairs$p.value[pairs$contrast == 'F - M']

new_data <- data.frame(cluster = i, 
                       chisq = av$`Pr(>Chisq)`[2],
                       d_f_p.value,
                       d_m_p.value,
                       f_m_p.value
)
gen_data <- rbind(new_data, gen_data)
}

gen_data$issignif = ifelse(gen_data$chisq <0.05 | gen_data$d_f_p.value<0.05 | gen_data$d_m_p.value<0.05| gen_data$f_m_p.value<0.05,
                           '*',
                           NA)

cluster_data_joined$Status.y <- factor(cluster_data_joined$Status.y , levels = c('M','D','F'))
library(ggsignif)
ggplot(data = subset(cluster_data_joined,harmony.wnn_res0.4_clusters==27 & Status.y %in% c('M','D','F')), aes(y= cells_in_cluster/cells_total, x= Status.y, shape = Status.y, group = Status.y))+
  geom_point(position = position_jitterdodge(dodge.width=.9), size = 4, aes(color = Status.y))+
  stat_summary(geom = 'bar', fun = "mean", position = 'dodge', fill = NA, size = 1, color = 'black')+
  stat_summary(geom = 'errorbar', fun.data = mean_se, width = 0.4, position = position_dodge(.9), linewidth = 1)+
  theme_minimal()+
  geom_signif(xmin = 2, xmax = 3, y_position = 0.015, annotation = '*', color = 'black', textsize = 10, tip_length =c(0,0))+
  ylim(0,0.017)+
  labs(x = 'Sex', y = 'Proportion', title = 'Cluster 27', shape = 'Sex', color ='Sex')


cluster_27 <- subset(cluster_data_joined,harmony.wnn_res0.4_clusters==27 & Status.y %in% c('M','D','F'))
                     

cluster_27_markers <- FindMarkers(obj, 27)
