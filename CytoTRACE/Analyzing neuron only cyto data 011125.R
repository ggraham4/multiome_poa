###analyzing neuron only cyto data
library(lme4)
library(emmeans)
library(ggsignif)

cyto_data <- read.csv('X:/Fish Lab/Experiments/sex change single nuc POA/Seurat Outputs/011125 all neurons cyto.csv')

cyto_model <- lmer(cyto~cluster+(1|individual), data = cyto_data)
summary(cyto_model)
#of course they significantly differ
car::Anova(cyto_model, type = 'III')

cyto_model_status <- lmer(cyto~as.factor(cluster)*status+(1|individual), data = cyto_data)
car::Anova(cyto_model_status, type = 'III')

cyto_model_output <- as.data.frame(pairs(emmeans(cyto_model_status, 'status', by ='cluster'), adjust ='none'))
cyto_model_output$q.value <- p.adjust(cyto_model_output$p.value, 'fdr', nrow(cyto_model_output))

cyto_model_output_signif <- subset(cyto_model_output, q.value<0.05)

#differences in cluster 15 and 31

cyto_data$status <- factor(cyto_data$status, levels = c('NRM','M','D','E','NF','F'))
ggplot(subset(cyto_data, cluster ==15), aes(x = status, y = cyto, fill = status))+
  geom_violin()+
  geom_boxplot(color= 'white', fill = NA, width = .3, linewidth =1.5)+
  geom_point(position = position_jitterdodge(0.5))+
  geom_signif(xmin = 2, xmax =3, y_position = 1.05, annotation = c('**'), color = "black", textsize = 10, tip_length = c(0,0))+
  geom_signif(xmin = 1, xmax =3, y_position = 1.15, annotation = c('*'), color = "black", textsize = 10, tip_length = c(0,0))+
  ylim(0,1.2)+
  labs(title ='15')+
  theme(plot.title = element_text(hjust = 0.5))


ggplot(subset(cyto_data, cluster ==31), aes(x = status, y = cyto, fill = status))+
  geom_violin()+
  geom_boxplot(color= 'white', fill = NA, width = .3, linewidth =1.5)+
  geom_point(position = position_jitterdodge(0.5))+
  geom_signif(xmin = 3, xmax =6, y_position = 1.05, annotation = c('**'), color = "black", textsize = 10, tip_length = c(0,0))+
  geom_signif(xmin = 2, xmax =6, y_position = 1.15, annotation = c('**'), color = "black", textsize = 10, tip_length = c(0,0))+
  ylim(0,1.2)+
  labs(title ='31')+
  theme(plot.title = element_text(hjust = 0.5))

  

