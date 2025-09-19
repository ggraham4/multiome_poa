retinoic_acid = c('LOC111582092',
                'rdh8a',
                'rdh12',
                'rdh20',
                'rdh1',
                'LOC111572826',                'LOC111587579',
                'LOC111587581','LOC111579690','si:dkey-23o4.6','LOC111575778','LOC111587582','LOC111575523','LOC111566272','LOC111569321','rdh14b','rdh5','rdh10a','LOC111570850')

obj <- AddModuleScore(obj, list(retinoic_acid),
                      name = 'Retinoic')
DotPlot(obj, 'Retinoic1')+
  coord_flip()

retinoic_data = obj@meta.data%>%
  group_by(individual, Status, final_clusters)%>%
  summarize(mean_retinoic = mean(Retinoic1),
            se_retinoic = sd(Retinoic1)/sqrt(n()))

retinoic_data$Status <- factor(retinoic_data$Status, levels= c('NRM','M','D','E','NF','F'))

ggplot(subset(retinoic_data, final_clusters ==8), aes(x = Status, y = mean_retinoic))+
  geom_boxplot()+
    geom_point()

ggplot(subset(retinoic_data, final_clusters ==22), aes(x = Status, y = mean_retinoic))+
  geom_boxplot()+
    geom_point()

ggplot(subset(retinoic_data, final_clusters ==24), aes(x = Status, y = mean_retinoic))+
  geom_boxplot()+
    geom_point()

ggplot(subset(retinoic_data, final_clusters ==23), aes(x = Status, y = mean_retinoic))+
  geom_boxplot()+
    geom_point()

dme_data = data.frame()
for(i in 0:26){
  dat = subset(obj@meta.data, final_clusters ==i & Status %in% c('M','D',"F"))
  mod = lmer(Retinoic1~Status+(1|individual), data = dat)
  av = car::Anova(mod, 3)
  p = av$`Pr(>Chisq)`[2]
  dme_data = rbind(dme_data, data.frame(cluster = i, 
                                        p = p,
                                        singular = isSingular(model)))
  
}

ggplot(subset(retinoic_data, final_clusters ==14), aes(x = Status, y = mean_retinoic))+
  geom_boxplot()+
    geom_point()

