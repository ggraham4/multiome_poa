library(lme4)
library(emmeans)
library(ggsignif)
library(car)

cyto_data <- read.csv('/Volumes/jrhodes/Fish Lab/Experiments/sex change single nuc POA/Seurat Outputs/011125 all neurons cyto.csv')

mean_cyto = mean(cyto_data$cyto)

cyto_cluster_signif <- data.frame()
for(i in unique(cyto_data$cluster)){
  if(i ==30 | i == 15){next}
  print(i)
  
  model_data <- subset(cyto_data, cluster ==i)
  
  model = lmer((cyto - mean_cyto)~1 + (1|individual), data = model_data)
  sum =summary(model)$coefficients%>%as.data.frame()
  
  av = car::Anova(model, type = 'III')%>% as.data.frame()
  
  p.value = av$`Pr(>Chisq)`
  
  newd = data.frame(cluster = i,
                    singular = isSingular(model),
                    estimate = sum$Estimate,
                    p.value = p.value)
  
  cyto_cluster_signif = rbind(cyto_cluster_signif, newd)
}
cyto_cluster_signif$q.value = p.adjust(cyto_cluster_signif$p.value, 'fdr', nrow(cyto_cluster_signif))

cyto_cluster_signif$issignif = ifelse(cyto_cluster_signif$q.value<0.05, '*', NA)

sum(cyto_cluster_signif$q.value<0.05 & cyto_cluster_signif$estimate>0)/sum(cyto_cluster_signif$q.value<0.05)
sum(cyto_cluster_signif$q.value<0.05 & cyto_cluster_signif$estimate<0)/sum(cyto_cluster_signif$q.value<0.05)


ggplot(cyto_data, aes(x = as.factor(cluster), y = cyto, group = cluster))+
  geom_boxplot()


sex_differences = data.frame()
for(i in unique(cyto_data$cluster)){
  if(i ==30 | i == 15){next}
  print(i)
  
  model_data <- subset(cyto_data, cluster ==i & status%in%c('M','D','F'))
  
  model = lmer(cyto ~status + (1|individual), data = model_data)
  
  pairs <- emmeans(model, 'status')%>%
    pairs(adjust = 'none')%>%
    as.data.frame()
  
  av = car::Anova(model, type = 'III')%>% as.data.frame()
  
  p.value = av$`Pr(>Chisq)`[2]
  
  newd = data.frame(cluster = i,
                    singular = isSingular(model),
                    d_f_p.value = pairs$p.value[pairs$contrast== 'D - F'],
                    d_m_p.value = pairs$p.value[pairs$contrast== 'D - M'],
                    f_m_p.value = pairs$p.value[pairs$contrast== 'F - M'],
                    p.value = p.value)
  
  sex_differences = rbind(sex_differences, newd)
}

sex_differences$d_f_q.value = p.adjust(sex_differences$d_f_p.value, 'fdr', nrow(sex_differences))
sex_differences$d_m_q.value = p.adjust(sex_differences$d_m_p.value, 'fdr', nrow(sex_differences))
sex_differences$f_m_q.value = p.adjust(sex_differences$f_m_p.value, 'fdr', nrow(sex_differences))

sex_differences$issignif = ifelse(sex_differences$d_f_q.value<0.05|
                                    sex_differences$d_m_q.value<0.05|
                                    sex_differences$f_m_q.value<0.05,
                                  '*',
                                  NA)
#nothing



