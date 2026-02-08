library(tidyverse)
library(ggplot2)

data = read_csv("/Users/ggraham/Desktop/multiome_poa/Sex Classifier and Linearity/progress_classifier/logistic_individual_probabilities_07_05_2025.csv")

ggplot(data, aes(x = as.factor(cluster), y = prediction, color = status, shape = status))+
  geom_point(size = 4)


mean_dominants = data%>%
  subset(status == 'D')%>%
  group_by(cluster)%>%
  summarize(mean_dom =mean(prediction),
            n_degs = mean(n_degs))

data_with_dom_mean = data%>%
  dplyr::select(-c(n_degs))%>%
  right_join(mean_dominants, by = 'cluster')

data_grouped = data%>%
  group_by(cluster, status)%>%
  summarize(mean = mean(prediction),
            se = sd(prediction)/sqrt(n()),
            mean_dominant = mean(prediction[status =='D']))%>%
  fill(mean_dominant, .direction = 'down')%>%
  right_join(mean_dominants, by = 'cluster')

status =data%>%
  group_by( status)%>%
  summarize(mean = mean(prediction),
            se = sd(prediction)/sqrt(n()))

cluster= data%>%
  group_by( cluster)%>%
  summarize(mean = mean(prediction),
            se = sd(prediction)/sqrt(n()))

ggplot(data_with_dom_mean, aes(x = fct_reorder(as.factor(cluster), mean_dom), y = prediction, color = status, shape = status))+
  geom_point(aes(size = 4), alpha = 0.5)+
  theme_minimal()+
  labs(x = 'Cluster', y = 'Prediction', title = 'Individuals', shape = 'Status', color = 'Status')

ggplot(data_with_dom_mean, aes(x = fct_reorder(as.factor(cluster), mean_dom), y = prediction, color = status, shape = status))+
  stat_summary(geom= 'point', fun='mean', size=3)+
  theme_minimal()+
  labs(x = 'Cluster', y = 'Prediction', title = 'Individuals', shape = 'Status', color = 'Status')


ggplot(data_grouped, aes(x = fct_reorder(as.factor(cluster), mean_dominant), y = mean, color = status, shape = status))+
  geom_point(aes(size = 4), position = position_jitterdodge(0.1,seed =1), alpha = 0.5)+
  geom_errorbar(aes(x = fct_reorder(as.factor(cluster), mean), y = mean, ymin = mean-se, ymax =mean+se), position = position_jitterdodge(0.1, seed = 1), alpha = 0.5)+
    theme_minimal()+
    labs(x = 'Cluster', y = 'Mean Prediction +/- SE', title = 'Groups', shape = 'Status', color = 'Status')

data_subset = subset(data_with_dom_mean, n_degs >15)

ggplot(data_subset, aes(x = fct_reorder(as.factor(cluster), mean_dom), y = prediction, color = status, shape = status))+
  geom_point(aes(size = 4), alpha = 0.5)+
  theme_minimal()+
  labs(x = 'Cluster', y = 'Prediction', title = 'Individuals - Clusters with many DEGs', shape = 'Status', color = 'Status')

data_grouped_subset =subset(data_grouped, n_degs > 15)

ggplot(data_grouped_subset, aes(x = fct_reorder(as.factor(cluster), mean_dominant), y = mean, color = status, shape = status))+
  geom_point(aes(size = 4), position = position_jitterdodge(0.1,seed =1), alpha = 0.5)+
  geom_errorbar(aes(x = fct_reorder(as.factor(cluster), mean), y = mean, ymin = mean-se, ymax =mean+se), position = position_jitterdodge(0.1, seed = 1), alpha = 0.5)+
    theme_minimal()+
    labs(x = 'Cluster', y = 'Mean Prediction +/- SE', title = 'Groups - Clusters with many DEGs', shape = 'Status', color = 'Status')


mf_data = read.csv('/Users/ggraham/Desktop/multiome_poa/Sex Classifier and Linearity/progress_classifier/logistic_mf_validation_07_07_2025.csv')


correct=mf_data_grouped = mf_data%>%
  group_by(cluster, status)%>%
  summarize(mean_proba = mean(proba),
            se = sd(proba)/sqrt(n()))

ggplot(mf_data_grouped, aes(x = as.factor(cluster), y = mean_proba,shape = status, color = status))+
  geom_point()+
  geom_errorbar(aes(x = as.factor(cluster), y = mean_proba, ymin = mean_proba-se, ymax = mean_proba+se))+
  labs(x = 'Cluster', y = 'Mean +/- SE Probability', shape = 'Sex', color = "Sex")+
  theme_minimal()


ggplot(subset(mf_data_grouped, cluster %in% c(1,6,11,0)), aes(x = as.factor(cluster), y = mean_proba,shape = status, color = status))+
  geom_point()+
  geom_errorbar(aes(x = as.factor(cluster), y = mean_proba, ymin = mean_proba-se, ymax = mean_proba+se))+
  labs(x = 'Cluster', y = 'Mean +/- SE Probability', shape = 'Sex', color = "Sex")+
  theme_minimal()

ggplot(subset(mf_data, cluster %in% c(1,6,11,0)), aes(x = as.factor(cluster), y = proba,shape = status, color = status))+
  geom_point()+
  labs(x = 'Cluster', y = 'Mean +/- SE Probability', shape = 'Sex', color = "Sex")+
  theme_minimal()

## I think I can use this data as long as we have good qc, so m < 0.1, f > 0.9

ggplot(subset(mf_data_grouped, mean_proba>0.85 | mean_proba<0.15), aes(x = as.factor(cluster), y = mean_proba,shape = status, color = status))+
  geom_point()+
  geom_errorbar(aes(x = as.factor(cluster), y = mean_proba, ymin = mean_proba-se, ymax = mean_proba+se))+
  labs(x = 'Cluster', y = 'Mean +/- SE Probability', shape = 'Sex', color = "Sex")+
  theme_minimal()

# maybe not?
# maybe I do a 1 sample t test to see if the predicted males are different to 0
# and f different to 1

data =mf_data%>%subset(cluster == 6 & status == 'm')

t.test(data$proba) # not different
data2 =mf_data%>%subset(cluster == 6 & status == 'f')

t.test(1-data2$proba) # not different

#ok maybe this is the way to do it

t_data = data.frame()
for(i in unique(mf_data$cluster)){
  
  data =mf_data%>%subset(cluster == i & status == 'm')

male_p = t.test(data$proba)$p.value
data2 =mf_data%>%subset(cluster == i & status == 'f')

female_p=t.test(1-data2$proba)$p.value 
newd = data.frame(cluster =i, 
                  male_p = male_p,
                  female_p = female_p)
t_data = rbind(t_data, newd)

}
t_data$different = ifelse(t_data$male_p<0.05| t_data$female_p<0.05, '*', NA)
confident_clusters = t_data$cluster[is.na(t_data$different)]

ggplot(subset(data_with_dom_mean, cluster %in% confident_clusters), aes(x = fct_reorder(as.factor(cluster), mean_dom), y = prediction, color = status, shape = status))+
  stat_summary(geom= 'point', fun='mean', size=3)+
  theme_minimal()+
  labs(x = 'Cluster', y = 'Prediction', title = 'Individuals', shape = 'Status', color = 'Status')

plott =subset(data_with_dom_mean, cluster %in% confident_clusters)%>%
  group_by(status, cluster)%>%
  summarize(mean_prediction = mean(prediction),
            se_prediction = sd(prediction)/sqrt(n()))%>%
  subset(status == 'D')

library(forcats)
ggplot(plott, aes(x = fct_reorder(as.factor(cluster), mean_prediction, .desc= T), y = mean_prediction))+
  geom_point()+
  geom_errorbar(aes(x =fct_reorder(as.factor(cluster), mean_prediction,.desc= T), y =mean_prediction ,ymin =mean_prediction -se_prediction,ymax = mean_prediction+se_prediction ), width = 0.4)+
  labs(y = 'Prediction', x= 'Cluster')+
  theme_minimal()

