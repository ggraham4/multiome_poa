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


ggplot(data_with_dom_mean, aes(x = fct_reorder(as.factor(cluster), mean_dom), y = prediction, color = status, shape = status))+
  geom_point(aes(size = 4), alpha = 0.5)+
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

