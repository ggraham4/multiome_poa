library(ggplot2)
library(tidyverse)
library(ggsignif)
library(forcats)

data = read_csv("/Users/ggraham/Desktop/multiome_poa/Sex Classifier and Linearity/progress_classifier/logistic_individual_probabilities_06_30_2025.csv")

### first, lets group by sex and cluster and gather means

data_grouped_cluster_sex = data%>%
  group_by(cluster, status)%>%
  summarize(mean_pred = mean(prediction),
            se_pred = sd(prediction)/sqrt(n()))

plot_1 = ggplot(data_grouped_cluster_sex, aes(x = as.factor(cluster), mean_pred, y = mean_pred, shape = status, color =status))+
   geom_errorbar(aes(x = as.factor(cluster), y = mean_pred, ymin = mean_pred-se_pred, ymax = mean_pred+se_pred), alpha = 0.5)+
   geom_point(size = 3)+
  scale_shape_manual(values = c(1,2,3))+
  labs(x ='Cluster', y= 'Mean +/- SE Prediction')
plot_1


ggplot(data_grouped_cluster_sex, aes(x = fct_reorder(as.factor(cluster), mean_pred), mean_pred, y = mean_pred, shape = status, color =status))+
   geom_point(size = 3)+
  scale_shape_manual(values = c(1,2,3))+
  labs(x ='Cluster', y= 'Mean +/- SE Prediction')

data_grouped_sex = data%>%
  group_by( status)%>%
  summarize(mean_pred = mean(prediction),
            se_pred = sd(prediction)/sqrt(n()))

data_grouped_cluster = data%>%
  group_by(cluster)%>%
  summarize(mean_pred = mean(prediction),
            se_pred = sd(prediction)/sqrt(n()))

data_grouped_cluster



