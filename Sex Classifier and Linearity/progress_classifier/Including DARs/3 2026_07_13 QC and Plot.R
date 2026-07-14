library(Seurat)
library(patchwork)
library(tidyverse)
library(ggplot2)
library(Polychrome)
library(emmeans)
library(ggsignif)
  clown_go = readRDS("Functions/clown_go2")  
library(clusterProfiler)
  library(forcats)



##### Classifier #####
{
mf_data = read.csv("Sex Classifier and Linearity/progress_classifier/Including DARs/logistic_mf_validation_07_13_2026.csv")
data = read_csv("Sex Classifier and Linearity/progress_classifier/Including DARs/logistic_individual_probabilities_07_13_2026.csv")
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


classifier = ggplot(plott, aes(x = fct_reorder(as.factor(cluster), mean_prediction, .desc= T), y = mean_prediction))+
  geom_errorbar(aes(x =fct_reorder(as.factor(cluster), mean_prediction,.desc= T), y =mean_prediction ,ymin =mean_prediction -se_prediction,ymax = mean_prediction+se_prediction ), width = 0.4)+
    geom_point(aes( color = neurotransmitter, shape = neurotransmitter), size = 2)+
  labs(y = 'Prediction', x= 'Cluster')+
  theme_minimal()+
  theme(legend.position = 'none')+
  scale_color_manual(values = c('purple', 'black', '#5bb450', 'maroon'))+
  scale_shape_manual(values =c(15, 16, 17, 18))
}
classifier

ggsave(plot = classifier,
       file = "classifier with dars_c_10.svg",
       device = "svg",
       units = "in",
       width = 2.5,
       height = 1.5,
       path = "Manuscript/Plots/Manuscript v1.3/")

data_with_dom_mean%>%
  summarize(mean_pred = mean(prediction),
            se_pred = sd(prediction)/sqrt(n()))

data_with_dom_mean$neurotransmitter = ifelse(data_with_dom_mean$cluster %in% c(24,0, 6, 9), 'Mixed',NA)
data_with_dom_mean$neurotransmitter = ifelse(data_with_dom_mean$cluster %in% c(1,2, 11, 13, 15, 20, 26), 'Glial',data_with_dom_mean$neurotransmitter)
data_with_dom_mean$neurotransmitter = ifelse(data_with_dom_mean$cluster %in% c(19, 3, 14, 25), 'GABAergic',data_with_dom_mean$neurotransmitter)
data_with_dom_mean$neurotransmitter = ifelse(is.na(data_with_dom_mean$neurotransmitter), 'Glutamatergic', data_with_dom_mean$neurotransmitter)

data_with_dom_mean%>%
  group_by(neurotransmitter)%>%
  summarize(mean_pred = mean(prediction),
            se_pred = sd(prediction)/sqrt(n()))

write.csv(data_with_dom_mean, '/Users/ggraham/Desktop/multiome_poa/Manuscript/Supplementary Tables/classifier_scores.csv')
