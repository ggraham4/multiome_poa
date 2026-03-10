library(ggplot2)
library(tidyverse)
library(ggpattern)  

#cmat_observed_all = data.frame()
#cmat_null_all = data.frame()
predict_observed_all = data.frame()
predict_null_all = data.frame()
for(i in c('Male', "Intermediate",'Female')){

  #cmat_observed = read.csv(paste0("~/Desktop/multiome_poa/Python/support_vector_classifier/2026_03_09_Run_9_phase/train_clusters_predict_new_neurons_df_cmat_100_restarts_observed_",i,".csv"))
  #cmat_null = read.csv(paste0("~/Desktop/multiome_poa/Python/support_vector_classifier/2026_03_09_Run_9_phase/train_clusters_predict_new_neurons_df_cmat_100_restarts_shuffled_",i,".csv"))

  predict_observed = read.csv(paste0("~/Desktop/multiome_poa/Python/support_vector_classifier/2026_03_09_Run_9_phase/train_clusters_predict_new_neurons_pred_norm_100_restarts_observed_",i,".csv"))
  predict_null = read.csv(paste0("~/Desktop/multiome_poa/Python/support_vector_classifier/2026_03_09_Run_9_phase/train_clusters_predict_new_neurons_pred_norm_100_restarts_shuffled_",i,".csv"))
  
  pivot_observed=predict_observed%>%
    select(-X)%>%
    pivot_longer(starts_with('X'), values_to = 'percent', names_to = 'cluster')%>%
    group_by(cluster)%>%
    summarize(mean_perc = mean(percent),
              sd_perc = sd(percent))%>%
    mutate(Phase =i)
  
    pivot_null=predict_null%>%
    select(-X)%>%
    pivot_longer(starts_with('X'), values_to = 'percent', names_to = 'cluster')%>%
    group_by(cluster)%>%
    summarize(mean_perc = mean(percent),
              sd_perc = sd(percent))%>%
    mutate(Phase =i)
    
    predict_observed_all = rbind(predict_observed_all, pivot_observed)
      predict_null_all = rbind(predict_null_all, pivot_null)

  
  
}

predict_observed_all_null = predict_observed_all%>%
  right_join(predict_null_all, by = c('cluster', 'Phase'))

predict_observed_all_null$Phase= factor(predict_observed_all_null$Phase, levels = c('Male','Intermediate','Female'))

predict_observed_all_null$cluster = gsub('X','', predict_observed_all_null$cluster)

p = ggplot(predict_observed_all_null, aes(x = fct_reorder(cluster, mean_perc.x, .desc = T), y = mean_perc.x, color=Phase,
                                      shape = Phase))+
  geom_point(position = position_dodge(width = .9),size =2,
           aes(y = mean_perc.y), color = 'black')+
    geom_errorbar(aes(x = cluster, y=mean_perc.y, ymin = mean_perc.y-sd_perc.y, ymax = mean_perc.y+sd_perc.y), width = 0, position = position_dodge(width = .9), color= 'black')+
geom_point(position = position_dodge(width = .9), size = 2)+
  geom_errorbar(aes(x = cluster, y=mean_perc.x, ymin = mean_perc.x-sd_perc.x, ymax = mean_perc.x+sd_perc.x), position = 'dodge')+
    labs(x = 'Cluster', y = 'Mean Percent +/- SD')+
  theme_minimal()+
    theme(legend.position = 'none')
p # maybe a ribbon for the null?


  predict_observed_male = read.csv(paste0("~/Desktop/multiome_poa/Python/support_vector_classifier/2026_03_09_Run_9_phase/train_clusters_predict_new_neurons_pred_norm_100_restarts_observed_",'male',".csv"))
  predict_observed_female = read.csv(paste0("~/Desktop/multiome_poa/Python/support_vector_classifier/2026_03_09_Run_9_phase/train_clusters_predict_new_neurons_pred_norm_100_restarts_observed_",'female',".csv"))
  predict_observed_intermediate = read.csv(paste0("~/Desktop/multiome_poa/Python/support_vector_classifier/2026_03_09_Run_9_phase/train_clusters_predict_new_neurons_pred_norm_100_restarts_observed_",'intermediate',".csv"))


# have to do kruskal wallis test becasue the data are not perfectly IID
clust_6_data = data.frame(clust_6 = c(predict_observed_male$X6, predict_observed_female$X6, predict_observed_intermediate$X6),
                     phase =c(rep('M', 100),rep('F', 100), rep('I', 100)) )
test = kruskal.test(clust_6 ~ phase, data = clust_6_data)
test
pairwise.wilcox.test(clust_6_data$clust_6, clust_6_data$phase)

clust_0_data = data.frame(clust_0 = c(predict_observed_male$X0, predict_observed_female$X0, predict_observed_intermediate$X0),
                     phase =c(rep('M', 100),rep('F', 100), rep('I', 100)) )

test2 = kruskal.test(clust_0 ~ phase, data = clust_0_data)
test2
pairwise.wilcox.test(clust_0_data$clust_0, clust_0_data$phase)

# it is very significant but it would probably be more precise to do it by 
# individual? talk to justin
# the individuals do not have enough cells for the bootstrap approach, so this
# is probably the better way to do it

### claude prettier plot
predict_observed_all$Phase <- factor(predict_observed_all$Phase, levels = c('Male','Intermediate','Female'))
predict_null_all$Phase     <- factor(predict_null_all$Phase,     levels = c('Male','Intermediate','Female'))

predict_observed_all$cluster <- gsub('X','', predict_observed_all$cluster)
predict_null_all$cluster     <- gsub('X','', predict_null_all$cluster)

# Order clusters by observed mean (Male, for consistency)
cluster_order <- predict_observed_all %>%
  filter(Phase == 'Male') %>%
  arrange(desc(mean_perc)) %>%
  pull(cluster)%>%unique()

predict_observed_all$cluster <- factor(predict_observed_all$cluster, levels = cluster_order)
predict_null_all$cluster     <- factor(predict_null_all$cluster,     levels = cluster_order)

# Colorblind-safe palette (Okabe-Ito)
phase_patterns <- c('Male' = 'stripe', 'Intermediate' = 'circle', 'Female' = 'crosshatch')

p1 <- ggplot() +
  geom_bar_pattern(
    data = predict_observed_all,
    aes(x = cluster, y = mean_perc, fill = Phase, group = Phase, pattern = Phase),
    stat = 'identity',
    position = position_dodge(width = 0.9),
    colour = 'black',
    pattern_colour = 'black',
    pattern_fill = 'black',
    pattern_density = 0.3,
    alpha = 0.5,
    pattern_alpha = 0.5,
    width = 0.7
  ) +
  geom_bar(
    data = predict_null_all,
    aes(x = cluster, y = mean_perc, group = Phase),
    stat = 'identity',
    position = position_dodge(width = 0.9),
    fill = 'black',
    colour = 'black',
    alpha = 0.5,
    width = 0.7
  ) +
  geom_errorbar(
    data = predict_observed_all,
    aes(x = cluster, ymin = mean_perc - sd_perc, ymax = mean_perc + sd_perc, group = Phase),
    position = position_dodge(width = 0.9),
    width = 0.2, linewidth = 0.5, colour = 'black'
  ) +
  geom_errorbar(
    data = predict_null_all,
    aes(x = cluster, ymin = mean_perc - sd_perc, ymax = mean_perc + sd_perc, group = Phase),
    position = position_dodge(width = 0.9),
    width = 0.2, linewidth = 0.5, colour = 'black'
  ) +
  scale_pattern_manual(values = phase_patterns) +
  labs(x = 'Cluster', y = 'Mean Percent ± SD', fill = 'Phase', pattern = 'Phase') +
  theme_minimal(base_size = 10) +
  theme(
    legend.position = c(0.95, 0.95),
    legend.justification = c(1, 1),
    legend.background = element_rect(fill = 'white', colour = 'grey80', linewidth = 0.3),
    legend.margin = margin(4, 4, 4, 4),
    panel.grid.major.x = element_blank()
  )

p1


ggsave(plot = p1,
       file = "svc_9_phase.svg",
       device = "svg",
       units = "in",
       width = 3.4,
       height = 3.4,
       path = "/Users/ggraham/Desktop/multiome_poa/Manuscript/Plots/SVC/")


