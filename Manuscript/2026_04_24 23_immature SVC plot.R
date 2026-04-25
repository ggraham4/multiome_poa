library(ggplot2)
library(tidyverse)
library(ggpattern)  

#cmat_observed_all = data.frame()
#cmat_null_all = data.frame()
predict_observed_all = data.frame()
predict_null_all = data.frame()
for(i in c('male', "intermediate",'female')){

  #cmat_observed = read.csv(paste0("~/Desktop/multiome_poa/Python/support_vector_classifier/2026_03_09_Run_9_phase/train_clusters_predict_new_neurons_df_cmat_100_restarts_observed_",i,".csv"))
  #cmat_null = read.csv(paste0("~/Desktop/multiome_poa/Python/support_vector_classifier/2026_03_09_Run_9_phase/train_clusters_predict_new_neurons_df_cmat_100_restarts_shuffled_",i,".csv"))
  base = '/Users/ggraham/Desktop/multiome_poa/Python/support_vector_classifier/2026_04_22 other immature populations/Clusters/23'
  predict_observed = read.csv(paste0(base,'/train_clusters_predict_new_neurons_pred_norm_100_restarts_observed_',i,'.csv'))
  predict_null = read.csv(paste0(base,"/train_clusters_predict_new_neurons_pred_norm_100_restarts_shuffled_",i,".csv"))
  
  emp_p_vals <- colMeans(predict_null >= predict_observed)
  #names(emp_p_vals) = str_remove(names(emp_p_vals), "X")
  emp_p_vals = as.data.frame(emp_p_vals)
  emp_p_vals$cluster = rownames(emp_p_vals)
    
  pivot_observed=predict_observed%>%
    select(-X)%>%
    pivot_longer(starts_with('X'), values_to = 'percent', names_to = 'cluster')%>%
    group_by(cluster)%>%
    summarize(mean_perc = mean(percent),
              sd_perc = sd(percent))%>%
    mutate(Phase =i)%>%
    right_join(emp_p_vals, by = 'cluster') %>%
  mutate(sig = ifelse(emp_p_vals < 0.05, "*", ""))
  
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

predict_observed_all_null$Phase= factor(predict_observed_all_null$Phase, levels = c('male','intermediate','female'))

predict_observed_all_null$cluster = gsub('X','', predict_observed_all_null$cluster)

predict_observed_all_null$mean_null = predict_observed_all_null$mean_perc.y
predict_observed_all_null$sd_null = predict_observed_all_null$sd_perc.y



p = ggplot(predict_observed_all_null, aes(x = fct_reorder(cluster, mean_perc.x, .desc = T), y = mean_perc.x, color=Phase,
                                      shape = Phase))+
  geom_point(position = position_dodge(width = .9),size =2,
           aes(y = mean_perc.y), color = 'black')+
    geom_errorbar(aes(x = cluster, y=mean_perc.y, ymin = mean_perc.y-sd_perc.y, ymax = mean_perc.y+sd_perc.y), width = 0, position = position_dodge(width = .9), color= 'red')+
geom_point(position = position_dodge(width = .9), size = 2)+
  geom_errorbar(aes(x = cluster, y=mean_perc.x, ymin = mean_perc.x-sd_perc.x, ymax = mean_perc.x+sd_perc.x), position = 'dodge')+
    labs(x = 'Cluster', y = 'Mean Percent +/- SD')+
  theme_minimal()+
    theme(legend.position = 'none')
p # maybe a ribbon for the null?


  predict_observed_male = read.csv(paste0(base,"/train_clusters_predict_new_neurons_pred_norm_100_restarts_observed_",'male',".csv"))
  predict_observed_female = read.csv(paste0(base,"/train_clusters_predict_new_neurons_pred_norm_100_restarts_observed_",'female',".csv"))
  predict_observed_intermediate = read.csv(paste0(base,"/train_clusters_predict_new_neurons_pred_norm_100_restarts_observed_",'intermediate',".csv"))



### claude prettier plot
predict_observed_all$Phase <- factor(predict_observed_all$Phase, levels = c('male','intermediate','female'))
predict_null_all$Phase     <- factor(predict_null_all$Phase,     levels = c('male','intermediate','female'))

predict_observed_all$cluster <- gsub('X','', predict_observed_all$cluster)
predict_null_all$cluster     <- gsub('X','', predict_null_all$cluster)

# Order clusters by observed mean (Male, for consistency)
cluster_order <- predict_observed_all %>%
  filter(Phase == 'male') %>%
  arrange(desc(mean_perc)) %>%
  pull(cluster)%>%unique()

predict_observed_all$cluster <- factor(predict_observed_all$cluster, levels = cluster_order)
predict_null_all$cluster     <- factor(predict_null_all$cluster,     levels = cluster_order)

# Colorblind-safe palette (Okabe-Ito)
predict_observed_all$Phase = ifelse(predict_observed_all$Phase == 'male', 'M', 
                                    ifelse(predict_observed_all$Phase == 'female', 'F',
                                           ifelse(predict_observed_all$Phase == 'intermediate', 'I', NA)
                                    )
)
phase_patterns <- c('M' = 'stripe', 'I' = 'circle', 'F' = 'crosshatch')

col = c('#1965B0', '#4EB265', '#DC050C')

predict_observed_all$Phase = factor(predict_observed_all$Phase, levels = 
                                      c('M',
                                        'I',
                                        'F'))

p1 <- ggplot(predict_observed_all%>%na.omit()) +
  geom_bar(
    aes(x = cluster, y = mean_perc, fill = Phase, group = Phase),
    stat = 'identity',
    position = position_dodge(width = 0.9)
  ) +
  geom_errorbar(
    aes(x = cluster, ymin = mean_perc - sd_perc, ymax = mean_perc + sd_perc, group = Phase),
    position = position_dodge(width = 0.9),
    width = 0.4, linewidth = 1, colour = 'black'
  ) +
  geom_text(
    aes(
      x = cluster,
      y = mean_perc + sd_perc + 2,   # offset above error bar
      label = sig,
      group = Phase
    ),
    position = position_dodge(width = 0.9),
    size = 10
  ) +
  labs(x = 'Cluster', y = 'Mean Percent ± SD', fill = 'Phase') +
  theme_minimal() +
  scale_fill_manual(values = col)+
theme(
    legend.position = c(0.95, 0.95),
    legend.justification = c(1, 1),
    legend.background = element_rect(fill = 'white', colour = 'grey80', linewidth = 0.3),
    legend.margin = margin(4, 4, 4, 4),
    panel.grid.major.x = element_blank()
  )
p1

ggsave(plot = p1,
       file = 'predictions_23_imm.svg',
       device = "svg",
       units = "in",
       width = 4,
       height = 2,
       path = paste0("Manuscript/Plots/Manuscript v1.2/SVC_preds"))


