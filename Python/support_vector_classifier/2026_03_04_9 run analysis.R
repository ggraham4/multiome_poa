cmat_observed = read.csv("~/Desktop/multiome_poa/Python/support_vector_classifier/2026_03_04_Run_9/train_clusters_predict_new_neurons_df_cmat_100_restarts_observed.csv")
cmat_null = read.csv("~/Desktop/multiome_poa/Python/support_vector_classifier/2026_03_04_Run_9/train_clusters_predict_new_neurons_df_cmat_100_restarts_shuffled.csv")

predict_observed = read.csv("~/Desktop/multiome_poa/Python/support_vector_classifier/2026_03_04_Run_9/train_clusters_predict_new_neurons_pred_norm_100_restarts_observed.csv")
predict_null = read.csv("~/Desktop/multiome_poa/Python/support_vector_classifier/2026_03_04_Run_9/train_clusters_predict_new_neurons_pred_norm_100_restarts_shuffled.csv")

predict_observed_means = colMeans(predict_observed[,-1]) # remove indexing col
predict_null_means = colMeans(predict_null[-1])

predict_observed_means%>%sort() ##### cluster 0 and cluster 6 interesting
predict_null_means%>%sort()

hist(predict_observed_means)
hist(predict_null_means)

# lets see if I can make this a pretty plot
#Fold change pointrange char

# lets also do some statistics to see which are significantly different from their nulls

p_one_tail = 1-(colMeans(predict_observed[-1] >= predict_null[-1]))
p_two_tail <- 2 * pmin(p_one_tail, 1 - p_one_tail)
signif = ifelse(p_two_tail<0.05, '*', '')

newd = data.frame(observed = predict_observed_means,
                  null = predict_null_means, 
                  observed_sd = apply(FUN='sd', X=predict_observed[,-1], MARGIN = 2),
                  null_sd = apply(FUN='sd', X=predict_null[,-1], MARGIN = 2),
                  cluster = gsub('X','',names(predict_observed_means)),
                  empirical_p = p_two_tail,
                  signif = signif)

newd$cluster = fct_reorder(newd$cluster, newd$observed, .desc= T)


svc_9 = ggplot(newd, aes(x = cluster, y = observed))+
  geom_pointrange(aes(
    x = cluster, y = observed, ymin = observed -observed_sd, ymax = observed+observed_sd, color = 'Observed'),
    position = position_nudge(x = -0.15, y = 0), size =0.2)+
  geom_pointrange(aes(x = cluster, y = null, ymin = null -null_sd, ymax = null+null_sd, color = 'Null'),
                  position = position_nudge(x = 0.15, y = 0), shape = 17, size =0.2)+
  geom_text(aes(x = cluster, y = (observed+observed_sd+1), label = signif), size = 7,position = position_nudge(x = -0.15, y = 0))+
  scale_color_manual(values = c('red', 'blue'))+
  theme_minimal()+
  labs(x = 'Cluster', y = '% of Cells', color = '', title = 'Prediction from 9')+
  theme(legend.position = 'none')
svc_9

ggsave(plot = svc_9,
       file = "svc_9.svg",
       device = "svg",
       units = "in",
       width = 3.5,
       height = 1.5,
       path = "/Users/ggraham/Desktop/multiome_poa/Manuscript/Plots/SVC/")


