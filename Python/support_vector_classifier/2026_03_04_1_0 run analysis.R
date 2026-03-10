cmat_observed = read.csv("~/Desktop/multiome_poa/Python/support_vector_classifier/2026_03_04_Run_1_0/train_clusters_predict_new_neurons_df_cmat_100_restarts_observed.csv")
cmat_null = read.csv("~/Desktop/multiome_poa/Python/support_vector_classifier/2026_03_04_Run_1_0/train_clusters_predict_new_neurons_df_cmat_100_restarts_shuffled.csv")

predict_observed = read.csv("~/Desktop/multiome_poa/Python/support_vector_classifier/2026_03_04_Run_1_0/train_clusters_predict_new_neurons_pred_norm_100_restarts_observed.csv")
predict_null = read.csv("~/Desktop/multiome_poa/Python/support_vector_classifier/2026_03_04_Run_1_0/train_clusters_predict_new_neurons_pred_norm_100_restarts_shuffled.csv")

predict_observed_means = colMeans(predict_observed[,-1]) # remove indexing col
predict_null_means = colMeans(predict_null[-1])

predict_observed_means%>%sort() ##### cluster 9 ding ding ding
predict_null_means%>%sort() ## basically a uniform distribution

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


svc_1_0a = ggplot(newd, aes(x = cluster, y = observed))+
  geom_pointrange(aes(
    x = cluster, y = observed, ymin = observed -observed_sd, ymax = observed+observed_sd, color = 'Observed'),
    position = position_nudge(x = -0.15, y = 0), size =0.2)+
  geom_pointrange(aes(x = cluster, y = null, ymin = null -null_sd, ymax = null+null_sd, color = 'Null'),
                  position = position_nudge(x = 0.15, y = 0), shape = 17, size =0.2)+
  geom_text(aes(x = cluster, y = (observed+observed_sd+1), label = signif), size = 7,position = position_nudge(x = -0.15, y = 0))+
  scale_color_manual(values = c('red', 'blue'))+
  theme_minimal()+
  labs(x = 'Cluster', y = '% of Cells', color = '', title = 'Prediction from 1_0')+
  theme(legend.position = 'none')
svc_1_0a

newd$cluster <- fct_reorder(newd$cluster, newd$observed, .desc = TRUE)

bar_patterns <- c('Observed' = 'stripe', 'Null' = 'none')

svc_1_0 <- ggplot(newd, aes(x = cluster)) +
  geom_bar_pattern(
    aes(y = observed, pattern = 'Observed'),
    stat = 'identity',
    fill = '#4B9CD3',
    colour = 'black',
    pattern_colour = 'black',
    pattern_fill = 'black',
    pattern_density = 0.1,
    alpha = 0.5,
    width = 0.7
  ) +
  geom_bar(
    aes(y = null),
    stat = 'identity',
    fill = 'black',
    colour = 'black',
    alpha = 0.5,
    width = 0.7
  ) +
  geom_errorbar(
    aes(ymin = observed - observed_sd, ymax = observed + observed_sd),
    width = 0.2, linewidth = 0.5, colour = 'black'
  ) +
  geom_errorbar(
    aes(ymin = null - null_sd, ymax = null + null_sd),
    width = 0.2, linewidth = 0.5, colour = 'black'
  ) +
  geom_text(
    aes(y = observed + observed_sd + 1, label = signif),
    size = 5
  ) +
  scale_pattern_manual(values = bar_patterns, labels = c('Observed', 'Null')) +
  labs(x = 'Cluster', y = '% of Cells', pattern = '', title = 'Prediction from 1_0') +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid.major.x = element_blank(),
    legend.position = 'none'
  )

svc_1_0

ggsave(plot = svc_1_0,
       file = "svc_1_0_a.svg",
       device = "svg",
       units = "in",
       width = 3.4,
       height = 1.5,
       path = "/Users/ggraham/Desktop/multiome_poa/Manuscript/Plots/SVC/")


####3 also make this for the trajectory analysis
library(ggplot2)

df <- data.frame(x = 1, y = factor(1:18))

ggplot(df, aes(x = x, y = y, color = y)) +
  geom_tile(aes(fill = y), width = 0.3, height = 0.5) +
  scale_fill_discrete() +
  theme_void() +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    legend.position = "none"
  )

ggsave("/Users/ggraham/Desktop/multiome_poa/Manuscript/Plots/lineage colors.svg", width = 1.15, height = 2.7)
