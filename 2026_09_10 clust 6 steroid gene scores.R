library(ggplot2)
library(tidyverse)
library(ggsignif)

clust_6 = read.csv("Manuscript/cluster_6_steroid_receptor_activity_scores_WITH_SPECIFIC.csv")

clust_6_grouped = clust_6%>%
  group_by(individual, group)%>%
  summarize(mean_esr2b = mean(ESR2B_specific_score),
            se_esr2b = sd(ESR2B_specific_score)/sqrt(n()),
            mean_ar = mean(AR_specific_score),
            se_ar = sd(AR_specific_score)/sqrt(n()),
            mean_pgr = mean(PGR_specific_score),
            se_pgr = sd(PGR_specific_score)/sqrt(n()))

clust_6_grouped$group = factor(clust_6_grouped$group, levels = c('M',
                                                                 'I',
                                                                 'LI',
                                                                 'NF',
                                                                 'F'))
  colors = c('#1965B0', '#4EB265', '#F7F056', '#7BAFDE', '#DC050C')

esr2b =ggplot(clust_6_grouped, aes(x = group, 
                            y = mean_esr2b,
                            fill = group))+
  geom_boxplot(outlier.shape = NA)+
    geom_point(size     = 1,
               position = position_jitter(width = 0.1, seed = 42)) +
    scale_x_discrete(drop = FALSE) +
    scale_fill_manual(values = colors) +
    theme_classic() +
      theme(legend.position = 'none')+
  labs(y = "ESR2B Specific Gene Z-Score",x = 'Phase')+
  geom_signif(xmin = 1, xmax = 5,
                        y_position = 1.25,
                        annotation = '*',
                        color      = 'black', tip_length = c(0, 0),
                        textsize   = 6)+

  geom_signif(xmin = 1, xmax = 1.9,
                        y_position = 1,
                        annotation = 'p = 0.435',
                        color      = 'black', tip_length = c(0, 0),
                        textsize   = 3)+
    geom_signif(xmin = 2.1, xmax = 5,
                        y_position = 1,
                        annotation = 'p = 0.153',
                        color      = 'black', tip_length = c(0, 0),
                        textsize   = 3)
  

"
I-F - 0.1534015

M-F 0.0174705

M-I  0.4347995
"

ggsave(plot = esr2b,
       file = "esr2b_specific_targets_6.svg",
       device = "svg",
       units = "in",
       width = 2,
       height = 2,
       path = '/Users/ggraham/Desktop/multiome_poa/Manuscript/Plots/')


ar =ggplot(clust_6_grouped, aes(x = group, 
                            y = mean_ar,
                            fill = group))+
  geom_boxplot(outlier.shape = NA)+
    geom_point(size     = 1,
               position = position_jitter(width = 0.1, seed = 42)) +
    scale_x_discrete(drop = FALSE) +
    scale_fill_manual(values = colors) +
    theme_classic() +
      theme(legend.position = 'none')+
  labs(y = "AR Specific Gene Z-Score",x = 'Phase')+
  geom_signif(xmin = 1, xmax = 5,
                        y_position = 1.25,
                        annotation = '*',
                        color      = 'black', tip_length = c(0, 0),
                        textsize   = 6)+

  geom_signif(xmin = 1, xmax = 1.9,
                        y_position = 1,
                        annotation = 'p = 0.085',
                        color      = 'black', tip_length = c(0, 0),
                        textsize   = 3)+
    geom_signif(xmin = 2.1, xmax = 5,
                        y_position = 1,
                        annotation = 'p = 0.590',
                        color      = 'black', tip_length = c(0, 0),
                        textsize   = 3)

ar
"
I-F  0.08504682

M-F - 0.01534851

M-I  0.58965971
"

ggsave(plot = ar,
       file = "ar_specific_targets_6.svg",
       device = "svg",
       units = "in",
       width = 2,
       height = 2,
       path = '/Users/ggraham/Desktop/multiome_poa/Manuscript/Plots/')

pgr =ggplot(clust_6_grouped, aes(x = group, 
                            y = mean_pgr,
                            fill = group))+
  geom_boxplot(outlier.shape = NA)+
    geom_point(size     = 1,
               position = position_jitter(width = 0.1, seed = 42)) +
    scale_x_discrete(drop = FALSE) +
    scale_fill_manual(values = colors) +
    theme_classic() +
      theme(legend.position = 'none')+
  labs(y = "PGR Specific Gene Z-Score",x = 'Phase')+
  geom_signif(xmin = 1, xmax = 5,
                        y_position = 1.25,
                        annotation = 'p = 0.722',
                        color      = 'black', tip_length = c(0, 0),
                        textsize   = 3)+

  geom_signif(xmin = 1, xmax = 1.9,
                        y_position = 1,
                        annotation = 'p = 0.114',
                        color      = 'black', tip_length = c(0, 0),
                        textsize   = 3)+
    geom_signif(xmin = 2.1, xmax = 5,
                        y_position = 1,
                        annotation = '*',
                        color      = 'black', tip_length = c(0, 0),
                        textsize   = 6)
pgr
"
I-F 0.02448413

M-F 0.72228253

M-I   0.11362804
"

ggsave(plot = pgr,
       file = "pgr_specific_targets_6.svg",
       device = "svg",
       units = "in",
       width = 2,
       height = 2,
       path = '/Users/ggraham/Desktop/multiome_poa/Manuscript/Plots/')


