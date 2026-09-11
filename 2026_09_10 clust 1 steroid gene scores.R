library(ggplot2)
library(tidyverse)
library(ggsignif)

clust_1 = read.csv("Manuscript/cluster_1_steroid_receptor_activity_scores_WITH_SPECIFIC.csv")

clust_1_grouped = clust_1%>%
  group_by(individual, group)%>%
  summarize(mean_esr2b = mean(ESR2B_specific_score),
            se_esr2b = sd(ESR2B_specific_score)/sqrt(n()),
            mean_ar = mean(AR_specific_score),
            se_ar = sd(AR_specific_score)/sqrt(n()),
            mean_pgr = mean(PGR_specific_score),
            se_pgr = sd(PGR_specific_score)/sqrt(n()))

clust_1_grouped$group = factor(clust_1_grouped$group, levels = c('M',
                                                                 'I',
                                                                 'LI',
                                                                 'NF',
                                                                 'F'))
  colors = c('#1965B0', '#4EB265', '#F7F056', '#7BAFDE', '#DC050C')

esr2b =ggplot(clust_1_grouped, aes(x = group, 
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
                        annotation = 'p - 0.924',
                        color      = 'black', tip_length = c(0, 0),
                        textsize   = 3)+

  geom_signif(xmin = 1, xmax = 1.9,
                        y_position = 1,
                        annotation = 'p = 0.297',
                        color      = 'black', tip_length = c(0, 0),
                        textsize   = 3)+
    geom_signif(xmin = 2.1, xmax = 5,
                        y_position = 1,
                        annotation = 'p = 0.498',
                        color      = 'black', tip_length = c(0, 0),
                        textsize   = 3)
  

"
I-F  0.4601519 -0.3081420 1.2284459 0.2971102

M-F  0.1176685 -0.6796278 0.9149648 0.9235498

M-I -0.3424834 -1.1107774 0.4258105 0.4984857
"

ggsave(plot = esr2b,
       file = "esr2b_specific_targets_1.svg",
       device = "svg",
       units = "in",
       width = 2,
       height = 2,
       path = '/Users/ggraham/Desktop/multiome_poa/Manuscript/Plots/')


ar =ggplot(clust_1_grouped, aes(x = group, 
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
                        y_position = 1.75,
                        annotation = 'p = 0.495',
                        color      = 'black', tip_length = c(0, 0),
                        textsize   = 3)+

  geom_signif(xmin = 1, xmax = 1.9,
                        y_position = 1.5,
                        annotation = 'p = 0.998',
                        color      = 'black', tip_length = c(0, 0),
                        textsize   = 3)+
    geom_signif(xmin = 2.1, xmax = 5,
                        y_position = 1.5,
                        annotation = 'p = 0.502',
                        color      = 'black', tip_length = c(0, 0),
                        textsize   = 3)

ar
"
I-F -0.35923292 -1.1689728 0.4505069 0.5016598

M-F -0.37681384 -1.2171206 0.4634929 0.4945459

M-I -0.01758092 -0.8273208 0.7921589 0.9982713
"

ggsave(plot = ar,
       file = "ar_specific_targets_1.svg",
       device = "svg",
       units = "in",
       width = 2,
       height = 2,
       path = '/Users/ggraham/Desktop/multiome_poa/Manuscript/Plots/')

pgr =ggplot(clust_1_grouped, aes(x = group, 
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
                        y_position = 2.25,
                        annotation = 'p = 0.762',
                        color      = 'black', tip_length = c(0, 0),
                        textsize   = 3)+

  geom_signif(xmin = 1, xmax = 1.9,
                        y_position = 2,
                        annotation = 'p = 0.998',
                        color      = 'black', tip_length = c(0, 0),
                        textsize   = 3)+
    geom_signif(xmin = 2.1, xmax = 5,
                        y_position = 2,
                        annotation = 'p - 0.708',
                        color      = 'black', tip_length = c(0, 0),
                        textsize   = 3)
pgr
"
I-F -0.2020691 -0.8538907 0.4497526 0.7084080

M-F -0.1855536 -0.8619808 0.4908737 0.7624221

M-I  0.0165155 -0.6353061 0.6683371 0.9976465
"

ggsave(plot = pgr,
       file = "pgr_specific_targets_1.svg",
       device = "svg",
       units = "in",
       width = 2,
       height = 2,
       path = '/Users/ggraham/Desktop/multiome_poa/Manuscript/Plots/')


