library(ggplot2)
library(tidyverse)
library(ggsignif)

clust_6 = read.csv("Manuscript/updatedcluster_6_steroid_receptor_SPECIFIC_PROMOTER_ZSCORES.csv")

clust_6_grouped = clust_6%>%
  group_by(individual, group)%>%
  summarize(mean_esr2b = mean(ESR2B_score),
            se_esr2b = sd(ESR2B_score)/sqrt(n()),
            mean_ar = mean(AR_score),
            se_ar = sd(AR_score)/sqrt(n()),
            mean_pgr = mean(PGR_score),
            se_pgr = sd(PGR_score)/sqrt(n()))

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
                        y_position = 0.25,
                        annotation = '**',
                        color      = 'black', tip_length = c(0, 0),
                        textsize   = 6)+

  geom_signif(xmin = 1, xmax = 1.9,
                        y_position = 0.2,
                        annotation = 'p = 0.256',
                        color      = 'black', tip_length = c(0, 0),
                        textsize   = 3)+
    geom_signif(xmin = 2.1, xmax = 5,
                        y_position = 0.2,
                        annotation = 'p = 0.069',
                        color      = 'black', tip_length = c(0, 0),
                        textsize   = 3)
esr2b

## statistical test
clust_6_grouped_mdf = subset(clust_6_grouped, group %in% c("M",'I',"F"))
clust_6_grouped_mdf$group = factor(clust_6_grouped_mdf$group, levels = c('M', "I","F"))
mod_esr2b = lm(mean_esr2b~group, data = clust_6_grouped_mdf)
anova(mod_esr2b, test = 'Chisq')
pairs(emmeans(mod_esr2b, 'group'), adjust = 'none')

range(clust_6_grouped_mdf$mean_esr2b)

ggsave(plot = esr2b,
       file = "esr2b_targets_1.svg",
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
                        y_position = 0.25,
                        annotation = 'p = 0.100',
                        color      = 'black', tip_length = c(0, 0),
                        textsize   = 3)+

  geom_signif(xmin = 1, xmax = 1.9,
                        y_position = 0.2,
                        annotation = 'p = 0.243',
                        color      = 'black', tip_length = c(0, 0),
                        textsize   = 3)+
    geom_signif(xmin = 2.1, xmax = 5,
                        y_position = 0.2,
                        annotation = 'p = 0.558',
                        color      = 'black', tip_length = c(0, 0),
                        textsize   = 3)
ar

mod_ar = lm(mean_ar~group, data = clust_6_grouped_mdf)
anova(mod_ar, test = 'Chisq')
pairs(emmeans(mod_ar, 'group'), adjust = 'none')
range(clust_6_grouped_mdf$mean_ar)

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
                        y_position = 0.2,
                        annotation = '**',
                        color      = 'black', tip_length = c(0, 0),
                        textsize   = 6)+

  geom_signif(xmin = 1, xmax = 1.9,
                        y_position = 0.15,
                        annotation = 'p = 0.310',
                        color      = 'black', tip_length = c(0, 0),
                        textsize   = 3)+
    geom_signif(xmin = 2.1, xmax = 5,
                        y_position = 0.15,
                        annotation = 'p = 0.059',
                        color      = 'black', tip_length = c(0, 0),
                        textsize   = 3)
pgr

mod_pgr = lm(mean_pgr~group, data = clust_6_grouped_mdf)
anova(mod_pgr, test = 'Chisq')
pairs(emmeans(mod_pgr, 'group'), adjust = 'none')

ggsave(plot = pgr,
       file = "pgr_specific_targets_6.svg",
       device = "svg",
       units = "in",
       width = 2,
       height = 2,
       path = '/Users/ggraham/Desktop/multiome_poa/Manuscript/Plots/')
