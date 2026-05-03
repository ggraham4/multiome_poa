library(ggplot2)
library(lme4)
library(tidyverse)
library(emmeans)
library(ggsignif)
library(multcomp)
`%notin%` = Negate(`%in%`)
set.seed(0)

measures = read.csv('Measures/all_data.csv')
measures$Time_Day_2= as.numeric(measures$Time_Day_2)
measures$length_final_cm= as.numeric(measures$length_final_cm)
measures$Log_11KT = as.numeric(measures$Log_11KT)

status_to_phase = list('D'='I',
                       'E' = 'LI',
                       'EP' = 'LIP',
                       'S' = 'IP',
                       'M'='M',
                       'F'='F',
                       'NF'='NF',
                       'NM'='NM')
measures$Phase = unlist(status_to_phase[measures$Status])
measures$Phase = factor(measures$Phase, levels =c('F','M','I','IP','LI','LIP','NF','NM'))
measures$Dominance = ifelse(measures$Phase %in%c('F','I','LI','NF'), 'Dominant', 'Subordinate')
measures$length_final_cm = as.numeric(measures$length_final_cm)

### length_final_cm  ----
mass_model = lm(length_final_cm~Phase, data = measures)
anova(mass_model, test = 'Chisq')
p1 = pairs(emmeans(mass_model, 'Phase'), adjust = 'none')
cld_mas = cld(emmeans(mass_model, 'Phase'), Letters = letters, adjust = "none", alpha = 0.05)
cld_mas_df = as.data.frame(cld_mas)

mean(measures$length_final_cm[measures$Phase=='F'])/mean(measures$length_final_cm[measures$Phase=='M'])
mean(measures$length_final_cm[measures$Phase=='I'])/mean(measures$length_final_cm[measures$Phase=='IP'])
mean(measures$length_final_cm[measures$Phase=='LI'])/mean(measures$length_final_cm[measures$Phase=='LIP'])
mean(measures$length_final_cm[measures$Phase=='NF'])/mean(measures$length_final_cm[measures$Phase=='NM'])

mas_mean_sd = measures%>%
  group_by(Phase, Dominance)%>%
  summarize(mean  =mean(length_final_cm),
            se = sd(length_final_cm)/sqrt(n()))

mas_plot = ggplot(measures, aes(y = length_final_cm, x = Phase, color = Dominance, shape = Dominance))+
  stat_summary(data = subset(measures, !Phase %in% c('NM', 'NF')),
               geom = 'crossbar', fun = 'mean') +
  geom_errorbar(data = mas_mean_sd, aes(x = Phase, y = mean, ymin = mean-se, ymax = mean+se),
                color ='black',
                width = 0.4, linewidth = 0.5)+
  labs(x = 'Phase', y = 'Mass (g)') +
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(linewidth=1, color = 'black'),
        legend.key = element_blank())+
  geom_vline(xintercept = 2.5)+
    geom_vline(xintercept = 4.5,)+
      geom_vline(xintercept = 6.5)+
  annotate("segment", 
           x = 1.5, xend = 1.5, 
           y = mean(measures$length_final_cm[measures$Phase=='F']), 
           yend = mean(measures$length_final_cm[measures$Phase=='M']),
           arrow = arrow(ends = "both", length = unit(0.2, "cm")))+
  geom_text(x = 2, y = mean(c(mean(measures$length_final_cm[measures$Phase=='M']),
                              mean(measures$length_final_cm[measures$Phase=='F']))),
            label = '3.12x', 
            color = 'black',
            size = 2.5,
             legend.key = element_blank())+
    annotate("segment", 
           x = 3.5, xend = 3.5, 
           y = mean(measures$length_final_cm[measures$Phase=='I']), 
           yend = mean(measures$length_final_cm[measures$Phase=='IP']),
           arrow = arrow(ends = "both", length = unit(0.2, "cm")))+
    geom_text(x = 4, y = mean(c(mean(measures$length_final_cm[measures$Phase=='I']),
                              mean(measures$length_final_cm[measures$Phase=='IP']))),
            label = '2.33x', 
            color = 'black',
            size = 2.5,
             legend.key = element_blank())+
      annotate("segment", 
           x = 5.5, xend = 5.5, 
           y = mean(measures$length_final_cm[measures$Phase=='LI']), 
           yend = mean(measures$length_final_cm[measures$Phase=='LIP']),
           arrow = arrow(ends = "both", length = unit(0.2, "cm")))+
      geom_text(x = 6, y = mean(c(mean(measures$length_final_cm[measures$Phase=='LI']),
                              mean(measures$length_final_cm[measures$Phase=='LIP']))),
            label = '2.19x', 
            color = 'black',
            size = 2.5,
             legend.key = element_blank())+
        annotate("segment", 
           x = 7.5, xend = 7.5, 
           y = mean(measures$length_final_cm[measures$Phase=='NF']), 
           yend = mean(measures$length_final_cm[measures$Phase=='NM']),
           arrow = arrow(ends = "both", length = unit(0.2, "cm")))+
        geom_text(x = 8, y = mean(c(mean(measures$length_final_cm[measures$Phase=='NF']),
                              mean(measures$length_final_cm[measures$Phase=='NM']))),
            label = '2.52x', 
            color = 'black',
            size = 2.5,
             legend.key = element_blank())+
    geom_point( size = 2, position = position_jitterdodge(0.25)) +
  geom_text(data = cld_mas_df, aes(x = Phase, y = 1.1*max(measures$length_final_cm), label = .group), 
          size = 3, inherit.aes = FALSE)+
  theme(legend.position = 'none')
mas_plot

#ggsave(plot = mas_plot,
 #      file = "mass_plot.svg",
  #     device = "svg",
   #    units = "in",
    #   width = 3,
     #  height = 2,
      # path = '/Users/ggraham/Desktop/multiome_poa/Manuscript/Plots/Supp 1')


### length final _cm ####
measures$length_final_cm = as.numeric(measures$length_final_cm)

len_model = lm(length_final_cm~Phase, data = measures)
anova(len_model, test = 'Chisq')
p1 = pairs(emmeans(len_model, 'Phase'), adjust = 'none')
cld_len = cld(emmeans(len_model, 'Phase'), Letters = letters, adjust = "none", alpha = 0.05)
cld_len_df = as.data.frame(cld_len)

mean(measures$length_final_cm[measures$Phase=='F'])/mean(measures$length_final_cm[measures$Phase=='M'])
mean(measures$length_final_cm[measures$Phase=='I'])/mean(measures$length_final_cm[measures$Phase=='IP'])
mean(measures$length_final_cm[measures$Phase=='LI'])/mean(measures$length_final_cm[measures$Phase=='LIP'])
mean(measures$length_final_cm[measures$Phase=='NF'])/mean(measures$length_final_cm[measures$Phase=='NM'])

len_mean_sd = measures%>%
  group_by(Phase, Dominance)%>%
  summarize(mean  =mean(length_final_cm),
            se = sd(length_final_cm)/sqrt(n()))

len_plot = ggplot(measures, aes(y = length_final_cm, x = Phase, color = Dominance, shape = Dominance))+
  stat_summary(data = subset(measures, !Phase %in% c('NM', 'NF')),
               geom = 'crossbar', fun = 'mean') +
  geom_errorbar(data = len_mean_sd, aes(x = Phase, y = mean, ymin = mean-se, ymax = mean+se),
                color ='black',
                width = 0.4, linewidth = 0.5)+
  labs(x = 'Phase', y = 'Length (cm)') +
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(linewidth=1, color = 'black'),
        legend.key = element_blank())+
  geom_vline(xintercept = 2.5)+
    geom_vline(xintercept = 4.5,)+
      geom_vline(xintercept = 6.5)+
  annotate("segment", 
           x = 1.5, xend = 1.5, 
           y = mean(measures$length_final_cm[measures$Phase=='F']), 
           yend = mean(measures$length_final_cm[measures$Phase=='M']),
           arrow = arrow(ends = "both", length = unit(0.2, "cm")))+
  geom_text(x = 2, y = mean(c(mean(measures$length_final_cm[measures$Phase=='M']),
                              mean(measures$length_final_cm[measures$Phase=='F']))),
            label = '1.44x', 
            color = 'black',
            size = 2.5,
             legend.key = element_blank())+
    annotate("segment", 
           x = 3.5, xend = 3.5, 
           y = mean(measures$length_final_cm[measures$Phase=='I']), 
           yend = mean(measures$length_final_cm[measures$Phase=='IP']),
           arrow = arrow(ends = "both", length = unit(0.2, "cm")))+
    geom_text(x = 4, y = mean(c(mean(measures$length_final_cm[measures$Phase=='I']),
                              mean(measures$length_final_cm[measures$Phase=='IP']))),
            label = '1.31x', 
            color = 'black',
            size = 2.5,
             legend.key = element_blank())+
      annotate("segment", 
           x = 5.5, xend = 5.5, 
           y = mean(measures$length_final_cm[measures$Phase=='LI']), 
           yend = mean(measures$length_final_cm[measures$Phase=='LIP']),
           arrow = arrow(ends = "both", length = unit(0.2, "cm")))+
      geom_text(x = 6, y = mean(c(mean(measures$length_final_cm[measures$Phase=='LI']),
                              mean(measures$length_final_cm[measures$Phase=='LIP']))),
            label = '1.28x', 
            color = 'black',
            size = 2.5,
             legend.key = element_blank())+
        annotate("segment", 
           x = 7.5, xend = 7.5, 
           y = mean(measures$length_final_cm[measures$Phase=='NF']), 
           yend = mean(measures$length_final_cm[measures$Phase=='NM']),
           arrow = arrow(ends = "both", length = unit(0.2, "cm")))+
        geom_text(x = 8, y = mean(c(mean(measures$length_final_cm[measures$Phase=='NF']),
                              mean(measures$length_final_cm[measures$Phase=='NM']))),
            label = '1.28x', 
            color = 'black',
            size = 2.5,
             legend.key = element_blank())+
    geom_point( size = 2, position = position_jitterdodge(0.25)) +
  geom_text(data = cld_len_df, aes(x = Phase, y = 1.1*max(measures$length_final_cm), label = .group), 
          size = 3, inherit.aes = FALSE)+
  theme(legend.position = 'none')
len_plot

ggsave(plot = len_plot,
       file = "length_plot.svg",
       device = "svg",
       units = "in",
       width = 3,
       height = 2,
       path = '/Users/ggraham/Desktop/multiome_poa/Manuscript/Plots/Supp 1')





