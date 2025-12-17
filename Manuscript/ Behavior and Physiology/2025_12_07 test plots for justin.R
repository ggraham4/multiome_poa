library(ggplot2)
library(lme4)
library(tidyverse)
library(emmeans)
library(ggsignif)

measures = read.csv('Measures/all_data.csv')
measures$Time_Day_2= as.numeric(measures$Time_Day_2)
measures$Behaviors_Day_2= as.numeric(measures$Behaviors_Day_2)
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

time_model = lm(Time_Day_2~Phase, data = measures)
anova(time_model, test = 'Chisq')
p = pairs(emmeans(time_model, 'Phase'), adjust = 'none')

set.seed(0)
time_plot_minimal = ggplot(measures, aes(y = Time_Day_2, x = Phase)) +
  geom_point(data = subset(measures, !Phase %in% c('NM', 'NF')),
             size = 1, position = position_jitterdodge(0.5)) +
  geom_point(data = subset(measures, Phase %in% c('NM', 'NF')),
             size = 1) +
  geom_boxplot(data = subset(measures, !Phase %in% c('NM', 'NF')),
               alpha = 0) +
  labs(x = 'Phase', y = 'Time in nest (s)') +
  theme_minimal()
time_plot_minimal

ggsave(plot = time_plot_minimal,
       file = "time_plot_minimal.svg",
       device = "svg",
       units = "in",
       width = 2.5,
       height = 2,
       path = "Manuscript/Plots/Fig.1/Time")

time_plot_minimal_signif = ggplot(measures, aes(y = Time_Day_2, x = Phase)) +
  geom_point(data = subset(measures, !Phase %in% c('NM', 'NF')),
             size = 1, position = position_jitterdodge(0.5)) +
  geom_point(data = subset(measures, Phase %in% c('NM', 'NF')),
             size = 1) +
  geom_boxplot(data = subset(measures, !Phase %in% c('NM', 'NF')),
               alpha = 0) +
  labs(x = 'Phase', y = 'Time in nest (s)') +
  theme_minimal()+
  geom_signif(xmin = c(1), xmax = c(1.9),
              y_position = c(max(measures$Time_Day_2*1.05)),
              annotation =c("*"),
              color = "black", tip_length = c(0,0),
              textsize=6,             
              vjust = 0.5)+
    geom_signif(xmin = c(1), xmax = c(3),
              y_position = c(max(measures$Time_Day_2*1.2)),
              annotation =c("p=0.066"),
              color = "black", tip_length = c(0,0),
              textsize=3, 
              vjust = -0.2)+
    geom_signif(xmin = c(2.1), xmax = c(4),
              y_position = c(max(measures$Time_Day_2*1.05)),
              annotation =c("**"),
              color = "black", tip_length = c(0,0),
              textsize=6,
              vjust = 0.5)+
  ylim(0, 1.25*max(measures$Time_Day_2))

time_plot_minimal_signif

ggsave(plot = time_plot_minimal_signif,
       file = "time_plot_minimal_signif.svg",
       device = "svg",
       units = "in",
       width = 2.5,
       height = 2,
       path = "Manuscript/Plots/Fig.1/Time")

library(multcomp)
cld_results = cld(emmeans(time_model, 'Phase'), Letters = letters, adjust = "none", alpha = 0.05)
cld_df = as.data.frame(cld_results)

cld_df$.group = trimws(cld_df$.group)

time_plot_minimal_letters = ggplot(measures, aes(y = Time_Day_2, x = Phase)) +
  geom_point(data = subset(measures, !Phase %in% c('NM', 'NF')),
             size = 1, position = position_jitterdodge(0.5)) +
  geom_point(data = subset(measures, Phase %in% c('NM', 'NF')),
             size = 1) +
  geom_boxplot(data = subset(measures, !Phase %in% c('NM', 'NF')),
               alpha = 0) +
  labs(x = 'Phase', y = 'Time in nest (s)') +
  theme_minimal()+
  geom_text(data = cld_df, 
            aes(x = Phase, y = 650, label = .group),
            size = 4) +
  ylim(0, 1.25*max(measures$Time_Day_2))

time_plot_minimal_letters

ggsave(plot = time_plot_minimal_letters,
       file = "time_plot_minimal_letters.svg",
       device = "svg",
       units = "in",
       width = 2.5,
       height = 2,
       path = "Manuscript/Plots/Fig.1/Time")



measures$Dominance = ifelse(measures$Phase %in%c('F','I','LI','NF'), 'Dominant', 'Subordinate')

time_plot_minimal_dominance = ggplot(measures, aes(y = Time_Day_2, x = Phase, color = Dominance, shape = Dominance) )+
  geom_point(data = subset(measures, !Phase %in% c('NM', 'NF')),
             size = 1, position = position_jitterdodge(0.5)) +
  geom_point(data = subset(measures, Phase %in% c('NM', 'NF')),
             size = 1) +
  geom_boxplot(data = subset(measures, !Phase %in% c('NM', 'NF')),
               alpha = 0) +
  labs(x = 'Phase', y = 'Time in nest (s)') +
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(linewidth=1, color = 'black'),
        legend.key = element_blank())+
  geom_vline(xintercept = 2.5)+
    geom_vline(xintercept = 4.5,)+
      geom_vline(xintercept = 6.5)
time_plot_minimal_dominance

ggsave(plot = time_plot_minimal_dominance,
       file = "time_plot_minimal_dominance.svg",
       device = "svg",
       units = "in",
       width = 4.5,
       height = 2,
       path = "Manuscript/Plots/Fig.1/Time")

mean(measures$Time_Day_2[measures$Phase=='M'])/mean(measures$Time_Day_2[measures$Phase=='F'])
mean(measures$Time_Day_2[measures$Phase=='I'])/mean(measures$Time_Day_2[measures$Phase=='IP'])
mean(measures$Time_Day_2[measures$Phase=='LIP'])/mean(measures$Time_Day_2[measures$Phase=='LI'])
mean(measures$Time_Day_2[measures$Phase=='NM'])/mean(measures$Time_Day_2[measures$Phase=='NF'])

time_plot_jennystyle_dominance = ggplot(measures, aes(y = Time_Day_2, x = Phase, color = Dominance, shape = Dominance) )+
  geom_point(
             size = 2, position = position_jitterdodge(0.25)) +
  stat_summary(data = subset(measures, !Phase %in% c('NM', 'NF')),
               geom = 'crossbar', fun = 'mean') +
  labs(x = 'Phase', y = 'Time in nest (s)') +
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
           y = mean(measures$Time_Day_2[measures$Phase=='F']), 
           yend = mean(measures$Time_Day_2[measures$Phase=='M']),
           arrow = arrow(ends = "both", length = unit(0.2, "cm")))+
  geom_text(x = 2, y = mean(c(mean(measures$Time_Day_2[measures$Phase=='F']),
                              mean(measures$Time_Day_2[measures$Phase=='M']))),
            label = '1.69x', 
            color = 'black',
            size = 2.5,
             legend.key = element_blank())+
    annotate("segment", 
           x = 3.5, xend = 3.5, 
           y = mean(measures$Time_Day_2[measures$Phase=='I']), 
           yend = mean(measures$Time_Day_2[measures$Phase=='IP']),
           arrow = arrow(ends = "both", length = unit(0.2, "cm")))+
    geom_text(x = 4, y = mean(c(mean(measures$Time_Day_2[measures$Phase=='I']),
                              mean(measures$Time_Day_2[measures$Phase=='IP']))),
            label = '2.14x', 
            color = 'black',
            size = 2.5,
             legend.key = element_blank())+
      annotate("segment", 
           x = 5.5, xend = 5.5, 
           y = mean(measures$Time_Day_2[measures$Phase=='LI']), 
           yend = mean(measures$Time_Day_2[measures$Phase=='LIP']),
           arrow = arrow(ends = "both", length = unit(0.2, "cm")))+
      geom_text(x = 6, y = mean(c(mean(measures$Time_Day_2[measures$Phase=='LI']),
                              mean(measures$Time_Day_2[measures$Phase=='LIP']))),
            label = '1.84x', 
            color = 'black',
            size = 2.5,
             legend.key = element_blank())+
        annotate("segment", 
           x = 7.5, xend = 7.5, 
           y = mean(measures$Time_Day_2[measures$Phase=='NF']), 
           yend = mean(measures$Time_Day_2[measures$Phase=='NM']),
           arrow = arrow(ends = "both", length = unit(0.2, "cm")))+
        geom_text(x = 8, y = mean(c(mean(measures$Time_Day_2[measures$Phase=='NF']),
                              mean(measures$Time_Day_2[measures$Phase=='NM']))),
            label = '10.67x', 
            color = 'black',
            size = 2.5,
             legend.key = element_blank())


time_plot_jennystyle_dominance

ggsave(plot = time_plot_jennystyle_dominance,
       file = "time_plot_jennystyle_dominance.svg",
       device = "svg",
       units = "in",
       width = 4.5,
       height = 2,
       path = "Manuscript/Plots/Fig.1/Time")


