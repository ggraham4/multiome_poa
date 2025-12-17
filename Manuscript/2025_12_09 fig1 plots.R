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
measures$Dominance = ifelse(measures$Phase %in%c('F','I','LI','NF'), 'Dominant', 'Subordinate')
measures$Behaviors_Day_2 = as.numeric(measures$Behaviors_Day_2)
### Sum Behaviors ----
beh_model = lm(Behaviors_Day_2~Phase, data = measures)
anova(beh_model, test = 'Chisq')
p1 = pairs(emmeans(beh_model, 'Phase'), adjust = 'none')
cld_beh = cld(emmeans(beh_model, 'Phase'), Letters = letters, adjust = "none", alpha = 0.05)
cld_beh_df = as.data.frame(cld_beh)

mean(measures$Behaviors_Day_2[measures$Phase=='M'])/mean(measures$Behaviors_Day_2[measures$Phase=='F'])
mean(measures$Behaviors_Day_2[measures$Phase=='I'])/mean(measures$Behaviors_Day_2[measures$Phase=='IP'])
mean(measures$Behaviors_Day_2[measures$Phase=='LIP'])/mean(measures$Behaviors_Day_2[measures$Phase=='LI'])
mean(measures$Behaviors_Day_2[measures$Phase=='NM'])/mean(measures$Behaviors_Day_2[measures$Phase=='NF'])

beh_mean_sd = measures%>%
  group_by(Phase, Dominance)%>%
  summarize(mean  =mean(Behaviors_Day_2),
            se = sd(Behaviors_Day_2)/sqrt(n()))

beh_plot = ggplot(measures, aes(y = Behaviors_Day_2, x = Phase, color = Dominance, shape = Dominance))+
  stat_summary(data = subset(measures, !Phase %in% c('NM', 'NF')),
               geom = 'crossbar', fun = 'mean') +
  geom_errorbar(data = beh_mean_sd, aes(x = Phase, y = mean, ymin = mean-se, ymax = mean+se),
                color ='black',
                width = 0.4, linewidth = 0.5)+
  labs(x = 'Phase', y = 'Parental Behaviors') +
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
           y = mean(measures$Behaviors_Day_2[measures$Phase=='F']), 
           yend = mean(measures$Behaviors_Day_2[measures$Phase=='M']),
           arrow = arrow(ends = "both", length = unit(0.2, "cm")))+
  geom_text(x = 2, y = mean(c(mean(measures$Behaviors_Day_2[measures$Phase=='F']),
                              mean(measures$Behaviors_Day_2[measures$Phase=='M']))-20),
            label = '2.21x', 
            color = 'black',
            size = 2.5,
             legend.key = element_blank())+
    annotate("segment", 
           x = 3.5, xend = 3.5, 
           y = mean(measures$Behaviors_Day_2[measures$Phase=='I']), 
           yend = mean(measures$Behaviors_Day_2[measures$Phase=='IP']),
           arrow = arrow(ends = "both", length = unit(0.2, "cm")))+
    geom_text(x = 4, y = mean(c(mean(measures$Behaviors_Day_2[measures$Phase=='I']),
                              mean(measures$Behaviors_Day_2[measures$Phase=='IP']))),
            label = '3.88x', 
            color = 'black',
            size = 2.5,
             legend.key = element_blank())+
      annotate("segment", 
           x = 5.5, xend = 5.5, 
           y = mean(measures$Behaviors_Day_2[measures$Phase=='LI']), 
           yend = mean(measures$Behaviors_Day_2[measures$Phase=='LIP']),
           arrow = arrow(ends = "both", length = unit(0.2, "cm")))+
      geom_text(x = 6, y = mean(c(mean(measures$Behaviors_Day_2[measures$Phase=='LI']),
                              mean(measures$Behaviors_Day_2[measures$Phase=='LIP']))),
            label = '7.64x', 
            color = 'black',
            size = 2.5,
             legend.key = element_blank())+
        annotate("segment", 
           x = 7.5, xend = 7.5, 
           y = mean(measures$Behaviors_Day_2[measures$Phase=='NF']), 
           yend = mean(measures$Behaviors_Day_2[measures$Phase=='NM']),
           arrow = arrow(ends = "both", length = unit(0.2, "cm")))+
        geom_text(x = 8, y = mean(c(mean(measures$Behaviors_Day_2[measures$Phase=='NF']),
                              mean(measures$Behaviors_Day_2[measures$Phase=='NM']))),
            label = '11.08x', 
            color = 'black',
            size = 2.5,
             legend.key = element_blank())+
    geom_point( size = 2, position = position_jitterdodge(0.25)) +
  geom_text(data = cld_beh_df, aes(x = Phase, y = 1.1*max(measures$Behaviors_Day_2), label = .group), 
          size = 3, inherit.aes = FALSE)+
  theme(legend.position = 'none')
beh_plot

ggsave(plot = beh_plot,
       file = "beh_plot.svg",
       device = "svg",
       units = "in",
       width = 3,
       height = 2,
       path = "Manuscript/Plots/Fig.1/Behaviors")




### Time ----
tim_mean_sd = measures%>%
  group_by(Phase, Dominance)%>%
  summarize(mean  =mean(Time_Day_2),
            se = sd(Time_Day_2)/sqrt(n()))

time_model = lm(Time_Day_2~Phase, data = measures)
anova(time_model, test = 'Chisq')
p2 = pairs(emmeans(time_model, 'Phase'), adjust = 'none')
cld_results = cld(emmeans(time_model, 'Phase'), Letters = letters, adjust = "none", alpha = 0.05)
cld_df = as.data.frame(cld_results)

mean(measures$Time_Day_2[measures$Phase=='M'])/mean(measures$Time_Day_2[measures$Phase=='F'])
mean(measures$Time_Day_2[measures$Phase=='I'])/mean(measures$Time_Day_2[measures$Phase=='IP'])
mean(measures$Time_Day_2[measures$Phase=='LIP'])/mean(measures$Time_Day_2[measures$Phase=='LI'])
mean(measures$Time_Day_2[measures$Phase=='NM'])/mean(measures$Time_Day_2[measures$Phase=='NF'])

time_plot = ggplot(measures, aes(y = Time_Day_2, x = Phase, color = Dominance, shape = Dominance))+
  stat_summary(data = subset(measures, !Phase %in% c('NM', 'NF')),
               geom = 'crossbar', fun = 'mean') +
    geom_errorbar(data = tim_mean_sd, aes(x = Phase, y = mean, ymin = mean-se, ymax = mean+se),
                color ='black',
                width = 0.4, linewidth = 0.5)+
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
             legend.key = element_blank())+
    geom_point( size = 2, position = position_jitterdodge(0.25)) +
  geom_text(data = cld_df, aes(x = Phase, y = 1.1*max(measures$Time_Day_2), label = .group), 
          size = 3, inherit.aes = FALSE)+
  
  theme(legend.position = 'none')
time_plot

ggsave(plot = time_plot,
       file = "time_plot_nolegend.svg",
       device = "svg",
       units = "in",
       width = 3,
       height = 2,
       path = "Manuscript/Plots/Fig.1/Time")


### 11KT ----
kt = measures%>%
  group_by(Phase, Dominance)%>%
  subset(!is.na(Log_11KT))%>%
  summarize(mean  =mean(Log_11KT,),
            se = sd(Log_11KT)/sqrt(n()))

kt_model = lm(Log_11KT~Phase, data = measures)
anova(kt_model, test = 'Chisq')
p3 = pairs(emmeans(kt_model, 'Phase'), adjust = 'none')
cld_kt = cld(emmeans(kt_model, 'Phase'), Letters = letters, adjust = "none", alpha = 0.05)
cld_df_kt = as.data.frame(cld_kt)

  mean(measures$Log_11KT[measures$Phase=='M']%>%na.omit())/mean(measures$Log_11KT[measures$Phase=='F']%>%na.omit())
mean(measures$Log_11KT[measures$Phase=='I'])/mean(measures$Log_11KT[measures$Phase=='IP'])
mean(measures$Log_11KT[measures$Phase=='LIP'])/mean(measures$Log_11KT[measures$Phase=='LI'])
mean(measures$Log_11KT[measures$Phase=='NM'])/mean(measures$Log_11KT[measures$Phase=='NF'])

kt_plot = ggplot(measures, aes(y = Log_11KT, x = Phase, color = Dominance, shape = Dominance))+
  stat_summary(data = subset(measures, !Phase %in% c('NM', 'NF')),
               geom = 'crossbar', fun = 'mean') +
    geom_errorbar(data = kt, aes(x = Phase, y = mean, ymin = mean-se, ymax = mean+se),
                color ='black',
                width = 0.4, linewidth = 0.5)+
  labs(x = 'Phase', y = 'Log10 11KT (pg/ml)') +
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
           y = mean(measures$Log_11KT[measures$Phase=='F']%>%na.omit()), 
           yend = mean(measures$Log_11KT[measures$Phase=='M']%>%na.omit()),
           arrow = arrow(ends = "both", length = unit(0.2, "cm")))+
  geom_text(x = 2, y = mean(c(mean(measures$Log_11KT[measures$Phase=='F' & !is.na(measures$Log_11KT)]),
                              mean(measures$Log_11KT[measures$Phase=='M'& !is.na(measures$Log_11KT)]))),
            label = '1.85x', 
            color = 'black',
            size = 2.5,
             legend.key = element_blank())+
    annotate("segment", 
           x = 3.5, xend = 3.5, 
           y = mean(measures$Log_11KT[measures$Phase=='I']), 
           yend = mean(measures$Log_11KT[measures$Phase=='IP']),
           arrow = arrow(ends = "both", length = unit(0.2, "cm")))+
    geom_text(x = 4, y = 0.05+mean(c(mean(measures$Log_11KT[measures$Phase=='I']),
                              mean(measures$Log_11KT[measures$Phase=='IP']))),
            label = '1.11x', 
            color = 'black',
            size = 2.5,
             legend.key = element_blank())+
      annotate("segment", 
           x = 5.5, xend = 5.5, 
           y = mean(measures$Log_11KT[measures$Phase=='LI']), 
           yend = mean(measures$Log_11KT[measures$Phase=='LIP']),
           arrow = arrow(ends = "both", length = unit(0.2, "cm")))+
      geom_text(x = 6, y = mean(c(mean(measures$Log_11KT[measures$Phase=='LI']),
                              mean(measures$Log_11KT[measures$Phase=='LIP']))),
            label = '1.69x', 
            color = 'black',
            size = 2.5,
             legend.key = element_blank())+
        annotate("segment", 
           x = 7.5, xend = 7.5, 
           y = mean(measures$Log_11KT[measures$Phase=='NF']), 
           yend = mean(measures$Log_11KT[measures$Phase=='NM']),
           arrow = arrow(ends = "both", length = unit(0.2, "cm")))+
        geom_text(x = 8, y = mean(c(mean(measures$Log_11KT[measures$Phase=='NF']),
                              mean(measures$Log_11KT[measures$Phase=='NM']))),
            label = '1.35x', 
            color = 'black',
            size = 2.5,
             legend.key = element_blank())+
    geom_point( size = 2, position = position_jitterdodge(0.25)) +
  geom_text(data = cld_df_kt, aes(x = Phase, y = 1.1*max(measures$Log_11KT%>%na.omit()), label = .group), 
          size = 3, inherit.aes = FALSE)+
  
  theme(legend.position = 'none')
kt_plot

ggsave(plot = kt_plot,
       file = "kt_plot.svg",
       device = "svg",
       units = "in",
       width = 3,
       height = 2,
       path = "Manuscript/Plots/Fig.1/kt")


###---- Percent Testicular -----

pct_test = measures%>%
  group_by(Phase, Dominance)%>%
  subset(!is.na(Percent_Testicular))%>%
  summarize(mean  =mean(Percent_Testicular,),
            se = sd(Percent_Testicular)/sqrt(n()))

test_model = lm(Percent_Testicular~Phase, data = measures)
anova(test_model, test = 'Chisq')
p4 = pairs(emmeans(test_model, 'Phase'), adjust = 'none')
cld_pct_test = cld(emmeans(test_model, 'Phase'), Letters = letters, adjust = "none", alpha = 0.05)
cld_pct_test_df = as.data.frame(cld_pct_test)

  #mean(measures$Percent_Testicular[measures$Phase=='M']%>%na.omit())/mean(measures$Log_11KT[measures$Phase=='F']%>%na.omit())
mean(measures$Percent_Testicular[measures$Phase=='I'])/mean(measures$Percent_Testicular[measures$Phase=='IP'])
mean(measures$Percent_Testicular[measures$Phase=='LIP'])/mean(measures$Percent_Testicular[measures$Phase=='LI'])
#mean(measures$Percent_Testicular[measures$Phase=='NM'])/mean(measures$Percent_Testicular[measures$Phase=='NF'])

test_plot = ggplot(subset(measures, !Phase %in%c('F','NF')), aes(y = Percent_Testicular, x = as.character(Phase)))+
  stat_summary(data = subset(measures, !Phase %in% c('NM', 'NF')),
               geom = 'crossbar', fun = 'mean') +
    geom_errorbar(data = pct_test, aes(x = as.character(Phase), y = mean, ymin = mean-se, ymax = mean+se),
                color ='black',
                width = 0.4, linewidth = 0.5)+
  labs(x = 'Phase', y = 'Testicular %') +
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(linewidth=1, color = 'black'),
        legend.key = element_blank())+
    geom_point( size = 2, position = position_jitterdodge(0.25), shape = 1) +
  geom_text(data = cld_pct_test_df, aes(x = Phase, y = 1.1*max(measures$Percent_Testicular%>%na.omit()), label = .group), 
          size = 3, inherit.aes = FALSE)+
  theme(legend.position = 'none')+
    scale_x_discrete(limits = c('M','I','IP','LI','LIP','NM')) +
        scale_y_continuous(labels = scales::percent)
test_plot

ggsave(plot = test_plot,
       file = "test_plot.svg",
       device = "svg",
       units = "in",
       width = 2.5,
       height = 2,
       path = "Manuscript/Plots/Fig.1/Testicular")

#### ovarian ----
pct_ov = measures%>%
  group_by(Phase, Dominance)%>%
  subset(!is.na(Percent_Ovarian))%>%
  summarize(mean  =mean(Percent_Ovarian,),
            se = sd(Percent_Ovarian)/sqrt(n()))

ov_model = lm(Percent_Ovarian~Phase, data = measures)
anova(ov_model, test = 'Chisq')
p5 = pairs(emmeans(ov_model, 'Phase'), adjust = 'none')
cld_pct_ov = cld(emmeans(ov_model, 'Phase'), Letters = letters, adjust = "none", alpha = 0.05)
cld_pct_ov_df = as.data.frame(cld_pct_ov)

  #mean(measures$Percent_Testicular[measures$Phase=='M']%>%na.omit())/mean(measures$Log_11KT[measures$Phase=='F']%>%na.omit())
mean(measures$Percent_Ovarian[measures$Phase=='I'])/mean(measures$Percent_Ovarian[measures$Phase=='IP'])
mean(measures$Percent_Ovarian[measures$Phase=='LIP'])/mean(measures$Percent_Ovarian[measures$Phase=='LI'])
#mean(measures$Percent_Testicular[measures$Phase=='NM'])/mean(measures$Percent_Testicular[measures$Phase=='NF'])

ov_plot = ggplot(subset(measures, !Phase %in%c('F','NF')), aes(y = Percent_Ovarian, x = as.character(Phase)))+
  stat_summary(data = subset(measures, !Phase %in% c('NM', 'NF')),
               geom = 'crossbar', fun = 'mean') +
    geom_errorbar(data = pct_ov, aes(x = as.character(Phase), y = mean, ymin = mean-se, ymax = mean+se),
                color ='black',
                width = 0.4, linewidth = 0.5)+
  labs(x = 'Phase', y = 'Ovarian %') +
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(linewidth=1, color = 'black'),
        legend.key = element_blank())+
    geom_point( size = 2, position = position_jitterdodge(0.25), shape = 1) +
  geom_text(data = cld_pct_ov_df, aes(x = Phase, y = 1.1*max(measures$Percent_Ovarian%>%na.omit()), label = .group), 
          size = 3, inherit.aes = FALSE)+
  theme(legend.position = 'none')+
    scale_x_discrete(limits = c('M','I','IP','LI','LIP','NM')) +
        scale_y_continuous(labels = scales::percent)
ov_plot

ggsave(plot = ov_plot,
       file = "ov_plot.svg",
       device = "svg",
       units = "in",
       width = 2.5,
       height = 2,
       path = "Manuscript/Plots/Fig.1/Ovarian")

### volume ----
measures$Log10_Volume <- as.numeric(measures$Log10_Volume)
ev_2.5x_model <- lm(Log10_Volume~Phase+length_final_cm, data = measures)
ev_2.5x_pairs <- as.data.frame(pairs(emmeans(ev_2.5x_model, 'Phase'), adjust = 'none'))

measures_gonad= subset(measures, Phase %notin% c('NF','F'))
  
residual_model <- lm(Log10_Volume~length_final_cm, data = measures_gonad)

measures_gonad$corrected_length =residual_model$residuals
anova(ev_2.5x_model,test = 'Chisq')

corrected_vol = measures_gonad%>%
  group_by(Phase)%>%
  subset(!is.na(corrected_length))%>%
  summarize(mean  =mean(corrected_length,),
            se = sd(corrected_length)/sqrt(n()))

cld_pct_vol = cld(emmeans(ev_2.5x_model, 'Phase'), Letters = letters, adjust = "none", alpha = 0.05)
cld_pct_vol_df = as.data.frame(cld_pct_vol)

vol_plot = ggplot(subset(measures_gonad, !Phase %in%c('F','NF')), aes(y = corrected_length, x = as.character(Phase)))+
  stat_summary(data = subset(measures_gonad, !Phase %in% c('NM', 'NF')),
               geom = 'crossbar', fun = 'mean') +
    geom_errorbar(data = corrected_vol, aes(x = as.character(Phase), y = mean, ymin = mean-se, ymax = mean+se),
                color ='black',
                width = 0.4, linewidth = 0.5)+
  labs(x = 'Phase', y = 'Corrected Volume') +
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(linewidth=1, color = 'black'),
        legend.key = element_blank())+
    geom_point( size = 2, position = position_jitterdodge(0.25), shape = 1) +
  geom_text(data = cld_pct_vol_df, aes(x = Phase,
                                      y = 1.15*max(measures_gonad$corrected_length%>%na.omit()),
                                      label = .group), 
          size = 3, inherit.aes = FALSE)+
  theme(legend.position = 'none')+
    scale_x_discrete(limits = c('M','I','IP','LI','LIP','NM')) 
vol_plot

ggsave(plot = vol_plot,
       file = "vol_plot.svg",
       device = "svg",
       units = "in",
       width = 2.5,
       height = 2,
       path = "Manuscript/Plots/Fig.1/Volume")



###
len = measures%>%
  group_by(Phase, Dominance)%>%
  subset(!is.na(length_final_cm))%>%
  summarize(mean  =mean(length_final_cm,),
            se = sd(length_final_cm)/sqrt(n()))

len_plot =ggplot(measures, aes(y = length_final_cm, x = Phase, color = Dominance, shape = Dominance))+
  stat_summary(data = subset(measures, !Phase %in% c('NM', 'NF')),
               geom = 'crossbar', fun = 'mean') +
    geom_errorbar(data = len, aes(x = Phase, y = mean, ymin = mean-se, ymax = mean+se),
                color ='black',
                width = 0.4, linewidth = 0.5)+
  labs(x = 'Phase', y = 'Length (cm)') +
  theme_bw()+
    geom_point( size = 2, position = position_jitterdodge(0.25)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(linewidth=1, color = 'black'),
        legend.key = element_blank())+
    geom_vline(xintercept = 2.5)+
    geom_vline(xintercept = 4.5,)+
      geom_vline(xintercept = 6.5)+
      theme(legend.position = 'none')


ggsave(plot = len_plot,
       file = "len_plot.svg",
       device = "svg",
       units = "in",
       width = 2.5,
       height = 2,
       path = "Manuscript/Plots/Fig.1/Length")



mass = measures%>%
  group_by(Phase, Dominance)%>%
  subset(!is.na(mass_final_cm))%>%
  summarize(mean  =mean(mass_final_cm,),
            se = sd(mass_final_cm)/sqrt(n()))

massplot = ggplot(measures, aes(y = mass_final_cm, x = Phase, color = Dominance, shape = Dominance))+
  stat_summary(data = subset(measures, !Phase %in% c('NM', 'NF')),
               geom = 'crossbar', fun = 'mean') +
    geom_errorbar(data = mass, aes(x = Phase, y = mean, ymin = mean-se, ymax = mean+se),
                color ='black',
                width = 0.4, linewidth = 0.5)+
  labs(x = 'Phase', y = 'Mass (g)') +
  theme_bw()+
    geom_point( size = 2, position = position_jitterdodge(0.25)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(linewidth=1, color = 'black'),
        legend.key = element_blank())+
    geom_vline(xintercept = 2.5)+
    geom_vline(xintercept = 4.5,)+
      geom_vline(xintercept = 6.5)+
    theme(legend.position = 'none')


ggsave(plot = massplot,
       file = "mass_plot.svg",
       device = "svg",
       units = "in",
       width = 2.5,
       height = 2,
       path = "Manuscript/Plots/Fig.1/mass")








