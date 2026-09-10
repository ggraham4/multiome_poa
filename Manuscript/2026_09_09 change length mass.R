#making change length change size plots

library(ggplot2)
library(tidyverse)
library(multcomp)
library(emmeans)

data =read.csv('~/Desktop/multiome_poa/Measures/2025_12_26 all_data.csv')

dat_experiment = subset(data, Condition =='Experiment')

phase = c(
  'F' = 'F',
  'M'  = 'M',
  'D'  = 'I',
  'S'  = 'IP',
  'E'  = 'LI',
  'EP' = 'LIP',
  'NM' = 'NM',
  'NF' = 'NF'
)

dat_experiment$Phase = phase[as.character(dat_experiment$Status)]



dat_experiment_group = dat_experiment%>%
  group_by(Phase)%>%
  summarize(mean_change_mass = mean(as.numeric(Change_Mass)/100),
            mean_change_length = mean(as.numeric(Change_Length)/100),
            se_length = sd(as.numeric(Change_Length)/100)/sqrt(n()),
            se_mass = sd(as.numeric(Change_Mass)/100)/sqrt(n()))



mass_model = lm(Change_Mass~Phase, data = dat_experiment)

anova(mass_model, test = 'Chisq')
p1 = pairs(emmeans(mass_model, 'Phase'), adjust = 'none')
cld_mas = cld(emmeans(mass_model, 'Phase'), Letters = letters, adjust = "none", alpha = 0.05)
cld_mas_df = as.data.frame(cld_mas)

change_mass = ggplot(dat_experiment, aes(x = Phase,
                                         y = as.numeric(Change_Mass)/100))+
  geom_crossbar(data = dat_experiment_group, 
                aes(x = Phase,
                    ymin = mean_change_mass,
                    ymax = mean_change_mass,
                    y = mean_change_mass), 
                inherit.aes = F)+
  geom_errorbar(data = dat_experiment_group, 
                aes(x = Phase, 
                    y = mean_change_mass,
                    ymin = mean_change_mass-se_mass,
                    ymax =mean_change_mass+se_mass),
                inherit.aes = F)+
  geom_point(position = position_jitterdodge(),
             size = 2, 
             shape =1)+
  scale_y_continuous(labels = scales::percent)+
  theme_minimal()+
  labs(x = 'Phase', y = 'Change Mass %')+
    geom_text(data = cld_mas_df, aes(x = Phase,
                                     y = 1.1*max(as.numeric(dat_experiment$Change_Mass)/100),
                                     label = .group), 
          size = 3, inherit.aes = FALSE)
change_mass

ggsave(plot = change_mass,
       file = "change_mass.svg",
       device = "svg",
       units = "in",
       width = 3,
       height = 2,
       path = "Manuscript/Plots/Fig.1/")

Length_model = lm(Change_Length~Phase, data = dat_experiment)

anova(Length_model, test = 'Chisq')
p1 = pairs(emmeans(Length_model, 'Phase'), adjust = 'none')
cld_len = cld(emmeans(Length_model, 'Phase'), Letters = letters, adjust = "none", alpha = 0.05)
cld_len_df = as.data.frame(cld_len)

change_length = ggplot(dat_experiment, aes(x = Phase,
                                         y = as.numeric(Change_Length)/100))+
  geom_crossbar(data = dat_experiment_group, 
                aes(x = Phase,
                    ymin = mean_change_length,
                    ymax = mean_change_length,
                    y = mean_change_length), 
                inherit.aes = F)+
  geom_errorbar(data = dat_experiment_group, 
                aes(x = Phase, 
                    y = mean_change_length,
                    ymin = mean_change_length-se_length,
                    ymax =mean_change_length+se_length),
                inherit.aes = F)+
  geom_point(position = position_jitterdodge(),
             size = 2, 
             shape =1)+
  scale_y_continuous(labels = scales::percent)+
  theme_minimal()+
  labs(x = 'Phase', y = 'Change Length %')+
    geom_text(data = cld_len_df, aes(x = Phase,
                                     y = 1.1*max(as.numeric(dat_experiment$Change_Length)/100),
                                     label = .group), 
          size = 3, inherit.aes = FALSE)
change_length

ggsave(plot = change_length,
       file = "change_length.svg",
       device = "svg",
       units = "in",
       width = 3,
       height = 2,
       path = "Manuscript/Plots/Fig.1/")

dat_volume = subset(data, !is.na(Log10_Volume))
dat_volume$Log10_Volume = as.numeric(dat_volume$Log10_Volume)

dat_volume$Phase = phase[as.character(dat_volume$Status)]

dat_volume_group = dat_volume%>%
  group_by(Phase)%>%
  summarize(mean = mean(Log10_Volume),
            se_vol = sd(Log10_Volume)/sqrt(n()))

vol_model = lm(Log10_Volume~Phase, data = dat_volume)

anova(vol_model, test = 'Chisq')
p1 = pairs(emmeans(vol_model, 'Phase'), adjust = 'none')
cld_vol = cld(emmeans(vol_model, 'Phase'), Letters = letters, adjust = "none", alpha = 0.05)
cld_vol_df = as.data.frame(cld_vol)

volume = ggplot(dat_volume, aes(x = Phase, y = Log10_Volume))+
  geom_point(position = position_jitterdodge(),
             size = 2,
             shape =1)+
    geom_crossbar(data = dat_volume_group, 
                aes(x = Phase,
                    ymin = mean,
                    ymax = mean,
                    y = mean), 
                inherit.aes = F)+
  geom_errorbar(data = dat_volume_group, 
                aes(x = Phase, 
                    y = mean,
                    ymin = mean-se_vol,
                    ymax =mean+se_vol),
                inherit.aes = F)+
      geom_text(data = cld_vol_df, aes(x = Phase,
                                     y = 1.1*max(as.numeric(dat_volume$Log10_Volume)),
                                     label = .group))+

  theme_minimal()+
  labs(y = 'Log10 Gonadal Volume')
volume

#ggsave(plot = volume,
#       file = "volume.svg",
#       device = "svg",
#       units = "in",
#       width = 3,
#       height = 2,
 #      path = "Manuscript/Plots/Fig.1/")
```r
# ============================
# TESTICULAR VOLUME
# ============================

dat_testis = subset(data, !is.na(Log10_Testicular_Estimate))
dat_testis$Log10_Testicular_Estimate = as.numeric(dat_testis$Log10_Testicular_Estimate)

dat_testis$Phase = phase[as.character(dat_testis$Status)]

dat_testis_group = dat_testis %>%
  group_by(Phase) %>%
  summarize(
    mean = mean(Log10_Testicular_Estimate),
    se_testis = sd(Log10_Testicular_Estimate) / sqrt(n())
  )

testis_model = lm(Log10_Testicular_Estimate ~ Phase, data = dat_testis)

anova(testis_model, test = 'Chisq')

p_testis = pairs(
  emmeans(testis_model, 'Phase'),
  adjust = 'none'
)

cld_testis = cld(
  emmeans(testis_model, 'Phase'),
  Letters = letters,
  adjust = "none",
  alpha = 0.05
)

cld_testis_df = as.data.frame(cld_testis)


testis = ggplot(
  dat_testis,
  aes(x = Phase, y = Log10_Testicular_Estimate)
) +
  geom_point(
    position = position_jitterdodge(),
    size = 2,
    shape = 1
  ) +
  geom_crossbar(
    data = dat_testis_group,
    aes(
      x = Phase,
      ymin = mean,
      ymax = mean,
      y = mean
    ),
    inherit.aes = FALSE
  ) +
  geom_errorbar(
    data = dat_testis_group,
    aes(
      x = Phase,
      y = mean,
      ymin = mean - se_testis,
      ymax = mean + se_testis
    ),
    inherit.aes = FALSE
  ) +
  geom_text(
    data = cld_testis_df,
    aes(
      x = Phase,
      y = 1.1 * max(dat_testis$Log10_Testicular_Estimate),
      label = .group
    ),
    inherit.aes = FALSE
  ) +
  theme_minimal() +
  labs(
    x = 'Phase',
    y = 'Log10 Testicular Volume'
  )

testis


# Save testicular plot

ggsave(
  plot = testis,
  file = "testicular_volume.svg",
  device = "svg",
  units = "in",
  width = 3,
  height = 2,
  path = "Manuscript/Plots/Fig.1/"
)



# ============================
# OVARIAN VOLUME
# ============================

dat_ovary = subset(data, !is.na(Log10_Ovarian_Estimate))
dat_ovary$Log10_Ovarian_Estimate = as.numeric(dat_ovary$Log10_Ovarian_Estimate)

dat_ovary$Phase = phase[as.character(dat_ovary$Status)]

dat_ovary_group = dat_ovary %>%
  group_by(Phase) %>%
  summarize(
    mean = mean(Log10_Ovarian_Estimate),
    se_ovary = sd(Log10_Ovarian_Estimate) / sqrt(n())
  )

ovary_model = lm(Log10_Ovarian_Estimate ~ Phase, data = dat_ovary)

anova(ovary_model, test = 'Chisq')

p_ovary = pairs(
  emmeans(ovary_model, 'Phase'),
  adjust = 'none'
)

cld_ovary = cld(
  emmeans(ovary_model, 'Phase'),
  Letters = letters,
  adjust = "none",
  alpha = 0.05
)

cld_ovary_df = as.data.frame(cld_ovary)


ovary = ggplot(
  dat_ovary,
  aes(x = Phase, y = Log10_Ovarian_Estimate)
) +
  geom_point(
    position = position_jitterdodge(),
    size = 2,
    shape = 1
  ) +
  geom_crossbar(
    data = dat_ovary_group,
    aes(
      x = Phase,
      ymin = mean,
      ymax = mean,
      y = mean
    ),
    inherit.aes = FALSE
  ) +
  geom_errorbar(
    data = dat_ovary_group,
    aes(
      x = Phase,
      y = mean,
      ymin = mean - se_ovary,
      ymax = mean + se_ovary
    ),
    inherit.aes = FALSE
  ) +
  geom_text(
    data = cld_ovary_df,
    aes(
      x = Phase,
      y = 1.1 * max(dat_ovary$Log10_Ovarian_Estimate),
      label = .group
    ),
    inherit.aes = FALSE
  ) +
  theme_minimal() +
  labs(
    x = 'Phase',
    y = 'Log10 Ovarian Volume'
  )

ovary


# Save ovarian plot

ggsave(
  plot = ovary,
  file = "ovarian_volume.svg",
  device = "svg",
  units = "in",
  width = 3,
  height = 2,
  path = "Manuscript/Plots/Fig.1/"
)



