
library(ggplot2)
library(tidyverse)
library(multcomp)
library(emmeans)
populate_statistics <- function(...) {

  models <- list(...)

  results <- lapply(
    names(models),
    function(model_name) {

      fit <- models[[model_name]]

      # ------------------------------------------------------
      # Get model formula
      # ------------------------------------------------------

      model_formula <- formula(fit)

      predictor <- all.vars(model_formula)[2]


      # ------------------------------------------------------
      # ANOVA
      # ------------------------------------------------------

      aov_tab <- anova(fit)


      # ------------------------------------------------------
      # Create output row
      # ------------------------------------------------------

      out <- tibble(
        model = paste0(
          "lm(",
          deparse(model_formula),
          ")"
        ),
        figure = NA_character_,
        anova_p.value = aov_tab$`Pr(>F)`[1],
        sum_of_squares = aov_tab$`Sum Sq`[1],
        Df = aov_tab$Df[1],
        f_value = aov_tab$`F value`[1],
        p.value_adjustment = "none",
        p.value_adjusted = aov_tab$`Pr(>F)`[1],
        pairwise_test_statistic = "t.ratio"
      )


      # ------------------------------------------------------
      # Initialize pairwise columns
      # ------------------------------------------------------

      for(pair in pairwise_names) {

        out[[paste0(pair, "_p.value")]] <- NA_real_

        out[[paste0(pair, "_statistic")]] <- NA_real_
      }


      # ------------------------------------------------------
      # EMMEANS
      # ------------------------------------------------------

      emm <- emmeans(
        fit,
        specs = predictor
      )


      # ------------------------------------------------------
      # Pairwise comparisons
      # ------------------------------------------------------

      pw <- pairs(
        emm,
        adjust = "none"
      ) %>%
        as.data.frame()


      # ------------------------------------------------------
      # Normalize contrast labels
      # ------------------------------------------------------

      pw$contrast_norm <- pw$contrast %>%
        str_replace_all("\\s+", "") %>%
        str_replace_all("−", "-")


      # ------------------------------------------------------
      # Fill pairwise columns
      # ------------------------------------------------------

      for(pair in pairwise_names) {

        parts <- str_split(
          pair,
          "_",
          simplify = TRUE
        )

        group1 <- level_map[parts[1]]
        group2 <- level_map[parts[2]]

        target <- paste0(
          group1,
          "-",
          group2
        )

        match_row <- pw %>%
          filter(
            contrast_norm == target
          )

        if(nrow(match_row) == 0) {
          next
        }

        out[[paste0(pair, "_p.value")]] <-
          match_row$p.value[1]

        out[[paste0(pair, "_statistic")]] <-
          match_row$t.ratio[1]
      }


      # ------------------------------------------------------
      # Return one row
      # ------------------------------------------------------

      out
    }
  )

  bind_rows(results)
}

### Fig S1 A- B ####
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

raw_mass_model = lm(mass_final_cm~Phase, data = measures)

raw_len_model = lm(length_final_cm~Phase, data = measures)



### Fig S1 C- F ####

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
dat_testis = subset(data, !is.na(Log10_Testicular_Estimate))

dat_testis$Phase = phase[as.character(dat_testis$Status)]

dat_testis_group = dat_testis %>%
  group_by(Phase) %>%
  summarize(
    mean = mean(Log10_Testicular_Estimate),
    se_testis = sd(Log10_Testicular_Estimate) / sqrt(n())
  )
dat_testis$Log10_Testicular_Estimate = as.numeric(dat_testis$Log10_Testicular_Estimate)

dat_ovary = subset(data, !is.na(Log10_Ovarian_Estimate))
dat_ovary$Log10_Ovarian_Estimate = as.numeric(dat_ovary$Log10_Ovarian_Estimate)

dat_ovary$Phase = phase[as.character(dat_ovary$Status)]

dat_ovary_group = dat_ovary %>%
  group_by(Phase) %>%
  summarize(
    mean = mean(Log10_Ovarian_Estimate),
    se_ovary = sd(Log10_Ovarian_Estimate) / sqrt(n())
  )

mass_model = lm(Change_Mass~Phase, data = dat_experiment)
Length_model = lm(Change_Length~Phase, data = dat_experiment)
testis_model = lm(Log10_Testicular_Estimate ~ Phase, data = dat_testis)
ovary_model = lm(Log10_Ovarian_Estimate ~ Phase, data = dat_ovary)



statistics_table <- populate_statistics(
  mass_model = mass_model,
  Length_model = Length_model,
  testis_model= testis_model, 
  ovary_model = ovary_model,
  raw_mass_model=raw_mass_model ,
  raw_len_model =raw_len_model
)

write.csv(statistics_table, '/Users/ggraham/Desktop/multiome_poa/Manuscript/stat_tables/Fig.S1A-F.csv')
