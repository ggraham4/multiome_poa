library(ggplot2)
library(tidyverse)
library(multcomp)
library(emmeans)


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

      model_terms <- attr(
        terms(model_formula),
        "term.labels"
      )


      # ------------------------------------------------------
      # Identify Phase
      # ------------------------------------------------------

      if(!"Phase" %in% model_terms) {

        stop(
          paste0(
            "Model '",
            model_name,
            "' does not contain Phase as a predictor."
          )
        )
      }


      # ------------------------------------------------------
      # Identify covariates
      #
      # Everything except Phase is treated as a covariate.
      # ------------------------------------------------------

      covariates <- setdiff(
        model_terms,
        "Phase"
      )


      # ------------------------------------------------------
      # Partial ANOVA tests
      #
      # drop1() tests each term while controlling for all
      # other terms in the model.
      # ------------------------------------------------------

      drop_tab <- drop1(
        fit,
        test = "F"
      )


      # ------------------------------------------------------
      # Get Phase statistics
      # ------------------------------------------------------

      phase_row <- which(
        rownames(drop_tab) == "Phase"
      )

      if(length(phase_row) == 0) {

        stop(
          paste0(
            "Could not find a Phase row for model '",
            model_name,
            "'."
          )
        )
      }


      phase_p <- drop_tab$`Pr(>F)`[phase_row]
      phase_ss <- drop_tab$`Sum of Sq`[phase_row]
      phase_df <- drop_tab$Df[phase_row]
      phase_f <- drop_tab$`F value`[phase_row]


      # ------------------------------------------------------
      # Get covariate statistics
      # ------------------------------------------------------

      if(length(covariates) > 0) {

        covariate <- covariates[1]

        cov_row <- which(
          rownames(drop_tab) == covariate
        )

        if(length(cov_row) > 0) {

          cov_p <- drop_tab$`Pr(>F)`[cov_row]
          cov_ss <- drop_tab$`Sum of Sq`[cov_row]
          cov_df <- drop_tab$Df[cov_row]
          cov_f <- drop_tab$`F value`[cov_row]

        } else {

          cov_p <- NA_real_
          cov_ss <- NA_real_
          cov_df <- NA_real_
          cov_f <- NA_real_
        }

      } else {

        covariate <- NA_character_
        cov_p <- NA_real_
        cov_ss <- NA_real_
        cov_df <- NA_real_
        cov_f <- NA_real_
      }


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

        # Phase statistics
        anova_p.value = phase_p,
        sum_of_squares = phase_ss,
        Df = phase_df,
        f_value = phase_f,

        # Covariate statistics
        covariate = covariate,
        covariate_p.value = cov_p,
        covariate_sum_of_squares = cov_ss,
        covariate_Df = cov_df,
        covariate_f_value = cov_f,

        # Pairwise statistics
        p.value_adjustment = "none",
        p.value_adjusted = phase_p,
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
      # EMMEANS for Phase
      #
      # For a covariate model, emmeans automatically evaluates
      # Phase while accounting for the covariate.
      # ------------------------------------------------------

      emm <- emmeans(
        fit,
        specs = ~ Phase
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

### Time ----
tim_mean_sd = measures%>%
  group_by(Phase, Dominance)%>%
  summarize(mean  =mean(Time_Day_2),
            se = sd(Time_Day_2)/sqrt(n()))


### 11KT ----
kt = measures%>%
  group_by(Phase, Dominance)%>%
  subset(!is.na(Log_11KT))%>%
  summarize(mean  =mean(Log_11KT,),
            se = sd(Log_11KT)/sqrt(n()))

kt_model = lm(Log_11KT~Phase, data = measures)
###---- Percent Testicular -----

pct_test = measures%>%
  group_by(Phase, Dominance)%>%
  subset(!is.na(Percent_Testicular))%>%
  summarize(mean  =mean(Percent_Testicular,),
            se = sd(Percent_Testicular)/sqrt(n()))


#### ovarian ----
pct_ov = measures%>%
  group_by(Phase, Dominance)%>%
  subset(!is.na(Percent_Ovarian))%>%
  summarize(mean  =mean(Percent_Ovarian,),
            se = sd(Percent_Ovarian)/sqrt(n()))


### volume ----
measures$Log10_Volume <- as.numeric(measures$Log10_Volume)


###
beh_model = lm(Behaviors_Day_2~Phase, data = measures)
ev_2.5x_model <- lm(Log10_Volume~Phase+length_final_cm, data = measures)
ov_model = lm(Percent_Ovarian~Phase, data = measures)
test_model = lm(Percent_Testicular~Phase, data = measures)
time_model = lm(Time_Day_2~Phase, data = measures)
beh_model = lm(Behaviors_Day_2~Phase, data = measures)

statistics_table <- populate_statistics(
  beh_model = beh_model,
  ev_2.5x_model = ev_2.5x_model,
  ov_model= ov_model, 
  test_model = test_model,
  time_model=time_model ,
  beh_model =beh_model
)

write.csv(statistics_table, '/Users/ggraham/Desktop/multiome_poa/Manuscript/stat_tables/Fig.1B-C.csv')







