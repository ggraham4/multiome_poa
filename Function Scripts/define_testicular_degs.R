define_testicular_degs <- function(data = Test_degs_cluster_0){
  
    data$issignif <- NA
  data$issignif <- ifelse(data$interaction_summary_q.value<0.05|
                            data$interaction_anova_q.value<0.05|
                            data$percent_testicular_summary_q.value<0.05|
                            data$percent_testicular_anova_q.value<0.05,
                          '*',NA)
  
  # Assign classes based on conditions
  data$class <- NA  # Initialize class column
  
  data$class[data$issignif == '*' &
               !is.na(data$issignif) &
              (data$interaction_anova_q.value>0.05|
                  data$interaction_summary_q.value>0.05)&
               data$percent_testicular_estimate<0] <- 'Negatively Correlated'
  
  
  data$class[data$issignif == '*' &
               !is.na(data$issignif) &
              (data$interaction_anova_q.value>0.05|
                  data$interaction_summary_q.value>0.05)&
               data$percent_testicular_estimate>0] <- 'Positively Correlated'
  
  data$class[data$issignif == '*' &
               !is.na(data$issignif) &
               (data$interaction_anova_q.value<0.05|
               data$interaction_summary_q.value<0.05)
               ] <- 'Interaction'

  return(data)

}
data <- define_testicular_degs(Test_degs_cluster_0)

        saveRDS(define_testicular_degs, 'Functions/define_testicular_degs')
define_testicular_degs<- readRDS('Functions/define_testicular_degs')
