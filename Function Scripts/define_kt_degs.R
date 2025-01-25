define_kt_degs <- function(data ){
  
    data$issignif <- NA
  data$issignif <- ifelse(data$interaction_summary_q.value<0.05|
                            data$interaction_anova_q.value<0.05|
                            data$KT_summary_q.value<0.05|
                            data$KT_anova_q.value<0.05,
                          '*',NA)
  
  # Assign classes based on conditions
  data$class <- NA  # Initialize class column
  
  data$class[data$issignif == '*' &
               !is.na(data$issignif) &
              (data$interaction_anova_q.value>0.05|
                  data$interaction_summary_q.value>0.05)&
               data$KT_estimate<0] <- 'Negatively Correlated'
  
  
  data$class[data$issignif == '*' &
               !is.na(data$issignif) &
              (data$interaction_anova_q.value>0.05|
                  data$interaction_summary_q.value>0.05)&
               data$KT_estimate>0] <- 'Positively Correlated'
  
  data$class[data$issignif == '*' &
               !is.na(data$issignif) &
               (data$interaction_anova_q.value<0.05|
               data$interaction_summary_q.value<0.05)
               ] <- 'Interaction'

  return(data)

}
data <- define_kt_degs(data)

        saveRDS(define_kt_degs, 'Functions/define_kt_degs')
define_kt_degs<- readRDS('Functions/define_kt_degs')
