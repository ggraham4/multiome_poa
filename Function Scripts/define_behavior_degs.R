define_behavior_degs <- function(data = beh_degs_cluster_0){
  
    data$issignif <- NA
  data$issignif <- ifelse(data$behavior_q.value<0.05|
                            data$behavior_x_status_q.value<0.05,
                          '*',NA)
  
  # Assign classes based on conditions
  data$class <- NA  # Initialize class column
  
  data$class[data$issignif == '*' &
               !is.na(data$issignif) &
               data$behavior_x_status_q.value>0.05&
               data$behavior_estimate<0] <- 'Negatively Correlated'
  
  
  data$class[data$issignif == '*' &
               !is.na(data$issignif) &
               data$behavior_x_status_q.value>0.05&
               data$behavior_estimate>0
               ] <- 'Positively Correlated'
  
  data$class[data$issignif == '*' &
               !is.na(data$issignif) &
               data$behavior_x_status_q.value<0.05
               ] <- 'Interaction'

  return(data)

}
        saveRDS(define_behavior_degs, 'Functions/define_behavior_degs')
define_behavior_degs<- readRDS('Functions/define_behavior_degs')
