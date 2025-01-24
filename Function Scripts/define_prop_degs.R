        define_degs_prop <- function(data, singular = TRUE) {
  if (!singular) {
    data <- data[data$singular == FALSE, ]
  }
  
  # Assign classes based on conditions
  data$class <- NA  # Initialize class column
  
  data$class[data$d_m_q.value < 0.05 & 
             data$prop_expressing_D >data$prop_expressing_M  & 
             data$d_f_q.value > 0.05] <- 'Early Upregulated'
  
  
  data$class[data$d_m_q.value > 0.05 & 
             data$d_f_q.value < 0.05 & 
             data$prop_expressing_D <data$prop_expressing_F] <- 'Late Upregulated'
  
  data$class[data$d_m_q.value < 0.05 & 
               data$prop_expressing_D <data$prop_expressing_M & 
             data$d_f_q.value > 0.05] <- 'Early Downregulated'
  
  
  data$class[data$d_m_q.value > 0.05 & 
             data$d_f_q.value < 0.05 & 
               data$prop_expressing_D >data$prop_expressing_F] <- 'Late Downregulated'
  
  data$class[data$d_m_q.value < 0.05 & 
             data$d_f_q.value < 0.05 & 
            data$prop_expressing_D >data$prop_expressing_F& 
            data$prop_expressing_D >data$prop_expressing_M] <- 'Transiently Upregulated'
  
  data$class[data$d_m_q.value < 0.05 & 
             data$d_f_q.value < 0.05 & 
             data$prop_expressing_D <data$prop_expressing_F & 
             data$prop_expressing_D <data$prop_expressing_M] <- 'Transiently Downregulated'
  
    data$class[data$d_m_q.value < 0.05 & 
             data$d_f_q.value < 0.05 & 
             data$prop_expressing_D <data$prop_expressing_M & 
             data$prop_expressing_D >data$prop_expressing_F] <- 'Progressively Downregulated'
    
      data$class[data$d_m_q.value < 0.05 & 
                data$d_f_q.value < 0.05 & 
            data$prop_expressing_D >data$prop_expressing_M & 
              data$prop_expressing_D <data$prop_expressing_F] <- 'Progressively Upregulated'
      
      data$class[data$f_m_q.value < 0.05 & 
                data$d_f_q.value > 0.05 & 
                data$d_m_q.value > 0.05 & 
              data$prop_expressing_F >data$prop_expressing_M] <- 'Terminally Upregulated'
      
  data$class[data$f_m_q.value < 0.05 & 
                data$d_f_q.value > 0.05 & 
                data$d_m_q.value > 0.05 & 
             data$prop_expressing_F <data$prop_expressing_M ] <- 'Terminally Downregulated'
  

  data$issignif <- NA
  data$issignif <- ifelse(data$f_m_q.value<0.05|
                            data$d_m_q.value<0.05|
                            data$d_f_q.value<0.05,
                          '*',NA)

  return(data)
        }

        saveRDS(define_degs_prop, 'Functions/define_degs_prop.rds')
define_degs_prop<- readRDS('Functions/define_degs_prop.rds')
