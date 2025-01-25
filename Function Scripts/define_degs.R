define_degs <- function(data, singular = TRUE) {
  if (!singular) {
    data <- data[data$singular == FALSE, ]
  }
  
  # Assign classes based on conditions
  data$class <- NA  # Initialize class column
  
  data$class[data$d_m_q.value < 0.05 & 
             data$d_m_estimate > 0 & 
             data$d_f_q.value > 0.05] <- 'Early Upregulated'
  
  
  data$class[data$d_m_q.value > 0.05 & 
             data$d_f_q.value < 0.05 & 
             data$d_f_estimate < 0] <- 'Late Upregulated'
  
  data$class[data$d_m_q.value < 0.05 & 
             data$d_m_estimate < 0 & 
             data$d_f_q.value > 0.05] <- 'Early Downregulated'
  
  
  data$class[data$d_m_q.value > 0.05 & 
             data$d_f_q.value < 0.05 & 
             data$d_f_estimate > 0] <- 'Late Downregulated'
  
  data$class[data$d_m_q.value < 0.05 & 
             data$d_f_q.value < 0.05 & 
             data$d_f_estimate > 0 & 
             data$d_m_estimate > 0] <- 'Transiently Upregulated'
  
  data$class[data$d_m_q.value < 0.05 & 
             data$d_f_q.value < 0.05 & 
             data$d_f_estimate < 0 & 
             data$d_m_estimate < 0] <- 'Transiently Downregulated'
  
    data$class[data$d_m_q.value < 0.05 & 
             data$d_f_q.value < 0.05 & 
             data$d_m_estimate < 0 & 
             data$d_f_estimate > 0] <- 'Progressively Downregulated'
    
      data$class[data$d_m_q.value < 0.05 & 
                data$d_f_q.value < 0.05 & 
             data$d_m_estimate > 0 & 
             data$d_f_estimate < 0] <- 'Progressively Upregulated'
      
      data$class[data$f_m_q.value < 0.05 & 
                data$d_f_q.value > 0.05 & 
                data$d_m_q.value > 0.05 & 
             data$f_m_estimate > 0 ] <- 'Terminally Upregulated'
      
  data$class[data$f_m_q.value < 0.05 & 
                data$d_f_q.value > 0.05 & 
                data$d_m_q.value > 0.05 & 
             data$f_m_estimate < 0  ] <- 'Terminally Downregulated'

  return(data)
}

        saveRDS(define_degs, 'Functions/define_degs')
define_degs<- readRDS('Functions/define_degs')
