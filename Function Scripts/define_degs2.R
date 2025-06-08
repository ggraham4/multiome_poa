define_degs2 <- function(data, singular = TRUE) {
  if (!singular) {
    data <- data[data$singular == FALSE, ]
  }
  
  # Assign classes based on conditions
  data$class <- NA  # Initialize class column
  
  data$class[data$d_m_p.value < 0.05 & data$av_q.value<0.05 &
               data$d_m_estimate > 0 & 
               data$d_f_p.value > 0.05] <- 'Early Upregulated'
  
  
  data$class[data$d_m_p.value > 0.05 & data$av_q.value<0.05 &
               data$d_f_p.value < 0.05 & 
               data$d_f_estimate < 0] <- 'Late Upregulated'
  
  data$class[data$d_m_p.value < 0.05 & data$av_q.value<0.05 &
               data$d_m_estimate < 0 & 
               data$d_f_p.value > 0.05] <- 'Early Downregulated'
  
  
  data$class[data$d_m_p.value > 0.05 & data$av_q.value<0.05 &
               data$d_f_p.value < 0.05 & 
               data$d_f_estimate > 0] <- 'Late Downregulated'
  
  data$class[data$d_m_p.value < 0.05 & data$av_q.value<0.05 &
               data$d_f_p.value < 0.05 & 
               data$d_f_estimate > 0 & 
               data$d_m_estimate > 0] <- 'Transiently Upregulated'
  
  data$class[data$d_m_p.value < 0.05 & data$av_q.value<0.05 &
               data$d_f_p.value < 0.05 & 
               data$d_f_estimate < 0 & 
               data$d_m_estimate < 0] <- 'Transiently Downregulated'
  
  data$class[data$d_m_p.value < 0.05 & data$av_q.value<0.05 &
               data$d_f_p.value < 0.05 & 
               data$d_m_estimate < 0 & 
               data$d_f_estimate > 0] <- 'Progressively Downregulated'
  
  data$class[data$d_m_p.value < 0.05 & data$av_q.value<0.05 &
               data$d_f_p.value < 0.05 & 
               data$d_m_estimate > 0 & 
               data$d_f_estimate < 0] <- 'Progressively Upregulated'
  
  data$class[data$f_m_p.value < 0.05 & data$av_q.value<0.05 &
               data$d_f_p.value > 0.05 & 
               data$d_m_p.value > 0.05 & 
               data$f_m_estimate > 0 ] <- 'Terminally Upregulated'
  
  data$class[data$f_m_p.value < 0.05 & data$av_q.value<0.05 &
               data$d_f_p.value > 0.05 & 
               data$d_m_p.value > 0.05 & 
               data$f_m_estimate < 0  ] <- 'Terminally Downregulated'
  
  data$class = ifelse(is.na(data$class) & (data$av_q.value<0.05), 'No Tukey', data$class)
  
  return(data)
}

saveRDS(define_degs2, 'Functions/define_degs2')
define_degs2<- readRDS('Functions/define_degs2')

### i think this is not very good, maybe there is a better way to do it?

