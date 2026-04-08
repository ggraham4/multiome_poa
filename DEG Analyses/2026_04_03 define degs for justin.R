define_degs3 = define_degs_2 = function(data, alpha = 0.05){
  ref_data = read.csv("DEG Analyses/Expression DEGs/pairwise_patterns.csv")
  
  significant = subset(data, av_q.value < alpha)
  if(nrow(significant)<1){return(NULL)}
  significant$full_label = NA
  significant$first_word = NA
  significant$second_word = NA
  significant$short_label = NA

  newd <- data.frame()
  
  for(index in 1:nrow(significant)){
    data_of_interest = significant[index, ]
    
    if(data_of_interest$d_m_p.value < 0.05){
      ref_locs_d_m = which(ref_data$d_m != 'NS')
    } else{
      ref_locs_d_m = which(ref_data$d_m == 'NS')
    }
    
    if(data_of_interest$f_m_p.value < 0.05){
      ref_locs_f_m = which(ref_data$f_m != 'NS')
    } else{
      ref_locs_f_m = which(ref_data$f_m == 'NS')
    }
    
    if(data_of_interest$d_f_p.value < 0.05){
      ref_locs_d_f = which(ref_data$d_f != 'NS')
    } else{
      ref_locs_d_f = which(ref_data$d_f == 'NS')
    }
    
    locs_with_correct_signif_pattern = Reduce(intersect, list(ref_locs_d_f, ref_locs_f_m, ref_locs_d_m))
    
    if(data_of_interest$d_m_p.value < 0.05){
      if(data_of_interest$d_m_estimate < 0){
        dir_locs_d_m = which(ref_data$d_m == '<')
      } else{
        dir_locs_d_m = which(ref_data$d_m == '>')
      }
    } else{
      dir_locs_d_m = which(ref_data$d_m == 'NS')
    }
    
    if(data_of_interest$f_m_p.value < 0.05){
      if(data_of_interest$f_m_estimate < 0){
        dir_locs_f_m = which(ref_data$f_m == '<')
      } else{
        dir_locs_f_m = which(ref_data$f_m == '>')
      }
    } else{
      dir_locs_f_m = which(ref_data$f_m == 'NS')
    }
    
    if(data_of_interest$d_f_p.value < 0.05){
      if(data_of_interest$d_f_estimate < 0){
        dir_locs_d_f = which(ref_data$d_f == '<')
      } else{
        dir_locs_d_f = which(ref_data$d_f == '>')
      }
    } else{
      dir_locs_d_f = which(ref_data$d_f == 'NS')
    }
    
    dir_locs = Reduce(intersect, list(dir_locs_d_f, dir_locs_d_m, dir_locs_f_m))
    
    final_loc = intersect(dir_locs, locs_with_correct_signif_pattern)
    
    if(length(final_loc) > 0){
      full_label = ref_data[final_loc, ]$Full_label
      first_word = ref_data[final_loc, ]$First_word
      second_word = ref_data[final_loc, ]$Second_word
      short_label = ref_data[final_loc, ]$Short_label

      data_of_interest$full_label = full_label
      data_of_interest$first_word = first_word
      data_of_interest$second_word = second_word
      data_of_interest$short_label = short_label

      newd = rbind(newd, data_of_interest)
    }
  }
  return(newd)     
}
