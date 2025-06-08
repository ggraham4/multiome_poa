### working on a better way to define degs

data = read.csv("DEG Outputs/05_31_2025 Neg Bin Anova First/cluster_6.csv")
data_old = read.csv("DEG Outputs/05_11_2025 Neg Bin w Doms New_clustering/cluster_0.csv")
data_test = read.csv("DEG Outputs/05_11_2025 Neg Bin w Doms New_clustering/cluster_14.csv")

define_degs_2 = function(data, alpha = 0.05, singular = F){
  ref_data = read.csv("DEG Analyses/Expression DEGs/pairwise_patterns.csv")
  
  significant = subset(data, (av_q.value < alpha) & (singular == singular))
  if(nrow(significant)<1){return(NULL)}
  significant$class = NA
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
      class = ref_data[final_loc, ]$Label
      data_of_interest$class = class
      newd = rbind(newd, data_of_interest)
    }
  }
  return(newd)     
}
test_case = define_degs_2(data_old)

saveRDS(define_degs_2, 'Functions/define_degs2')
define_degs_2 = readRDS('Functions/define_degs2')

### testing
path = '~/Desktop/multiome_poa/DEG Outputs/05_11_2025 Neg Bin w Doms New_clustering/'
dir = dir('~/Desktop/multiome_poa/DEG Outputs/05_11_2025 Neg Bin w Doms New_clustering')

dat = data.frame()
for(i in dir){
  print(i)
  data = read.csv(paste0(path, i))
  classed = define_degs_2(data)
  dat = rbind(dat, classed)
  
}

dat

unique(dat$class)
