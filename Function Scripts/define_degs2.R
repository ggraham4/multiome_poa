### working on a better way to define degs

data = read.csv("DEG Outputs/05_31_2025 Neg Bin Anova First/cluster_6.csv")
data_old = read.csv("DEG Outputs/05_11_2025 Neg Bin w Doms New_clustering/cluster_0.csv")
data_test = read.csv("DEG Outputs/05_11_2025 Neg Bin w Doms New_clustering/cluster_14.csv")

define_degs_2 = function(data, alpha = 0.05, singular_test = F){
  ref_data = read.csv("DEG Analyses/Expression DEGs/pairwise_patterns.csv")
  
  significant = subset(data, av_q.value < alpha & singular == singular_test)
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
test_case = define_degs_2(data_old)

saveRDS(define_degs_2, 'Functions/define_degs2')
define_degs_2 = readRDS('Functions/define_degs2')

### testing
path = '~/Desktop/multiome_poa/DEG Outputs/05_31_2025 Neg Bin Anova First/'
dir = dir('~/Desktop/multiome_poa/DEG Outputs/05_31_2025 Neg Bin Anova First')

dat = data.frame()
for(i in dir){
  print(i)
  data = read.csv(paste0(path, i))
  classed = define_degs_2(data)
  dat = rbind(dat, classed)
  
}

dat

unique(dat$full_label)

dat%>%
  group_by(cluster, short_label)%>%
  summarize(n_degs = n())%>%
  ggplot(aes(x = cluster, y = n_degs, fill = short_label))+
  geom_bar(stat = 'identity',position = 'stack')+
  scale_x_continuous(breaks = c(0:26))

mean_expression_cluster_plot = readRDS("Functions/mean_expression_cluster_plot.rds")
obj = readRDS("~/Desktop/optimal_clustering_rna_only.rds")

mean_expression_cluster_plot(obj, 
                             'LOC111577263',
                             1,
                             'final_clusters')

dat[dat$gene=='LOC111577263',]

mean_expression_cluster_plot(obj, 
                             'esr2b',
                             9,
                             'final_clusters')
dat[dat$gene=='esr2b',]

for(i in unique(dat$short_label)){
  genes = dat[dat$short_label==i,]
  
  first_gene = genes$gene[1]
  cluster = genes$cluster[1]
  
  print(mean_expression_cluster_plot(obj, 
                             first_gene,
                             cluster,
                             'final_clusters')+
    labs(subtitle = i)
)
  
  
  
}

for(i in unique(dat$full_label)){
  genes = dat[dat$full_label==i,]
  
  first_gene = genes$gene[1]
  cluster = genes$cluster[1]
  
  print(mean_expression_cluster_plot(obj, 
                             first_gene,
                             cluster,
                             'final_clusters')+
    labs(subtitle = i)
)
  
  
  
}

dat[dat$gene=='spns2',]
### so this isnt right, so I wonder if I need an estimate threshold

dat[dat$gene=='ehbp1l1a',]

for(i in 1:2){
  
  dap = dat[dat$full_label=='Transiently upregulated, reduced rebound',]
  gene = dap$gene[i]  
  cluster = dap$cluster[i]
  
  print(mean_expression_cluster_plot(obj, 
                             gene,
                             cluster,
                             'final_clusters')
)

  
}
