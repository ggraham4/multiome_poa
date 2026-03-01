### final defined DEGs fr fr this time
#setwd("C:/Users/Gabe/Desktop/multiome_poa")
define_degs2 = readRDS('Functions/define_degs2')
path_prefix = "DEG Outputs/05_11_2025 Neg Bin w Doms New_clustering/"
directory = dir("DEG Outputs/05_11_2025 Neg Bin w Doms New_clustering")
coalesced_data = data.frame()


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


for(i in directory){
  print(i)
  data = read.csv(paste0(path_prefix, i))
    data = subset(data, av_q.value<0.05)# edit here
  data2 = define_degs3(data)
  coalesced_data= rbind(coalesced_data, data2)
  
}

write.csv(coalesced_data, 'DEG Outputs/FINAL degs classified w singular.csv')

library(biomaRt)
ensembl <- useEnsembl(biomart = "genes", 
                      dataset = "aocellaris_gene_ensembl")

a <- listAttributes(ensembl)
biomart_basic <-
  getBM(
    mart = ensembl, #working mart 
    attributes = c("entrezgene_accession",
                   'entrezgene_description'))


namer = function(gene){
  gene_name = biomart_basic$entrezgene_description[biomart_basic$entrezgene_accession==gene]
  return(gene_name[1])
  
}


coalesced_data$gene_name= sapply(FUN= namer, X=coalesced_data$gene)
library(tidyverse)
coalesced_data_2  =coalesced_data%>%
  relocate(gene_name, .after = gene)%>%
  dplyr::select(-c(X, f_m_q.value, d_m_q.value, d_f_q.value))%>%
  distinct()

#write.csv(coalesced_data_2, 'DEG Outputs/FINAL degs classified w singular.csv')

