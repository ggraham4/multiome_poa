### final defined DEGs fr fr this time
setwd("C:/Users/Gabe/Desktop/multiome_poa")
define_degs2 = readRDS('Functions/define_degs2')
path_prefix = "DEG Outputs/05_11_2025 Neg Bin w Doms New_clustering/"
directory = dir("DEG Outputs/05_11_2025 Neg Bin w Doms New_clustering")
coalesced_data = data.frame()

for(i in directory){
  print(i)
  data = read.csv(paste0(path_prefix, i))
    data = subset(data, av_q.value<0.05 & singular == F)
  data2 = define_degs2(data)
  coalesced_data= rbind(coalesced_data, data2)
  
}

#write.csv(coalesced_data, 'DEG Outputs/FINAL_degs_classified_08_11_2025.csv')

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
  select(-c(X, f_m_q.value, d_m_q.value, d_f_q.value))

write.csv(coalesced_data_2, "DEG Outputs/FINAL_degs_classified_08_11_2025.csv")
