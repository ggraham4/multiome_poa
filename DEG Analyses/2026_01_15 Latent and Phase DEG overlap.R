library(tidyverse)
library(biomaRt)

clown_go2 = readRDS('Functions/clown_go2')
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



phase_degs = read.csv('~/Desktop/multiome_poa/DEG Outputs/FINAL degs classified w singular.csv')
latent_degs = read.csv('/Users/ggraham/Desktop/multiome_poa/DEG Outputs/2026_01_09 Neg Bin latent_all_clusters.csv')
latent_degs = subset(latent_degs, q.value<0.05)

latent_in_phase =latent_degs$gene[latent_degs$gene %in% phase_degs$gene]
latent_not_in_phase =latent_degs$gene[!latent_degs$gene %in% phase_degs$gene]
latent_degs$presence = ifelse(latent_degs$gene %in% latent_in_phase, 'Phase_Latent', 'Latent')
latent_degs$name = sapply(latent_degs$gene, namer)
latent_degs$direction = ifelse(latent_degs$estimate<0, 'down', 'up')

clown_go2(latent_not_in_phase)%>%dotplot() # lots of ion channels
clown_go2(latent_degs$gene)%>%dotplot() # oh ok its the same

clown_go2(latent_degs$gene[latent_degs$direction=='down'])%>%dotplot()
clown_go2(latent_degs$gene[latent_degs$direction=='up'])%>%dotplot()


phase_in_latent = phase_degs$gene[phase_degs$gene %in% latent_degs$gene]
phase_not_in_latent = phase_degs$gene[!phase_degs$gene %in% latent_degs$gene]

clown_go2(unique(c(phase_degs$gene, latent_degs$gene)))%>%dotplot()
#interesting

table(phase_degs$second_word)
table(phase_degs$cluster)
table(phase_degs$cluster, phase_degs$second_word)



