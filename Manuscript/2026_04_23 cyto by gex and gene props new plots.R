#libs
library(ggplot2)
library(ggsignif)
library(tidyverse)
library(emmeans)
library(Seurat)
library(CytoTRACE)

# read in data, subset to 6_NPO
obj  = readRDS("~/Desktop/optimal_clustering_rna_only.rds")
sub_6 = subset(obj,res0.8_50nn_40PC_45LSI==6)
degs = read.csv('/Users/ggraham/Desktop/multiome_poa/DEG Outputs/FINAL degs classified w singular.csv')
degs_6 = subset(degs, cluster == 6)


# define genes to examine
genes_interest = c(
  'drd3',
  'tacr3a',
  'npy7r',
  'cckb',
  'LOC111571064' # gnrh1
)

#loop through them to test for proportion differences
pep = readRDS('Functions/prop_expressing_plot.rds')
for(i in genes_interest){
  p = pep(sub_6, i)
  
    ggsave(plot = p,
       file = paste0(i, '_prop_6.svg'),
       device = "svg",
       units = "in",
       width = 2,
       height = 2,
       path = paste0("Manuscript/Plots/Manuscript v1.2/prop_expressing_6_NPO"))

  
}

sub_6$cyto = CytoTRACE(sub_6@assays$RNA$data%>%as.matrix())$CytoTRACE
# loop throgh genes interest and test for cyto differences
cdp = readRDS('Functions/cyto_differences_plot.rds')
for(i in genes_interest){
  c = cdp(sub_6, i)
  
    ggsave(plot = c,
       file = paste0(i, '_cyto_6.svg'),
       device = "svg",
       units = "in",
       width = 3,
       height = 2,
       path = paste0("Manuscript/Plots/Manuscript v1.2/cyto_genes_6_NPO"))

}



