library(Seurat)
library(ggplot2)
library(tidyverse)
latent_plotter = function(gene, cluster){
  gex =data.frame(exp =  obj@assays$RNA$data[gene,][obj@meta.data$res0.8_50nn_40PC_45LSI==cluster])
  gex$SexState = obj@meta.data$SexState[obj@meta.data$res0.8_50nn_40PC_45LSI==cluster]
  gex$Status = obj@meta.data$Status[obj@meta.data$res0.8_50nn_40PC_45LSI==cluster]
  gex$individual = obj@meta.data$individual[obj@meta.data$res0.8_50nn_40PC_45LSI==cluster]

  gex_grouped = gex%>%
    group_by(individual, Status)%>%
    summarize(mean_gex = mean(exp),
              SexState =mean(SexState))
  
 p= ggplot(gex_grouped, aes(x = SexState, y = mean_gex))+
   geom_point(aes(color = Status))+
   geom_smooth()
 return(p)
  
}


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
library(clusterProfiler)
clown_go = readRDS('Functions/clown_go.rds')
clown_go2 = readRDS('Functions/clown_go2')

########
# get the data
data = read.csv('/Users/ggraham/Desktop/multiome_poa/DEG Outputs/2026_01_09 Neg Bin latent_all_clusters.csv')


sub_data = subset(data, q.value < 0.05)
sub_data$direction = ifelse(sub_data$estimate>0, 'Up','Down')

sub_data$name = sapply(sub_data$gene, namer)

sub_data2 = subset(sub_data, q.value < 0.05 & cluster == 2)
sub_data1 = subset(sub_data, q.value < 0.05 & cluster == 1)
sub_data6 = subset(sub_data, q.value < 0.05 & cluster == 6)
sub_data11 = subset(sub_data, q.value < 0.05 & cluster == 11)

clown_go2(sub_data6$gene, p =0.1)%>%dotplot()
clown_go2(sub_data6$gene[sub_data$direction=='Up'])%>%dotplot()

clown_go2(sub_data6$gene[sub_data$direction=='Up'])$geneID[clown_go2(sub_data6$gene[sub_data$direction=='Up'])$Description == 'cellular response to estrogen stimulus']
clown_go2(sub_data6$gene[sub_data$direction=='Up'])$geneID[clown_go2(sub_data6$gene[sub_data$direction=='Up'])$Description == 'cellular response to testosterone stimulus']
#LOC111568069 pgr like

clown_go2(sub_data6$gene[sub_data$direction=='Down'], p = 0.1)%>%dotplot()


clown_go2(sub_data1$gene)%>%dotplot()
clown_go2(sub_data1$gene)$geneID[clown_go2(sub_data1$gene)$Description == 'neuroendocrine cell differentiation']
# nkx2.2

clown_go2(sub_data1$gene)$geneID[clown_go2(sub_data1$gene)$Description == 'neural crest cell migration involved in autonomic nervous system development']
# sema3ab

clown_go2(sub_data2$gene, p = 0.1)%>%dotplot()
clown_go2(sub_data11$gene, p = 0.1)%>%dotplot()


### by cluster

t = sub_data%>%group_by(cluster, direction)%>%
  summarize(n_genes = n())

ggplot(t, aes(x = as.factor(cluster), y = n_genes, fill = direction))+
  geom_bar(stat = 'identity', position = 'stack')+
  scale_y_continuous(breaks = 0:25)


obj <- readRDS('~/Desktop/optimal_clustering_rna_only.rds')

measures = read.csv('2025_12_26 all_data.csv')
measures$individual = measures$Fish
measures_sexstate =dplyr::select(measures,   c('individual', 'SexState'))

obj@meta.data = obj@meta.data%>%
  left_join(measures_sexstate, by = 'individual')

# plotting ###
latent_plotter('nkx2.2a', 1)

