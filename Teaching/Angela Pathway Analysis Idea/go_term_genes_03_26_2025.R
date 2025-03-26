library(biomaRt)
library(Seurat)
library(tidyverse)

### read in object
obj <- readRDS("C:/Users/Gabe/Desktop/RNA Object.rds")

all_genes = rownames(obj@assays$RNA$data)

#make a biomart
ocellaris_mart <- useEnsembl(biomart = 'genes',
                             dataset = 'aocellaris_gene_ensembl')

att = attributes(ocellaris_mart)$attributes

ocellaris_go <- getBM(mart = ocellaris_mart,
                                attributes = c('entrezgene_accession',
                                               'go_id',
                                               'name_1006'
                                               )) # i guess biomart doesnt have kegg for a ocellaris

enzymes_of_interest <- c(
  'hsd3b1',
  'hsd17b14',
  'hsd17b1',
  'LOC111577263' #brain aromatase
)
pathways_of_interest <- c(
  'steroid biosynthetic process',
  'glucocorticoid biosynthetic process',
  'hormone biosynthetic process',
  'steroid hydroxylase activity'
)

ocellaris_go$entrezgene_accession[ocellaris_go$name_1006%in%pathways_of_interest]
# it might be better to do this with human genes cause those are better annotated

## alternatively, we could use zack's IEG strat to find genes coexpressed


