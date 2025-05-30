### making clown_go_2 based on clown go

library(biomaRt)
human_to_ocellaris <- read.csv("C:/Users/Gabe/Desktop/multiome_poa/Reference/hsapiens_to_aocellaris.csv")

human_mart <- useEnsembl(biomart = 'genes',
                         dataset = 'hsapiens_gene_ensembl')

human_go <- getBM(mart = human_mart, 
                  attributes = c('entrezgene_accession', 
                                 'go_id',
                                 'name_1006' 
                  ))

library(dplyr)
joined <- human_to_ocellaris%>% 
  right_join(human_go, by =join_by('hsapiens_name' == 'entrezgene_accession'))%>% 
  subset(!is.na(hsapiens_name)& 
           hsapiens_name != ''&
           !is.na(aocellaris_name)&
           aocellaris_name != '')




term2gene_clown_go2 =joined[, c("go_id", "aocellaris_name")]
term2desc_clown_go2=joined[, c("go_id", "name_1006")]

saveRDS(term2gene_clown_go2, 'Function Scripts/Dependencies/Term2gene_clown_go2.rds')
saveRDS(term2desc_clown_go2, 'Function Scripts/Dependencies/Term2desc_clown_go2.rds')


clown_go2 <- function(significant.list, whole.list = NULL){ #list must be SYMBOL
  
  if (!require(biomaRt)) BiocManager::install('biomaRt')
  if (!require(clusterProfiler)) BiocManager::install('clusterProfiler')
  
  library(biomaRt)
  library(clusterProfiler)
  
  
  
  term2gene <- readRDS('Function Scripts/Dependencies/Term2gene_clown_go2.rds')
  term2desc<- readRDS( 'Function Scripts/Dependencies/Term2desc_clown_go2.rds')
  
  ego2 <- enricher(significant.list, 
                   TERM2GENE=term2gene, 
                   TERM2NAME=term2desc,
                   pAdjustMethod = "fdr", #What p adjust method should I use?
                   universe = whole.list,
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.05)
  return(ego2)
}

saveRDS(clown_go2, 'Functions/clown_go2')
clown_go2<- readRDS('Functions/clown_go2')

library(enrichR)
library(enrichplot)

clown_go2(subset_data$gene[subset_data$cluster==6])%>%dotplot()
