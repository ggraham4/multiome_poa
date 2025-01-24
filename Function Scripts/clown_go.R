#anemonefish_go

ensembl <- useEnsembl(biomart = "genes", 
                      dataset = "aocellaris_gene_ensembl")

a <- listAttributes(ensembl)
biomart_basic <-
  getBM(
    mart = ensembl, #working mart 
    attributes = c("entrezgene_accession",
                   'go_id',
                   'name_1006',
                   'namespace_1003'))%>%
  filter(namespace_1003 == 'biological_process') #change this to set which domain you care about

term2gene =biomart_basic[, c("go_id", "entrezgene_accession")]
term2desc=biomart_basic[, c("go_id", "name_1006")]

saveRDS(term2gene, '/volumes/jrhodes/Fish Lab/Experiments/sex change single nuc POA/R/Gabe/Term2gene.rds')
saveRDS(term2desc, '/volumes/jrhodes/Fish Lab/Experiments/sex change single nuc POA/R/Gabe/Term2desc.rds')

clown_go <- function(significant.list, whole.list = NULL){ #list must be SYMBOL
  term2gene <- readRDS('/volumes/jrhodes/Fish Lab/Experiments/sex change single nuc POA/R/Gabe/Term2gene.rds')
  term2desc<- readRDS( '/volumes/jrhodes/Fish Lab/Experiments/sex change single nuc POA/R/Gabe/Term2desc.rds')
  
  ego2 <- enricher(significant.list, 
                   TERM2GENE=term2gene, 
                   TERM2NAME=term2desc,
                   pAdjustMethod = "fdr", #What p adjust method should I use?
                   universe = whole.list,
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.05)
  return(ego2)
}

        saveRDS(clown_go, 'Functions/clown_go')
clown_go<- readRDS('Functions/clown_go')
  
  