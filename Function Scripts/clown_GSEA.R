### clown GSEA
library(biomaRt)
library(dplyr)

# 1. Fetch data with SYMBOLS (external_gene_name) to match your geneList
ensembl <- useEnsembl(biomart = "genes", dataset = "aocellaris_gene_ensembl")

biomart_basic <- getBM(
  mart = ensembl, 
  attributes = c("external_gene_name", # Changed from entrezgene_accession
                 "go_id", 
                 "name_1006", 
                 "namespace_1003")) %>%
  filter(namespace_1003 == 'biological_process') %>%
  filter(external_gene_name != "") # Remove entries with missing gene symbols

# 2. Create the mapping dataframes
# TERM2GENE: Column 1 = Term ID (GO), Column 2 = Gene ID (Symbol)
term2gene <- biomart_basic[, c("go_id", "external_gene_name")]

# TERM2NAME: Column 1 = Term ID (GO), Column 2 = Description
term2name <- biomart_basic[, c("go_id", "name_1006")]
# Remove duplicates in term2name so GSEA doesn't get confused
term2name <- unique(term2name)

#saveRDS(term2name, 'Function Scripts/Dependencies/Term2name.rds')

library(clusterProfiler)

  term2gene <- readRDS('Function Scripts/Dependencies/Term2gene.rds')
  term2name<- readRDS( 'Function Scripts/Dependencies/Term2name.rds')

clown_gsea <- function(ranked_list,pval = 0.05) {
    term2gene <- readRDS('Function Scripts/Dependencies/Term2gene.rds')
   term2name<- readRDS( 'Function Scripts/Dependencies/Term2name.rds')
  
  require(clusterProfiler)
  
  # Run GSEA
  # Note: minGSSize and maxGSSize are standard defaults to avoid 
  # analyzing gene sets that are too small (unreliable) or too big (vague).
  res <- GSEA(
    geneList      = ranked_list,    # Your sorted numeric vector
    TERM2GENE     = term2gene,      # The ID -> Symbol mapping
    TERM2NAME     = term2name,      # The ID -> Description mapping
    pvalueCutoff  = pval,
    minGSSize     = 10,             
    maxGSSize     = 500,
    seed          = 123             # Set seed for reproducibility
  )
  
  return(res)
}  
  
saveRDS(clown_gsea,'Functions/clown_gsea.rds')
clown_gsea = readRDS('Functions/clown_gsea.rds')



### example 
library(Seurat)
obj <- readRDS("~/Desktop/nemo.orig_harmony.integration_all_testd_clusters.rds") 
marks_6 = FindMarkers(obj, 6)


sorting_function = function(df, sorting_parameter = 'avg_log2FC', decreasing = TRUE){
    ranks = df[[sorting_parameter]]
    names(ranks) <- rownames(df)
    ranks_sorted <- sort(ranks, decreasing = TRUE)
return(ranks_sorted)

}
sorting_function = readRDS('Functions/gsea_sorter.rds')
clown_gsea = readRDS('Functions/clown_gsea.rds')

gsea_results <- clown_gsea(sorting_function(marks_6))

dotplot(gsea_results)
ridgeplot(gsea_results)

marks_1 = FindMarkers(obj, 1)
clown_gsea(sorting_function(marks_1))%>%dotplot()



