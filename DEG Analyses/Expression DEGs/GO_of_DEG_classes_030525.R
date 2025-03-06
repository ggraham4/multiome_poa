#Negative binomial lower stringency 
{
  library(parallel)

  library(Seurat)
  library(tidyr)
  library(lme4)
  library(dplyr)
  library(MASS)
  library(Signac)
  library('glmGamPoi')
  library(scran)
  library(emmeans)
  library(openxlsx)
  library(ggplot2)
  library(stringr)
  library(forcats)
  library(clusterProfiler)
library(biomaRt)
  library(Polychrome)
  P40 <- createPalette(40, c("#FF0000", "#00FF00", "#0000FF"), range = c(30, 80))
swatch(P40)
names(P40) <- NULL

  mean_expression_cluster_plot<- readRDS('Functions/mean_expression_cluster_plot.rds')
prop_cluster_plot<- readRDS( 'Functions/prop_cluster_plot.rds')
mean_expression_cluster_data<- readRDS('Functions/mean_expression_cluster_data.rds')
clown_go<- readRDS('Functions/clown_go')
define_degs<- readRDS('Functions/define_degs')

}

obj <- readRDS('/Users/ggraham/Desktop/snRNA-seq R Files 122524/RNA Object.rds')


### ANALYSIS ####
together_data <- data.frame()
for(i in 0:31){
  data <- read.csv(paste0('DEG Outputs/012425 Neg Bin w Doms Lower Stringency/cluster_', i, '.csv'))
  #data <- get(paste0('results_cluster',i))
  data$cluster <- i
  together_data <- rbind(together_data, data)
}

together_data_defined <- define_degs(together_data)
together_data_defined_only_signif <- together_data_defined[!is.na(together_data_defined$class),]


### GO of all DEGs ####
clown_go(together_data_defined_only_signif$gene)%>%dotplot()
all_degs_go <- clown_go(together_data_defined_only_signif$gene)
all_degs_go$geneID

enriched_genes_mon_ion_trans <- unlist(strsplit(all_degs_go$geneID[all_degs_go$Description == 'monoatomic ion transport'], "/"))

#name genes
ensembl <- useEnsembl(biomart = "genes", 
                      dataset = "aocellaris_gene_ensembl")

biomart_ocellaris <-
  getBM(
    mart = ensembl, #working mart 
    attributes = c("entrezgene_accession",
                   'entrezgene_description'))

biomart_ocellaris[biomart_ocellaris$entrezgene_accession %in% enriched_genes_mon_ion_trans,]


##### GENES THAT CHANGE BEFORE GONADS CHANGE -- early DEGs
early_degs <- together_data_defined_only_signif$gene[together_data_defined_only_signif$class %in% c('Early Upregulated', "Early Downregulated")]

#nothing
clown_go(early_degs)%>%dotplot()

### genes that change after gonads change -- late DEGs
late_degs <- together_data_defined_only_signif$gene[together_data_defined_only_signif$class %in% c('Late Upregulated', "Late Downregulated")]

clown_go(late_degs)%>%dotplot()

### genes that change as gonads change -- transient DEGs
transient_degs <- together_data_defined_only_signif$gene[together_data_defined_only_signif$class %in% c('Transiently Upregulated', "Transiently Downregulated")]

clown_go(transient_degs)%>%dotplot()

## genes that change at the ends of sex change -- terminal DEGs
terminal_degs <- together_data_defined_only_signif$gene[together_data_defined_only_signif$class %in% c('Terminally Upregulated', "Terminally Downregulated")]
clown_go(terminal_degs)%>%dotplot()

# progressive degs -- are there even any 
progressive_degs <- together_data_defined_only_signif$gene[together_data_defined_only_signif$class %in% c('Progressively Upregulated', "Progressively Downregulated")]

clown_go(progressive_degs)%>%dotplot()


### early up 
early_up_degs <- together_data_defined_only_signif$gene[together_data_defined_only_signif$class %in% c('Early Upregulated')]
clown_go(early_up_degs)%>%dotplot()

### early down
early_down_degs <- together_data_defined_only_signif$gene[together_data_defined_only_signif$class %in% c('Early Downregulated')]
clown_go(early_down_degs)%>%dotplot()

###late up 
late_up_degs <- together_data_defined_only_signif$gene[together_data_defined_only_signif$class %in% c('Late Upregulated')]
clown_go(late_up_degs)%>%dotplot()

### late down
late_down_degs <- together_data_defined_only_signif$gene[together_data_defined_only_signif$class %in% c('Late Downregulated')]
clown_go(late_down_degs)%>%dotplot()

### term up 
term_up_degs <- together_data_defined_only_signif$gene[together_data_defined_only_signif$class %in% c('Terminally Upregulated')]
clown_go(term_up_degs)%>%dotplot()

## term down
term_down_degs <- together_data_defined_only_signif$gene[together_data_defined_only_signif$class %in% c('Terminally Downregulated')]
clown_go(term_down_degs)%>%dotplot()



table(together_data_defined_only_signif$class )
length(together_data_defined_only_signif$class )

all_dpgs <- together_data_defined_only_signif$gene
clown_go(all_dpgs)%>%dotplot()







