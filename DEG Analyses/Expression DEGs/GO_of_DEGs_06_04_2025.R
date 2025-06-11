#GO of DEGs in enriched clusters
## here, I am going to use clown_go2 to test if the degs in a given cluster are
# enriched for anything

library(dplyr)
library(clusterProfiler)
library(ggplot2)
library(biomaRt)
clown_go2 <- readRDS("Functions/clown_go2")

### read in data
deg_csvs <- dir('DEG Outputs/05_11_2025 Neg Bin w Doms New_clustering')
deg_path <- 'DEG Outputs/05_11_2025 Neg Bin w Doms New_clustering/'

# read in DEGs list, I think first I will do all singular and nonsingular
all_degs <- data.frame()
for(i in deg_csvs){
  data <- read.csv(paste0(deg_path,i))
  all_degs <- rbind(all_degs, data)
}

all_degs_subset <- subset(all_degs, av_q.value <0.05)

### cluster 6 
#singular and non singular
clown_go2(all_degs_subset$gene[all_degs_subset$cluster==6])%>%dotplot()

#non-singular
clown_go2(all_degs_subset$gene[all_degs_subset$cluster==6& all_degs_subset$singular==F])%>%dotplot()
#ESTROGEN RESPONSE ELEMENT BINDING BABEY

### cluster 1
#singular and non singular
#clown_go2(all_degs_subset$gene[all_degs_subset$cluster==1])%>%dotplot()
#NA

#non-singular
#clown_go2(all_degs_subset$gene[all_degs_subset$cluster==1& all_degs_subset$singular==F])%>%dotplot()
#still NA

### cluster 11
#singular and non singular
#clown_go2(all_degs_subset$gene[all_degs_subset$cluster==11])%>%dotplot()
#NA

#non-singular
clown_go2(all_degs_subset$gene[all_degs_subset$cluster==11& all_degs_subset$singular==F])%>%dotplot()
#humoral immune response

### all degs
clown_go2(all_degs_subset$gene)%>%dotplot()
clown_go2(all_degs_subset$gene[all_degs_subset$singular==F])%>%dotplot()

#what are the genes
#clown_go2(all_degs_subset$gene[all_degs_subset$singular==F])
# the LOC is another PGR

clown_go2(all_degs_subset$gene[all_degs_subset$cluster==6& all_degs_subset$singular==F])
#pgr and myosin

### repeated DEGs
repeated_genes <- table(all_degs_subset$gene)%>%as.data.frame()
