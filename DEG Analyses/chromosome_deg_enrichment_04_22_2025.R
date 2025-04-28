library(edgeR)
library(Seurat)
library(tidyverse)
library(biomaRt)
define_degs<- readRDS('Functions/define_degs')
clown_go<- readRDS('Functions/clown_go')


obj <- readRDS('/Users/ggraham/Desktop/snRNA-seq R Files 122524/RNA Object.rds')

clown_mart <- useEnsembl(biomart = 'genes',
                             dataset = 'aocellaris_gene_ensembl')

clown_data <- getBM(mart = clown_mart, #call the human mart
                      attributes = c('entrezgene_accession', #entrezgene_acession is the NIH gene name
                                     'chromosome_name'
                      ))


## read in DEGs list
together_data <- data.frame()
for(i in 0:31){
  if(i %in% c(15,30)){next}
  data <- read.csv(paste0('DEG Outputs/012425 Neg Bin w Doms Lower Stringency/cluster_', i, '.csv'))
  #data <- get(paste0('results_cluster',i))
  data$cluster <- i
  together_data <- rbind(together_data, data)
}

together_data_defined <- define_degs(together_data)
together_data_defined_only_signif <- together_data_defined[!is.na(together_data_defined$class),]

all_degs <- unique(together_data_defined_only_signif$gene)


#### find which chromosomes are enriched
clown_data_subset <- subset(clown_data, entrezgene_accession %in% all_degs & chromosome_name != 'JAJUWX010000203.1')

ggplot(clown_data_subset, aes(x = as.factor(as.numeric(chromosome_name))))+
  geom_bar()

## perform fisher's exact test to compare number of genes present in a cluster vs # of degs

clown_data_subset_table <- clown_data%>%
  group_by(chromosome_name)%>%
  summarize(n_genes = n())

clown_data_deg_table <- clown_data_subset%>%
  group_by(chromosome_name)%>%
  summarize(n_degs = n())

clown_joined <- clown_data_subset_table%>%
  right_join(clown_data_deg_table, by = 'chromosome_name')

fisher_data <- data.frame()
for(i in 1:24){
  chromosome_data <- subset(clown_joined, chromosome_name == i)
  nonchromosome_data <- subset(clown_joined, chromosome_name != i)
  
  chromosome_genes <- sum(chromosome_data$n_genes - chromosome_data$n_degs)
  chromosome_degs <- sum(chromosome_data$n_degs)
  
  nonchromosome_genes <-  sum(nonchromosome_data$n_genes - nonchromosome_data$n_degs)
  nonchromosome_degs <-  sum(nonchromosome_data$n_degs)

  fisher_matrix <- matrix(NA, 2, 2)
  fisher_matrix[1,1] <- chromosome_degs
  fisher_matrix[1,2] <- chromosome_genes

  fisher_matrix[2,2] <- nonchromosome_genes
  fisher_matrix[2,1] <- nonchromosome_degs

  test <- fisher.test(fisher_matrix)
  
  newd <- data.frame(chromosome = i,
                     test_statistic = unname(test$estimate),
                     CI_low = test$conf.int[1],
                     CI_high = test$conf.int[2],
                     p_value = test$p.value)
  fisher_data <- rbind(fisher_data, newd)
  
}


### are they enriched for types of DEGs
#how do I ask this efficiently? I think fisher test woud be too clunky
# split into halves

split_names <- str_split(together_data_defined_only_signif$class, ' ')

together_data_defined_only_signif$first_word <- sapply(split_names,"[[",1)
together_data_defined_only_signif$second_word<- sapply(split_names,"[[",2)

together_data_defined_only_signif_chrm <- together_data_defined_only_signif%>%
  right_join(clown_data, by = join_by('gene'=='entrezgene_accession'))%>%
  subset(!is.na(class))

second_word_plot <- together_data_defined_only_signif_chrm%>%
  group_by(chromosome_name, second_word)%>%
  summarize(n = n())%>%
  subset(chromosome_name!= 'JAJUWX010000203.1')

ggplot(second_word_plot, aes(x = as.factor(as.numeric(chromosome_name)),y = n, fill = second_word))+
  geom_bar(stat = 'identity', position = 'fill')+
  geom_hline(yintercept = 0.5)


first_word_plot <- together_data_defined_only_signif_chrm%>%
  group_by(chromosome_name, first_word)%>%
  summarize(n = n())%>%
  subset(chromosome_name!= 'JAJUWX010000203.1')

ggplot(first_word_plot, aes(x = as.factor(as.numeric(chromosome_name)),y = n, fill = first_word))+
  geom_bar(stat = 'identity', position = 'fill')

### i think chisq is what I should use here

cont_table_first <- table(together_data_defined_only_signif_chrm$chromosome_name, together_data_defined_only_signif_chrm$first_word)

chisq.test(cont_table_first)

cont_table_second<- table(together_data_defined_only_signif_chrm$chromosome_name, together_data_defined_only_signif_chrm$second_word)

chisq.test(cont_table_second)
