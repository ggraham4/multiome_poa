#### testing clusters for differential enrichment of DEG types
library(dplyr)
library(ggplot2)
define_degs2 = readRDS('Functions/define_degs2')
### read in data and classify

path = 'DEG Outputs/05_31_2025 Neg Bin Anova First/'
dir = dir('DEG Outputs/05_31_2025 Neg Bin Anova First')

newd = data.frame()
for(i in dir){
  print(i)
  data = read.csv(paste0(path, i))
  classified = define_degs2(data)
  newd = rbind(newd, classified)
}

# here i am only working with non singular degs

### I think for this I use chisq but the interpretation is difficult ###

### or fisher's exact test
fisher_output = data.frame()
first_words <- unique(newd$first_word)
for(i in 0:26){ # for each cluster
  for(word in first_words){
    cluster_data = subset(newd, cluster ==i)
    degs_in_class_in_cluster = sum(cluster_data$first_word ==word)
    degs_not_in_class_in_cluster = sum(cluster_data$first_word !=word)
    
    non_cluster_data = subset(newd, cluster != i)
    degs_in_class_not_in_cluster = sum(non_cluster_data$first_word ==word)
    degs_not_in_class_not_in_cluster = sum(non_cluster_data$first_word !=word)
    
    fisher_matrix = matrix(c(degs_in_class_in_cluster,
                             degs_not_in_class_in_cluster,
                             degs_in_class_not_in_cluster,
                             degs_not_in_class_not_in_cluster),2,2)%>%t()
    
    test = fisher.test(fisher_matrix)
    output_data = data.frame(cluster = i,
                             term = word,
                             p.value = test$p.value, 
                             or = test$estimate)
   fisher_output = rbind(output_data, fisher_output) 
  }
}
fisher_output$issignif =ifelse(fisher_output$p.value<0.05, '*', NA)

### second word ###
second_word_data <- data.frame()
second_words <- unique(newd$second_word) # i think because there are only two terms I dont need to nest this one
for(i in 0:26){
    cluster_data = subset(newd, cluster ==i)
    non_cluster_data = subset(newd, cluster !=i)
    
    cluster_up = sum(cluster_data$second_word=='Upregulated')
    cluster_down = sum(cluster_data$second_word=='Downregulated')
    
    non_cluster_up = sum(non_cluster_data$second_word=='Upregulated')
    non_cluster_down = sum(non_cluster_data$second_word=='Downregulated')
    
    second_word_matrix = matrix(
      c(cluster_up,
        cluster_down,
        non_cluster_up,
        non_cluster_down),2,2)%>%t()
    
    test = fisher.test(second_word_matrix)
    
    output_data = data.frame(cluster = i,
                             p.value = test$p.value, 
                             or = test$estimate)
    
    second_word_data = rbind(second_word_data,output_data )
}
second_word_data$issignif = ifelse(second_word_data$p.value<0.05, "*", NA)

table(newd$second_word)


## what do individual genes look like
mean_expression_cluster_plot = readRDS("Functions/mean_expression_cluster_plot.rds")
obj = readRDS("~/Desktop/optimal_clustering_rna_only.rds")

mean_expression_cluster_plot(obj, 
                             'plp1b',
                             25)

mean_expression_cluster_plot(obj, 
                             'LOC111584347',
                             19)

table(newd$short_label)

table(newd$cluster, newd$first_word)%>%t()%>%heatmap()

table(newd$cluster, newd$first_word)

hclust(dist(table(newd$cluster, newd$first_word)%>%scale()))%>%plot()

hclust(dist(table(newd$cluster, newd$short_label)%>%scale()))%>%plot()

perc = table(newd$cluster, newd$second_word)%>%as.data.frame.matrix()
perc$degs = rowSums(perc)

perc$percent_down = perc$Downregulated/perc$degs
# seems like rg and 6 are mirrors

       
table(newd$second_word)
table(newd$first_word)

