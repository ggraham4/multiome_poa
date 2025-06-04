#### determining which clusters are significantly enriched for DEGs 
#> in the past, I used a chisq test to compare one cluster to all others,
#> i think instead, it makes sense to use a fisher's exact test and ask
#> if a cluster has a significantly larger proporiton of genes analyzed that 
#> were significant

#libraries
library(dplyr)

#set up DEG directory
deg_csvs <- dir('DEG Outputs/05_11_2025 Neg Bin w Doms New_clustering')
deg_path <- 'DEG Outputs/05_11_2025 Neg Bin w Doms New_clustering/'

# read in DEGs list, I think first I will do all singular and nonsingular
all_degs <- data.frame()
for(i in deg_csvs){
  data <- read.csv(paste0(deg_path,i))
  all_degs <- rbind(all_degs, data)
}

# perform fisher's exact test to see if any clusters are significantly different from the rest
# to do this, I will first make a table that summarizes all degs and all total gene in each cluster
#
all_degs$significant <- ifelse(all_degs$av_q.value<0.05, 1,0)

deg_data <- table(all_degs$cluster, all_degs$significant)%>%as.data.frame.matrix()
deg_data$non_degs = deg_data$`0`-  deg_data$`1`
deg_data$degs = deg_data$`1`
deg_data$cluster = rownames(deg_data)
deg_data$cluster<- as.numeric(deg_data$cluster)

#perform loop 
fisher_data_all_degs = data.frame()
for(cluster_of_interest in 0:26){
  cluster_data <- subset(deg_data, cluster_of_interest == cluster)
  other_data <- subset(deg_data, cluster_of_interest != cluster)
  
  fisher_matrix <- matrix(NA,2,2) #here, column 1 is successes (Degs) and column 2 is failures (non degs)
  
  fisher_matrix[1,1] <- cluster_data$degs
  fisher_matrix[1,2] <- cluster_data$non_degs
  
  fisher_matrix[2,1] <- sum(other_data$degs)
  fisher_matrix[2,2] <- sum(other_data$non_degs)
  
  test <- fisher.test(fisher_matrix)
  
  newd <- data.frame(cluster = cluster_of_interest,
                     p_value = test$p.value,
                     estimate = test$estimate)
  
  fisher_data_all_degs <- rbind(fisher_data_all_degs, newd)
  
}

#label signif
fisher_data_all_degs$issignif <- ifelse((fisher_data_all_degs$p_value <0.05 & fisher_data_all_degs$estimate>1),'*', NA)
fisher_data_all_degs$n_degs <- deg_data$degs

#plot the data
library(ggplot2)
ggplot(fisher_data_all_degs, aes(x = cluster, y = n_degs))+
  geom_bar(stat = 'identity')+
  geom_text(label = fisher_data_all_degs$issignif, aes(x = cluster, y = n_degs), size = 10)+
  scale_x_continuous(breaks = 0:26)

### alright now lets do non singular  
some_degs <- subset(all_degs, singular== F)
some_degs$significant <- ifelse(some_degs$av_q.value<0.05, 1,0)

some_deg_data <- table(some_degs$cluster, some_degs$significant)%>%as.data.frame.matrix()
some_deg_data$non_degs = some_deg_data$`0`-  some_deg_data$`1`
some_deg_data$degs = some_deg_data$`1`
some_deg_data$cluster = rownames(some_deg_data)
some_deg_data$cluster<- as.numeric(some_deg_data$cluster)

fisher_data_some_degs = data.frame()
for(cluster_of_interest in 0:26){
  cluster_data <- subset(some_deg_data, cluster_of_interest == cluster)
  other_data <- subset(some_deg_data, cluster_of_interest != cluster)
  
  fisher_matrix <- matrix(NA,2,2) #here, column 1 is successes (Degs) and column 2 is failures (non degs)
  
  fisher_matrix[1,1] <- cluster_data$degs
  fisher_matrix[1,2] <- cluster_data$non_degs
  
  fisher_matrix[2,1] <- sum(other_data$degs)
  fisher_matrix[2,2] <- sum(other_data$non_degs)
  
  test <- fisher.test(fisher_matrix)
  
  newd <- data.frame(cluster = cluster_of_interest,
                     p_value = test$p.value,
                     estimate = test$estimate)
  
  fisher_data_some_degs <- rbind(fisher_data_some_degs, newd)
  
}

#label signif
fisher_data_some_degs$issignif <- ifelse((fisher_data_some_degs$p_value <0.05 & fisher_data_some_degs$estimate>1),'*', NA)
fisher_data_some_degs$n_degs <- some_deg_data$degs

ggplot(fisher_data_some_degs, aes(x = cluster, y = n_degs))+
  geom_bar(stat = 'identity')+
  geom_text(label = fisher_data_some_degs$issignif, aes(x = cluster, y = n_degs), size = 10)+
  scale_x_continuous(breaks = 0:26)



