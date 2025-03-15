### Making a 3D Graph of Cell surface receptor expression to make spatial inference
library(igraph)
library(ggplot2)
library(MASS) # For MDS implementation
library(dplyr)

#read in cell cell communication matrix
example_data = read_tsv("A:/CellPhoneDB 030225/degs_analysis_relevant_interactions_03_03_2025_224407.txt")

example_data_filtered = example_data%>%
  dplyr::select(-c(1, 3:13))

example_data_filtered_long = example_data_filtered%>%
  pivot_longer(cols = colnames(example_data_filtered[,-1]), values_to = 'signif', names_to = 'cluster_pair')%>%
  subset(signif == 1)

example_data_filtered_long <- example_data_filtered_long %>%
  mutate(
    sender = sapply(strsplit(cluster_pair, "\\|"), `[`, 1),
    receiver = sapply(strsplit(cluster_pair, "\\|"), `[`, 2)
  )%>%
  subset(sender != '30' & sender !='15' & receiver != '30' & receiver!= '15')


example_data_filtered_long_grouped <- example_data_filtered_long%>%
  group_by(sender, receiver)%>%
  summarize(interactions =n())

#make interaction matrix
matrix = matrix(NA, 30, 30)
colnames(matrix) = c(0:14, 16:29,31)
rownames(matrix) = c(0:14, 16:29,31)


shortened_name = example_data_filtered_long_grouped

for(i in colnames(matrix)){
  for(j in colnames(matrix)){
    print(paste0(i,"|",j))
    sender = which(colnames(matrix)==i)
    reciever = which(rownames(matrix)==j)
    
    if(length(shortened_name$interactions[shortened_name$sender == i & shortened_name$receiver ==j])==0){
      matrix[i,j]=0
    }
    else{
      matrix[i,j] = shortened_name$interactions[shortened_name$sender ==i & shortened_name$receiver ==j]
    }
    
  }
} ### here, the matrix is organized such that rows are senders and columns are receivers

### doing stats

### I think I will do a permutation test to get the exact p value#
#I think I want to transpose the matrix such that columns are now senders
transposed_mat <- t(matrix)

hist(colSums(transposed_mat))

### I want to make it such that cells are shuffled around to various columns 

sum_vals = sum(transposed_mat)

distribution_data = c()
for(i in 1:100000){
  print(i)
  sampled_values <-sample(matrix)
  fake_matrix <- matrix(sampled_values, 30, 30)
  
  prop_perm <- sum(fake_matrix[,1])/sum_vals
  distribution_data <- append(distribution_data, prop_perm)
  
}

hist(distribution_data)
## incredibly normally distributed

# calculate the actual mean percent from each cluster
out_data <- data.frame()
for(i in colnames(transposed_mat)){
  cluster_index = which(colnames(transposed_mat)==i)
  prop_real_send <- sum(transposed_mat[,cluster_index])/sum(transposed_mat)
  
  newd = data.frame(cluster = i,
                    prop_real_send = prop_real_send)
  out_data <- rbind(newd, out_data)
}

### next calculate a null distribution

p_value = sum(distribution_data>prop_real)/length(distribution_data)

perm_p_sends<- c()
for(i in out_data$cluster){

  perm_p <- sum(distribution_data>out_data$prop_real[out_data$cluster==i])/length(distribution_data)

  perm_p_sends <- append(perm_p_sends, perm_p)
}

out_data$perm_p_send <- perm_p_sends

#im pretty sure I can use the same distribution for receptions
out_data_rec <- data.frame()
for(i in rownames(transposed_mat)){
  cluster_index = which(rownames(transposed_mat)==i)
  prop_real_rec <- sum(transposed_mat[cluster_index,])/sum(transposed_mat)
  
  newd = data.frame(cluster = i,
                    prop_real_rec = prop_real_rec)
  out_data_rec <- rbind(newd, out_data_rec)
}

out_data$prop_real_rec <- out_data_rec$prop_real_rec

perm_p_recs<- c()
for(i in out_data$cluster){
  
  perm_p <- sum(distribution_data>out_data$prop_real_rec[out_data$cluster==i])/length(distribution_data)
  
  perm_p_recs <- append(perm_p_recs, perm_p)
}
out_data$perm_p_rec <- perm_p_recs


out_data$sig_rec <- ifelse(out_data$perm_p_rec<0.05, '*', NA)
out_data$sig_send <- ifelse(out_data$perm_p_send<0.05, '*', NA)

#write.csv(out_data, 'Python/CellPhoneDb/Network/permutation_testing_senders_receivers_03_15_2025.csv')

### testing if a cluster_cluster pair is significant - do that later...

