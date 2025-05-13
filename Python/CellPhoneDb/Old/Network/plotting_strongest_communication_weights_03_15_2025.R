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
} ### here, the matrix is organized such that rows are senders and columns are receivers #this is the correct input for igraph directed



matrix_subset <- matrix[]
library(igraph)
library(Polychrome)
set.seed(935234)
P40 <- createPalette(40, c("#FF0000", "#00FF00", "#0000FF"), range = c(30, 80))
swatch(P40)


hist(matrix_subset) #only show top 5% of interactions
signif_interaction_value = max(matrix_subset)*0.79

matrix_subset_threshold <- ifelse(matrix_subset<signif_interaction_value, 0, matrix_subset)
#transpose such that columns are now senders and rows are now receivers

g <- graph_from_adjacency_matrix((matrix_subset_threshold), mode = "directed", weighted = TRUE) 

# Set edge width based on the number of interactions
E(g)$width <- ((E(g)$weight)-min((E(g)$weight)))/((max(E(g)$weight))-min(E(g)$weight))*4

out_data <- read.csv('Python/CellPhoneDb/Network/permutation_testing_senders_receivers_03_15_2025.csv')


Isolated = which(degree(g)==0)
g2 = delete_vertices(g, Isolated)

out_data$cluster <- as.character(out_data$cluster)
cluster_colors <- rep("white", length(V(g2)))
names(cluster_colors) <- V(g2)$name

# Assign colors based on conditions
red_clusters <- as.character(out_data$cluster[!is.na(out_data$sig_rec) & is.na(out_data$sig_send)])
blue_clusters <- as.character(out_data$cluster[is.na(out_data$sig_rec) & !is.na(out_data$sig_send)])
purple_clusters <- as.character(out_data$cluster[!is.na(out_data$sig_rec) & !is.na(out_data$sig_send)])

# Update colors for each group
cluster_colors[names(cluster_colors) %in% red_clusters] <- "orange"
cluster_colors[names(cluster_colors) %in% blue_clusters] <- "purple"
cluster_colors[names(cluster_colors) %in% purple_clusters] <- "gray"

# Assign to vertex colors
V(g2)$color <- cluster_colors

# Plot the graph 
set.seed(4)
plot(g2, 
     edge.color = 'black',
     edge.arrow.size =.5,
    edge.curved = 0.2,
    vertex.label.color='black')

# ok so the plot is right it just doesnt 'Look' right




