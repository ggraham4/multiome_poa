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

matrix_subset <- matrix[]
library(igraph)
library(Polychrome)
set.seed(935234)
P40 <- createPalette(40, c("#FF0000", "#00FF00", "#0000FF"), range = c(30, 80))
swatch(P40)


hist(matrix_subset) #only show top 5% of interactions
signif_interaction_value = max(matrix_subset)*.81

matrix_subset_threshold <- ifelse(matrix_subset<signif_interaction_value, 0, matrix_subset)

g <- graph_from_adjacency_matrix((matrix_subset_threshold), mode = "directed", weighted = TRUE) 

# Set edge width based on the number of interactions
E(g)$width <- ((E(g)$weight)-min((E(g)$weight)))/((max(E(g)$weight))-min(E(g)$weight))*4

V(g)$color <- P40[1:vcount(g)]


for (i in 1:vcount(g)) {
  expected_edges <- sum(matrix_subset[i,] > 0)
  actual_edges <- length(incident(g, i, mode="out"))
  if (expected_edges != actual_edges) {
    print(paste("Node", i, "expected", expected_edges, "edges but has", actual_edges))
  }
}

edge_colors <- c(rep("gray", ecount(g)))  # Default color

for (i in 1:ecount(g)) {
  # Get the endpoints of this edge
  endpoints <- as.numeric(ends(g, i))
  source_node <- endpoints[1]
  
  # Use the source node's color for this edge
  edge_colors[i] <- c(V(g)$color)[source_node+1]
}
edge_colors <- ifelse(is.na(edge_colors), 'darkgray', edge_colors)

# Set the edge colors
#E(g)$color <- edge_colors

Isolated = which(degree(g)==0)
g2 = delete_vertices(g, Isolated)

# Plot the graph 
set.seed(123)
plot(g2, 
     vertex.color ='lightblue', ### should color this based on major senders and receivers
     edge.color = 'black',
     edge.arrow.size =1)





