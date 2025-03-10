### Making a 3D Graph of Cell surface receptor expression to make spatial inference

library(igraph)
library(ggplot2)
library(MASS) # For MDS implementation
library(dplyr)

#read in cell cell communication matrix
example_data = read_tsv("A:/CellPhoneDB 030225/degs_analysis_relevant_interactions_03_03_2025_224407.txt")
unique(example_data$directionality)

example_data_adhesion = subset(example_data, directionality %in% c('Adhesion-Adhesion', 'Gap-Gap'))
example_data_adhesion_filtered = example_data_adhesion%>%
  dplyr::select(-c(1, 3:13))
  
example_data_adhesion_filtered_long = example_data_adhesion_filtered%>%
  pivot_longer(cols = colnames(example_data_adhesion_filtered[,-1]), values_to = 'signif', names_to = 'cluster_pair')%>%
  subset(signif == 1)

#make interaction matrix
matrix = matrix(NA, 30, 30)
colnames(matrix) = c(0:14, 16:29,31)
rownames(matrix) = c(0:14, 16:29,31)



#need to parse example_data_adhesion_filtered_long senders and recievers
example_data_adhesion_filtered_long_parsed <- example_data_adhesion_filtered_long %>%
  mutate(
    sender = sapply(strsplit(cluster_pair, "\\|"), `[`, 1),
    receiver = sapply(strsplit(cluster_pair, "\\|"), `[`, 2)
  )
example_data_adhesion_filtered_long_parsed_grouped <- example_data_adhesion_filtered_long_parsed%>%
  group_by(sender, receiver)%>%
  summarize(interactions =n())

shortened_name = example_data_adhesion_filtered_long_parsed_grouped

for(i in 0:30){
  for(j in 0:30){
    print(paste0(i,"|",j))
    sender = which(colnames(matrix)==i)
    reciever = which(rownames(matrix)==j)
    
    if(length(shortened_name$interactions[shortened_name$sender ==sender & shortened_name$receiver ==reciever])==0){
      matrix[i,j]=0
    }
    else{
  matrix[i,j] = shortened_name$interactions[shortened_name$sender ==sender & shortened_name$receiver ==reciever]
    }
  
  }
}
#dont care about a cluster interacting with itself
for(i in 0:30){
  matrix[i,i] =0
}

#summing across the diagonal because I dont care about the directionality
symmetric_matrix <- matrix + t(matrix) - diag(diag(matrix))


###now make strength matrix
dist_matrix <- 1 / (symmetric_matrix)  # Add small constant to avoid division by zero
#microglia and 31 are confusing everything, make them 0
dist_matrix <- ifelse(dist_matrix ==Inf, 0, dist_matrix)
dist_matrix <- dist_matrix[-30,-30] #remove 31
dist_matrix <- dist_matrix[-15,-15] #remove 

#convert to dist object
dist_obj <- as.dist(dist_matrix)

# Apply MDS to get 3D embedding
positions <- cmdscale(dist_obj, k = 3, eig = TRUE)

# Create a data frame with the positions
positions_df <- as.data.frame(positions$points)
colnames(positions_df) <- c("x", "y", "z")
positions_df$cluster <- rownames(positions_df)

library(plotly)
plot <- plot_ly(positions_df, x = ~x, y = ~y, z = ~z, 
        color = ~cluster, 
        type = "scatter3d", 
        mode = "markers+text",
        text = ~cluster,
        marker = list(size = 10))
plot

saveRDS(plot, 'Cell Type Annotations/Spatial Map Via Cell Adhesion Molecules/3D Plot of Spatial Relationship Betwen Clusters')
        

library(igraph)
g <- graph_from_adjacency_matrix(symmetric_matrix[-c(15,30),-c(15,30)], 
                                 mode = "undirected", 
                                 weighted = TRUE)


# Visualize the network
set.seed(123)
plot(g, 
     vertex.size = 15,
     vertex.label = V(g)$name,
     edge.width = E(g)$weight/max(E(g)$weight),  
     layout = layout_with_fr(g))        

