# COLTAN G PARKER
# RHODES LAB UIUC

# A ocellaris POA snRNAseq
# Neurotransmitter Interaction Network
# 1. calculate the percent of cells of a given NT identity expressing any receptor of a given NT
# 2. put data into a matrix and generate a heatmap 

# example row: GnRH x GnRH (%), x DA, x GnIH, x OXT, etc.


#== 0. LOAD LIBRARIES ==========================================================

library(tidyverse)
library(readxl)
#library(igraph)
library(tidygraph)
library(ggraph)

library(Seurat)


#== 1. LOAD DATA ===============================================================

# load data and subset neurons
setwd("S:/Fish Lab/Experiments/snRNAseq POA sex differences/GT DATA/")
GT_data <- readRDS("ocellaris2020_labeled_v02.rds")
GT_neurons <- subset(GT_data, class_two == "Neurons")
rm(GT_data)





# import genes of interest table
genes_interest <- read_excel("genes_interest.xlsx")
genes_interest <- filter(genes_interest, genes_interest$anycells == TRUE)

genes_network <- filter(genes_interest, genes_interest$genetype == "receptor")

genes_network <- append(unique(genes_network$system[genes_network$category == "neurotransmitter"]), 
                        unique(genes_network$system[genes_network$category == "neuropeptide"]))

# these are the genes that have receptors
# and can be included in the network
# in the future we should add Glu and GABA

# create objects for nodes
DA <- subset(GT_neurons, subset = th > 0 | th2 > 0)
Ach <- subset(GT_neurons, subset = CHAT > 0 | slc18a3b > 0)

GnRH <- subset(GT_neurons, subset = ENSAPEG00000007875 > 0)
GnIH <- subset(GT_neurons, subset = ENSAPEG00000019310 > 0)
Kiss <- subset(GT_neurons, subset = ENSAPEG00000002502 > 0)
Oxt <- subset(GT_neurons, subset = oxt > 0)
Avt <- subset(GT_neurons, subset = avp > 0)
Gal <- subset(GT_neurons, subset = galn > 0)
NPY <- subset(GT_neurons, subset = NPY > 0 | npy > 0)
Tac <- subset(GT_neurons, subset = tac1 > 0 | tac3a > 0)
CCK <- subset(GT_neurons, subset = ccka > 0 | cckb > 0)
PACAP <- subset(GT_neurons, subset = ADCYAP1 > 0 | adcyap1b > 0)
VIP <- subset(GT_neurons, subset = vip > 0)
GHRH <- subset(GT_neurons, subset = ghrh > 0)

nodes <- c(DA, Ach, GnRH, GnIH, Kiss, Oxt, Avt, Gal, NPY, Tac, CCK, PACAP, VIP, GHRH)


# list of receptors to look for
DA_R <- genes_interest$geneID[genes_interest$system == "dopamine" & genes_interest$genetype == "receptor"]
Ach_R <- genes_interest$geneID[genes_interest$system == "acetylcholine" & genes_interest$genetype == "receptor"]

GnRH_R <- genes_interest$geneID[genes_interest$system == "GnRH" & genes_interest$genetype == "receptor"]
GnIH_R <- genes_interest$geneID[genes_interest$system == "GnIH" & genes_interest$genetype == "receptor"]
Kiss_R <- genes_interest$geneID[genes_interest$system == "kisspeptin" & genes_interest$genetype == "receptor"]
Oxt_R <- genes_interest$geneID[genes_interest$system == "oxytocin" & genes_interest$genetype == "receptor"]
Avt_R <- genes_interest$geneID[genes_interest$system == "vasotocin" & genes_interest$genetype == "receptor"]
Gal_R <- genes_interest$geneID[genes_interest$system == "galanin" & genes_interest$genetype == "receptor"]
NPY_R <- genes_interest$geneID[genes_interest$system == "NPY" & genes_interest$genetype == "receptor"]
Tac_R <- genes_interest$geneID[genes_interest$system == "tachykinin" & genes_interest$genetype == "receptor"]
CCK_R <- genes_interest$geneID[genes_interest$system == "CCK" & genes_interest$genetype == "receptor"]
PACAP_R <- genes_interest$geneID[genes_interest$system == "PACAP" & genes_interest$genetype == "receptor"]
VIP_R <- genes_interest$geneID[genes_interest$system == "VIP" & genes_interest$genetype == "receptor"]
GHRH_R <- genes_interest$geneID[genes_interest$system == "GHRH" & genes_interest$genetype == "receptor"]

edges <- list(DA_R, Ach_R, GnRH_R, GnIH_R, Kiss_R, Oxt_R, Avt_R, Gal_R, NPY_R, Tac_R, CCK_R, PACAP_R, VIP_R, GHRH_R)


# general sketch of logic
# loop through node objects: for each node...
# pull the total number of cells of cells, and then the number of cells expressing any of one receptor
# maybe the solution is to subset again based on expression of each receptor set, then get cell number

# we have nodes and edges
# want an output matrix of 14x14
# where the value of each cell is the number of cells (%) in the object (row)
# expressing at least one of the receptors from each system

outputdf <- data.frame(matrix(ncol = 14, nrow = 0)) # need to give 14 columns

for (i in 1:14) {
  n_cells <- length(colnames(nodes[[i]])) # get number of cells in node
  test <- as.data.frame(t(as.matrix(GetAssayData(nodes[[i]], slot = "data")))) # pull out expression matrix to get cell count
  
  row <- c()
  
  for (j in 1:14) {
    edge_focal <- edges[[j]]
    
    test_working <- test
    
    for (k in edge_focal) {
      #print(k)
      gene_index <- which(colnames(test_working) == k)
      test_working <- test_working[test_working[,k] == 0, ] # keep only cells with 0 expression
      # we are left with test containing only the cells that are 100% negative
      # so n_cells - length(rownames(test)) will give the number of cells that are positive for at least one R
      # then that value / n_cells is the proportion of total cells expressing at least one R
    }
    
    n_negative <- length(rownames(test_working))
    n_positive <- n_cells - n_negative
    pct_positive <- (n_positive / n_cells) * 100
    
    row <- append(row, pct_positive)
  }
  
  write.csv(row, paste("./temp/row", i, ".csv"))
}


# get number of cells in node expressing at least one of receptor in set
# calculate A / B *100 = percent cells in node expressing at least one receptor


# read in dummy set 
comp <- read_excel("./temp/compilation.xlsx")
comp <- subset(comp, select = -c(node))
rownames(comp) <- colnames(comp)
heatmap(as.matrix(comp))

comp <- read_excel("./temp/compilation.xlsx")
comp <- pivot_longer(comp, cols = 2:13, names_to = "receptor", values_to = "percent")

# reformatting the data properly
# from https://www.jessesadler.com/post/network-analysis-with-r/

nodes <- comp %>%
  distinct(node) %>%
  rename(label = node) %>%
  rowid_to_column("id")

edges <- comp %>%
  left_join(nodes, by = c("node" = "label")) %>%
  rename(to = id) %>%
  left_join(nodes, by = c("receptor" = "label")) %>%
  rename(from = id)

edges <- select(edges, to, from, percent)

# okay let's try igraph
igraphtest <- graph_from_data_frame(d = edges, vertices = nodes, directed = TRUE)

plot(igraphtest,
     edge.arrow.size = 0.4)

# try tidygraph and ggraph
ggraphtest <- tbl_graph(nodes = nodes, 
                        edges = edges, 
                        directed = TRUE)
ggraph(ggraphtest) + 
  geom_node_point() + 
  geom_edge_link(aes(width = percent, alpha = percent)) +
  scale_edge_width(range = c(0.1, 4)) + 
  geom_node_text(aes(label = label), repel = TRUE) + 
  theme_graph()

# decide which edges may be worth removing
hist(comp$percent, breaks = 50)
hist(comp$percent, breaks = 100)
# remove all edges that are < 10 %

edges90 <- subset(edges, percent > 10)
edges80 <- subset(edges, percent > 20)

ggraph90 <- tbl_graph(nodes = nodes, 
                      edges = edges90, 
                      directed = TRUE)
ggraph80 <- tbl_graph(nodes = nodes, 
                      edges = edges80, 
                      directed = TRUE)

ggraph(ggraph80) + 
  geom_node_point(size = 4) + 
  #geom_edge_link(aes(width = percent, alpha = percent)) +
  geom_edge_fan(arrow = arrow(),
                end_cap = circle(0.2),
                aes(width = percent, 
                    alpha = percent)) +
  scale_edge_width(range = c(0.1, 1)) + 
  geom_node_text(aes(label = label), repel = TRUE) + 
  theme_graph()

plot <- ggraph(ggraph80, layout = "linear") + 
  geom_node_point(size = 16, shape = 1) + 
  #geom_edge_link(aes(width = percent, alpha = percent)) +
  geom_edge_arc(arrow = arrow(length = unit(4, "mm")),
                start_cap = circle(0.8),
                end_cap = circle(0.8),
                aes(width = percent, 
                    alpha = percent)) +
  scale_edge_width(range = c(0.1, 1)) + 
  geom_node_text(aes(label = label), repel = FALSE) + 
  theme_graph()
ggsave(plot, 
       filename = "network_v01.png", 
       device = "png", 
       units = "px", 
       width = 2400, 
       height = 2000, 
       path = "./results/")


#==== PLAN FOR V02 ====

#1 make a seurat object with columns in metadata for each cell type
#2 make two "lookup tibbles", one with id x node and one with id x receptor 
# 2.1. the second lookup tibble will have receptors beyond nodes (i.e. steroids,etc)
# so the heatmap will be asymetric but that is okay
#3 go through and do the same script but subsetting the main seurat object so there is not all that waste

# is there a way to devise two weight scores?
# maybe even just use two, wherein pct is used to weight arrow thickness
# and a second wherein avg expression level within expressing cells is color

# identify the clusters

# pull out metadata
meta <- GT_neurons@meta.data
cells <- rownames(meta)

# function for comparing the list of cells expressing a given gene
# with the list of all cells
makeTFlist <- function(x, y) {
  y %in% x
}

# now getting a vector of T/F for whether the cells are positive
ach <- WhichCells(GT_neurons, expression = CHAT > 0 | slc18a3b > 0) %>%
  makeTFlist(cells)
da <- WhichCells(GT_neurons, expression = th > 0 | th2 > 0 | slc6a3) %>%
  makeTFlist(cells)
ne <- WhichCells(GT_neurons, expression = dbh > 0 | slc6a2 > 0) %>%
  makeTFlist(cells)
sero <- WhichCells(GT_neurons, expression = TPH1 > 0 | tph2 > 0 | slc6a4a > 0) %>%
  makeTFlist(cells)
hist <- WhichCells(GT_neurons, expression = hdc > 0) %>%
  makeTFlist(cells)

gnrh <- WhichCells(GT_neurons, expression = ENSAPEG00000007875 > 0) %>%
  makeTFlist(cells)
gnih <- WhichCells(GT_neurons, expression = ENSAPEG00000019310 > 0) %>%
  makeTFlist(cells)
kiss <- WhichCells(GT_neurons, expression = ENSAPEG00000002502 > 0) %>%
  makeTFlist(cells)
oxt <- WhichCells(GT_neurons, expression = oxt > 0) %>%
  makeTFlist(cells)
avt <- WhichCells(GT_neurons, expression = avp > 0) %>%
  makeTFlist(cells)
gal <- WhichCells(GT_neurons, expression = galn > 0) %>%
  makeTFlist(cells)
npy <- WhichCells(GT_neurons, expression = NPY > 0 | npy > 0) %>%
  makeTFlist(cells)
ghrh <- WhichCells(GT_neurons, expression = ghrh > 0) %>%
  makeTFlist(cells)
cck <- WhichCells(GT_neurons, expression = ccka > 0 | cckb > 0) %>%
  makeTFlist(cells)





# add columns to metadata
meta <- mutate(meta,
               ach = ach)
GT_neurons@meta.data <- meta

# clusters <- list of names of clusters in metadata
test <- subset(GT_neurons, ach == TRUE)


# TRYING A DIFFERENT STRATEGY
id_key <- data.frame(1:15,
                     c("ach", "da", "ne"))

cell_lists <- list(WhichCells(GT_neurons, expression = CHAT > 0 | slc18a3b > 0), 
                   WhichCells(GT_neurons, expression = th > 0 | th2 > 0 | slc6a3))






for (i in clusters){
  #subset <- the seurat object by cluster column (cell identity) = TRUE
  #ncells <- number of cells in the subset
  
  #
  
}


for (i in 1:14) {
  n_cells <- length(colnames(nodes[[i]])) # get number of cells in node
  test <- as.data.frame(t(as.matrix(GetAssayData(nodes[[i]], slot = "data")))) # pull out expression matrix to get cell count
  
  row <- c()
  
  for (j in 1:14) {
    edge_focal <- edges[[j]]
    
    test_working <- test
    
    for (k in edge_focal) {
      #print(k)
      gene_index <- which(colnames(test_working) == k)
      test_working <- test_working[test_working[,k] == 0, ] # keep only cells with 0 expression
      # we are left with test containing only the cells that are 100% negative
      # so n_cells - length(rownames(test)) will give the number of cells that are positive for at least one R
      # then that value / n_cells is the proportion of total cells expressing at least one R
    }
    
    n_negative <- length(rownames(test_working))
    n_positive <- n_cells - n_negative
    pct_positive <- (n_positive / n_cells) * 100
    
    row <- append(row, pct_positive)
  }
  
  write.csv(row, paste("./temp/row", i, ".csv"))
}

  
  
