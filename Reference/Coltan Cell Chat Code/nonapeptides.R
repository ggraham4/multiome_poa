# cgp
# rhodes lab uiuc


#== 0. LOAD LIBRARIES ==========================================================
library(tidyverse)
library(readxl)
library(Seurat)


#== 1. LOAD DATA ===============================================================
# load data and subset neurons
setwd("S:/Fish Lab/Experiments/snRNAseq POA sex differences/GT DATA/")
GT_data <- readRDS("GT_data_wClasses_2022-02-07.rds")
GT_neurons <- subset(GT_data, class_two == "Neurons")

# import genes of interest table
genes_interest <- read_excel("genes_interest.xlsx")
genes <- genes_interest$geneID[genes_interest$anycells == TRUE]

genename_full <- genes_interest$gene_name[genes_interest$anycells == TRUE]
genetype <- genes_interest$genetype[genes_interest$anycells == TRUE]
genesystem <- genes_interest$system[genes_interest$anycells == TRUE]


#==== SUBSET JUST GENE POSITIVE CELLS ====

setwd("S:/Fish Lab/Experiments/snRNAseq POA sex differences/GT DATA/results/nonapeptides/")

target <- c("oxt", "avp")

# the loop
for (i in target) {
  print(i)
  subset <- GT_neurons[, GetAssayData(GT_neurons[["RNA"]])[i, ] > 0]
  
  # number and proportion per cluster for each gene
  c20 <- t(as.data.frame(rbind(table(subset$clusters20), prop.table(table(subset$clusters20)))))
  c49 <- t(as.data.frame(rbind(table(subset$clusters49), prop.table(table(subset$clusters49)))))
  
  write.csv(c20, paste("./distribution_coexpression_2022-02-07/", i, "_c20counts.csv", sep = ""))
  write.csv(c49, paste("./distribution_coexpression_2022-02-07/", i, "_c49counts.csv", sep = ""))
  
  # loop to generate the coexpression values
  genename <- c()
  count <- c()
  percent <- c()
  
  for (g in genes) {
    #print(g)
    
    c <- sum(GetAssayData(object = subset, slot = "data")[g,] > 0)
    p <- sum(GetAssayData(object = subset, slot = "data")[g,] > 0) / nrow(subset@meta.data)
    
    genename <- append(genename, g)
    count <- append(count, c)
    percent <- append(percent, p)
  }
  
  df <- data.frame(genename, genename_full, genesystem, genetype, count, percent)
  write.csv(df, paste("./distribution_coexpression_2022-02-07/", i, "_coexpression.csv", sep = ""))
}


#==== 2. SEX DIFF DEGs ====

setwd("S:/Fish Lab/Experiments/snRNAseq POA sex differences/GT DATA/results/nonapeptides/")

for (i in target) {
  print(i)
  subset <- GT_neurons[, GetAssayData(GT_neurons[["RNA"]])[i, ] > 0]
  
  markers <- FindMarkers(subset, 
                         group.by = "sex", 
                         ident.1 = "f", 
                         ident.2 = "m", 
                         slot = "data", 
                         logfc.threshold = 0, 
                         min.pct = 0, 
                         test.use = "wilcox")
  
  markers$p_val_fdr <- p.adjust(markers$p_val, method = "fdr")
  markers$gene <- rownames(markers)
  
  write.csv(markers, 
            file = paste("./DEGs_FvM_2022-02-07/", i, "_DEGs_FvM.csv", sep = ""))
}


#==== SEX DIFF PROPORTIONS ====
#---- FISHER'S EXACT TEST FOR EACH ----

# just need to compare hit (oxt+) miss (oxt-) by sex


# create empty vectors for values to place into
pop_gene <- c()

pRaw_DEGs_sig_ERE <- c()
pRaw_DEGs_sig_noERE <- c()
pRaw_DEGs_ns_ERE <- c()
pRaw_DEGs_ns_noERE <- c()

pRaw_fisher_pval <- c()
pRaw_fisher_OR <- c()
pRaw_fisher_OR_lower <- c()
pRaw_fisher_OR_upper <- c()

pFDR_DEGs_sig_ERE <- c()
pFDR_DEGs_sig_noERE <- c()
pFDR_DEGs_ns_ERE <- c()
pFDR_DEGs_ns_noERE <- c()

pFDR_fisher_pval <- c()
pFDR_fisher_OR <- c()
pFDR_fisher_OR_lower <- c()
pFDR_fisher_OR_upper <- c()

pBon_DEGs_sig_ERE <- c()
pBon_DEGs_sig_noERE <- c()
pBon_DEGs_ns_ERE <- c()
pBon_DEGs_ns_noERE <- c()

pBon_fisher_pval <- c()
pBon_fisher_OR <- c()
pBon_fisher_OR_lower <- c()
pBon_fisher_OR_upper <- c()


# loop over cluster names to create directory for each cluster
for (i in genes) {
  # print gene being worked on
  print(i)
  
  # add gene name to vector
  pop_gene <- append(pop_gene, i)
  
  # pull up DEG by sex output for each gene
  DEGs <- read.csv(paste("./unbiased_gene_defined_pop/", i, "/DEG_by_sex.csv", sep = ""))
  
  # for pRaw, calc fisher's exact test to see if ERE genes overrep
  pRaw_genes_sig <- DEGs$gene[DEGs$p_val < 0.05]
  pRaw_genes_ns <- DEGs$gene[DEGs$p_val > 0.05]
  
  # total number of significant vs. non (column totals)
  pRaw_count_sig <- as.numeric(length(pRaw_genes_sig))
  pRaw_count_ns <- as.numeric(length(pRaw_genes_ns))
  
  # calculate number of ERE-containing genes in each column
  pRaw_count_sig_ERE <- as.numeric(sum(ere_genes %in% pRaw_genes_sig))
  pRaw_count_sig_noERE <- pRaw_count_sig - pRaw_count_sig_ERE
  
  pRaw_count_ns_ERE <- as.numeric(sum(ere_genes %in% pRaw_genes_ns))
  pRaw_count_ns_noERE <- pRaw_count_ns - pRaw_count_ns_ERE
  
  # these are the 4 values for the exact test
  df_pRaw <- data.frame(
    pRaw_sig <- c(pRaw_count_sig_ERE, pRaw_count_sig_noERE),
    pRaw_ns <- c(pRaw_count_ns_ERE, pRaw_count_ns_noERE)
  )
  rownames(df_pRaw) <- c("ERE_present", "ERE_absent")
  colnames(df_pRaw) <- c("pRaw_sig", "pRaw_ns")
  
  png(paste("./ERE_in_DEGs_genedefinedtop100/", i, "_pRaw_mosaic.png", sep = ""))
  mosaicplot(df_pRaw,
             main = "mosaic pRaw",
             color = TRUE)
  dev.off()
  
  fisher_pRaw <- fisher.test(df_pRaw)
  
  # add values to output vectors
  pRaw_DEGs_sig_ERE <- append(pRaw_DEGs_sig_ERE, pRaw_count_sig_ERE)
  pRaw_DEGs_sig_noERE <- append(pRaw_DEGs_sig_noERE, pRaw_count_sig_noERE)
  pRaw_DEGs_ns_ERE <- append(pRaw_DEGs_ns_ERE, pRaw_count_ns_ERE)
  pRaw_DEGs_ns_noERE <- append(pRaw_DEGs_ns_noERE, pRaw_count_ns_noERE)
  
  pRaw_fisher_pval <- append(pRaw_fisher_pval, fisher_pRaw$p.value)
  pRaw_fisher_OR <- append(pRaw_fisher_OR, as.numeric(fisher_pRaw$estimate))
  pRaw_fisher_OR_lower <- append(pRaw_fisher_OR_lower, fisher_pRaw$conf.int[1])
  pRaw_fisher_OR_upper <- append(pRaw_fisher_OR_upper, fisher_pRaw$conf.int[2])
  
  
  # for pFDR, calc fisher's exact test to see if ERE genes overrep
  pFDR_genes_sig <- DEGs$gene[DEGs$p_val_fdr < 0.05]
  pFDR_genes_ns <- DEGs$gene[DEGs$p_val_fdr > 0.05]
  
  # total number of significant vs. non (column totals)
  pFDR_count_sig <- as.numeric(length(pFDR_genes_sig))
  pFDR_count_ns <- as.numeric(length(pFDR_genes_ns))
  
  # calculate number of ERE-containing genes in each column
  pFDR_count_sig_ERE <- as.numeric(sum(ere_genes %in% pFDR_genes_sig))
  pFDR_count_sig_noERE <- pFDR_count_sig - pFDR_count_sig_ERE
  
  pFDR_count_ns_ERE <- as.numeric(sum(ere_genes %in% pFDR_genes_ns))
  pFDR_count_ns_noERE <- pFDR_count_ns - pFDR_count_ns_ERE
  
  # these are the 4 values for the exact test
  df_pFDR <- data.frame(
    pFDR_sig <- c(pFDR_count_sig_ERE, pFDR_count_sig_noERE),
    pFDR_ns <- c(pFDR_count_ns_ERE, pFDR_count_ns_noERE)
  )
  rownames(df_pFDR) <- c("ERE_present", "ERE_absent")
  colnames(df_pFDR) <- c("pFDR_sig", "pFDR_ns")
  
  png(paste("./ERE_in_DEGs_genedefinedtop100/", i, "_pFDR_mosaic.png", sep = ""))
  mosaic_pFDR <- mosaicplot(df_pFDR,
                            main = "mosaic pFDR",
                            color = TRUE)
  dev.off()
  
  fisher_pFDR <- fisher.test(df_pFDR)
  
  # add to output vectors
  pFDR_DEGs_sig_ERE <- append(pFDR_DEGs_sig_ERE, pFDR_count_sig_ERE)
  pFDR_DEGs_sig_noERE <- append(pFDR_DEGs_sig_noERE, pFDR_count_sig_noERE)
  pFDR_DEGs_ns_ERE <- append(pFDR_DEGs_ns_ERE, pFDR_count_ns_ERE)
  pFDR_DEGs_ns_noERE <- append(pFDR_DEGs_ns_noERE, pFDR_count_ns_noERE)
  
  pFDR_fisher_pval <- append(pFDR_fisher_pval, fisher_pFDR$p.value)
  pFDR_fisher_OR <- append(pFDR_fisher_OR, as.numeric(fisher_pFDR$estimate))
  pFDR_fisher_OR_lower <- append(pFDR_fisher_OR_lower, fisher_pFDR$conf.int[1])
  pFDR_fisher_OR_upper <- append(pFDR_fisher_OR_upper, fisher_pFDR$conf.int[2])
  
  
  # for pBon, calc fisher's exact test to see if ERE genes overrep
  pBon_genes_sig <- DEGs$gene[DEGs$p_val_adj < 0.05]
  pBon_genes_ns <- DEGs$gene[DEGs$p_val_adj > 0.05]
  
  # total number of significant vs. non (column totals)
  pBon_count_sig <- as.numeric(length(pBon_genes_sig))
  pBon_count_ns <- as.numeric(length(pBon_genes_ns))
  
  # calculate number of ERE-containing genes in each column
  pBon_count_sig_ERE <- as.numeric(sum(ere_genes %in% pBon_genes_sig))
  pBon_count_sig_noERE <- pBon_count_sig - pBon_count_sig_ERE
  
  pBon_count_ns_ERE <- as.numeric(sum(ere_genes %in% pBon_genes_ns))
  pBon_count_ns_noERE <- pBon_count_ns - pBon_count_ns_ERE
  
  # these are the 4 values for the exact test
  df_pBon <- data.frame(
    pBon_sig <- c(pBon_count_sig_ERE, pBon_count_sig_noERE),
    pBon_ns <- c(pBon_count_ns_ERE, pBon_count_ns_noERE)
  )
  rownames(df_pBon) <- c("ERE_present", "ERE_absent")
  colnames(df_pBon) <- c("pBon_sig", "pBon_ns")
  
  png(paste("./ERE_in_DEGs_genedefinedtop100/", i, "_pBon_mosaic.png", sep = ""))
  mosaicplot(df_pBon,
             main = "mosaic pBon",
             color = TRUE)
  dev.off()
  
  fisher_pBon <- fisher.test(df_pBon)
  
  # add results to output vectors
  pBon_DEGs_sig_ERE <- append(pBon_DEGs_sig_ERE, pBon_count_sig_ERE)
  pBon_DEGs_sig_noERE <- append(pBon_DEGs_sig_noERE, pBon_count_sig_noERE)
  pBon_DEGs_ns_ERE <- append(pBon_DEGs_ns_ERE, pBon_count_ns_ERE)
  pBon_DEGs_ns_noERE <- append(pBon_DEGs_ns_noERE, pBon_count_ns_noERE)
  
  pBon_fisher_pval <- append(pBon_fisher_pval, fisher_pBon$p.value)
  pBon_fisher_OR <- append(pBon_fisher_OR, as.numeric(fisher_pBon$estimate))
  pBon_fisher_OR_lower <- append(pBon_fisher_OR_lower, fisher_pBon$conf.int[1])
  pBon_fisher_OR_upper <- append(pBon_fisher_OR_upper, fisher_pBon$conf.int[2])
  
}

# create dataframe
df <- data.frame(pop_gene, 
                 
                 pRaw_DEGs_sig_ERE,
                 pRaw_DEGs_sig_noERE,
                 pRaw_DEGs_ns_ERE,
                 pRaw_DEGs_ns_noERE,
                 
                 pRaw_fisher_pval,
                 pRaw_fisher_OR,
                 pRaw_fisher_OR_lower,
                 pRaw_fisher_OR_upper,
                 
                 pFDR_DEGs_sig_ERE,
                 pFDR_DEGs_sig_noERE,
                 pFDR_DEGs_ns_ERE,
                 pFDR_DEGs_ns_noERE,
                 
                 pFDR_fisher_pval,
                 pFDR_fisher_OR,
                 pFDR_fisher_OR_lower,
                 pFDR_fisher_OR_upper,
                 
                 pBon_DEGs_sig_ERE,
                 pBon_DEGs_sig_noERE,
                 pBon_DEGs_ns_ERE,
                 pBon_DEGs_ns_noERE,
                 
                 pBon_fisher_pval,
                 pBon_fisher_OR,
                 pBon_fisher_OR_lower,
                 pBon_fisher_OR_upper
)

df$pRaw_fisher_pval_FDR <- p.adjust(df$pRaw_fisher_pval, method = "fdr")
df$pRaw_fisher_pval_Bon <- p.adjust(df$pRaw_fisher_pval, method = "bonferroni")

df$pFDR_fisher_pval_FDR <- p.adjust(df$pFDR_fisher_pval, method = "fdr")
df$pFDR_fisher_pval_Bon <- p.adjust(df$pFDR_fisher_pval, method = "bonferroni")

df$pBon_fisher_pval_FDR <- p.adjust(df$pBon_fisher_pval, method = "fdr")
df$pBon_fisher_pval_Bon <- p.adjust(df$pBon_fisher_pval, method = "bonferroni")


# write csv
setwd("S:/Fish Lab/Experiments/snRNAseq POA sex differences/GT DATA/results/ERE_in_DEGs_genedefinedtop100/")

write.csv(df, file = "ERE_in_DEG_fisher_results.csv")












# plotting function to convince yourself that it's true
DimPlot(subset, 
        reduction = "umap", 
        label = FALSE)

FeaturePlot(GT_data,
            features = i, 
            cols = c("gray90", "blue"),
            order = TRUE)







#---- LIST OF TOP 100 GENES ----
# make a list of genes to loop over
genes <- genes_interest$gene[1:100]

# remove troublesome genes
genes <- genes[-c(13)]
genes <- genes[-c(51)]


#---- CREATE DIRECTORIES TO STORE OUTPUT ----
# first need to set working directory to create subdirectory for each cluster
setwd("S:/Fish Lab/Experiments/snRNAseq POA sex differences/GT DATA/results/unbiased_gene_defined_pop/")

# loop over cluster names to create directory for each cluster
for (i in genes) {
  dir.create(i)
}

# create a vector to store cell # in for inspection
cellnumber <- c()

# DEG analysis for marker+ cells compared against all other cells
# resource consulted for loop function: https://github.com/satijalab/seurat/issues/2619
for (i in genes) {
  print(i)
  subset <- GT_data[, GetAssayData(GT_data[["RNA"]])[i, ] > 0]
  
  print(table(subset$sex))
  table <- (table(subset$sex))
  cells <- as.numeric(sum(table))
  cellnumber <- append(cellnumber, cells)
  
  markers <- FindMarkers(subset, 
                         group.by = "sex", 
                         ident.1 = "f", 
                         ident.2 = "m", 
                         slot = "data", 
                         logfc.threshold = 0, 
                         min.pct = 0, 
                         test.use = "wilcox")
  markers$p_val_fdr <- p.adjust(markers$p_val, method = "fdr")
  markers$gene <- rownames(markers)
  
  write.csv(markers, file = paste(i, "/DEG_by_sex.csv", sep = ""))
}


#---- GENERATE SUMMARY SHEET ----

# set working directory and read in summary file
setwd("S:/Fish Lab/Experiments/snRNAseq POA sex differences/GT data/results/")

# create empty vectors for values to place into
pop_gene <- c()
sig_pRaw <- c()
sig_pBon <- c()
sig_pFDR <- c()

# loop over cluster names to create directory for each cluster
for (i in genes) {
  # print gene being worked on
  print(i)
  
  # add gene name to vector
  pop_gene <- append(pop_gene, i)
  
  # pull up DEG by sex output for each gene
  DEGs <- read.csv(paste("./unbiased_gene_defined_pop/", i, "/DEG_by_sex.csv", sep = ""))
  
  # grab number of significant p-values for each gene
  pRaw <- as.numeric(DEGs$p_val < 0.05)
  pRaw <- sum(pRaw, na.rm = TRUE)
  sig_pRaw <- append(sig_pRaw, pRaw)
  
  pBon <- as.numeric(DEGs$p_val_adj < 0.05)
  pBon <- sum(pBon, na.rm = TRUE)
  sig_pBon <- append(sig_pBon, pBon)
  
  pFDR <- as.numeric(DEGs$p_val_fdr < 0.05)
  pFDR <- sum(pFDR, na.rm = TRUE)
  sig_pFDR <- append(sig_pFDR, pFDR)
}

# create dataframe
df <- data.frame(pop_gene, sig_pRaw, sig_pBon, sig_pFDR)


# create empty vectors for cell counts
female <- c()
male <- c()
cells <- c()

# add cell number to each
for (i in genes) {
  print(i)
  subset <- GT_data[, GetAssayData(GT_data[["RNA"]])[i, ] > 0]
  
  print(table(subset$sex))
  table <- (table(subset$sex))
  f <- as.numeric(table[1])
  m <- as.numeric(table[2])
  c <- sum(f, m)
  
  female <- append(female, f)
  male <- append(male, m)
  cells <- append(cells, c)
}

# create dataframe
df <- data.frame(pop_gene, sig_pRaw, sig_pBon, sig_pFDR, cells, female, male)

df$female_cells_pct <- df$female/df$cells * 100
df$male_cells_pct <- df$male/df$cells * 100

df$pct_diff <- df$female_cells_pct - df$male_cells_pct

# write csv
setwd("S:/Fish Lab/Experiments/snRNAseq POA sex differences/GT DATA/results/unbiased_gene_defined_pop/")

write.csv(df, file = "DEG_by_sex_unbiased_pop_summary.csv")


# ---- STATISTICS ----
df$rank <- as.numeric(rownames(df))

cor.test(df$sig_pFDR, df$rank)
cor.test(df$sig_pFDR, df$cells)


summary(lm(df$sig_pFDR ~ df$rank * df$cells * df$pct_diff))
Anova(lm(df$sig_pFDR ~ df$rank * df$cells * df$pct_diff), type = "III")


Anova(lm(df$sig_pFDR ~ df$rank + df$cells + df$pct_diff), type = "III")





# double checking that the genes to be tested have enough + cells
for (i in genes) {
  print(i)
  genename <- genes_interest$gene_name[genes_interest$geneID == i]
  subset <- GT_data[, GetAssayData(GT_data[["RNA"]])[i, ] > 0]
  print(table(subset$sex))
}

# DEG analysis for marker+ cells compared against all other cells
# resource consulted for loop function: https://github.com/satijalab/seurat/issues/2619
for (i in genes) {
  print(i)
  genename <- genes_interest$gene_name[genes_interest$geneID == i] #
  subset <- GT_data[, GetAssayData(GT_data[["RNA"]])[i, ] > 0]
  
  print(table(subset$sex))
  table <- (table(subset$sex))
  if ((table[1] > 3) & (table[2] > 3)) {
    markers <- FindMarkers(subset, 
                           group.by = "sex", 
                           ident.1 = "m", 
                           ident.2 = "f", 
                           slot = "data", 
                           logfc.threshold = 0, 
                           min.pct = 0, 
                           test.use = "wilcox")
    markers$p_val_fdr <- p.adjust(markers$p_val, method = "fdr")
    write.csv(markers, file = paste("./results/DGE_marker_by_sex/", genename, ".csv", sep = ""))
  }
}






# BELOW IS A TEST CODE FOR A LARGER LOOP

# testing out with embedded "if" statement
# DEG analysis for marker+ cells compared against all other cells
# resource consulted: https://github.com/satijalab/seurat/issues/2619

# first create the folders
setwd("S:/Fish Lab/Experiments/snRNAseq POA sex differences/GT DATA/results/by_gene_marker/")
for (i in genes) {
  dir.create(i)
}

# then fill the folders with data
for (i in genes) {
  # print the gene being worked on
  print(i)
  genename <- genes_interest$gene_name[genes_interest$geneID == i]
  
  # setwd
  setwd(paste0("S:/Fish Lab/Experiments/snRNAseq POA sex differences/GT DATA/results/by_gene_marker/", i))
  
  # subset just cells expressing the gene
  subset <- GT_data[, GetAssayData(GT_data[["RNA"]])[i, ] > 0]
  
  # if there is enough gene+ cells do DEG analysis, else pass
  print(table(subset$sex))
  table <- (table(subset$sex))
  if ((table[1] > 3) & (table[2] > 3)) {
    markers <- FindMarkers(subset, 
                           group.by = "sex", 
                           ident.1 = "m", 
                           ident.2 = "f", 
                           slot = "data", 
                           logfc.threshold = 0, 
                           min.pct = 0, 
                           test.use = "wilcox")
    markers$p_val_fdr <- p.adjust(markers$p_val, method = "fdr")
    write.csv(markers, file = paste(i, "_DEG-by-sex.csv", sep = ""))
  } else {
    print("not enough cells for DEG")
  }
  
  # making plots
  print("plotting")
  temp_plot <- FeaturePlot(GT_data, features = i, order = TRUE)
  ggsave(filename = paste(i, "_umap.pdf", sep = ""), 
         plot = temp_plot, device = "pdf", 
         units = "cm", 
         width = 25, 
         height = 20)
}


ar <- read.csv("./results/by_gene_marker/ar/ar_DEG-by-sex.csv")

keyvals <- ifelse((ar$avg_log2FC < 0) & (ar$p_val_fdr < 0.05), 'purple',
                  ifelse((ar$avg_log2FC > 0) & (ar$p_val_fdr < 0.05), 'yellow', 'gray'))

keyvals[is.na(keyvals)] <- 'gray'
names(keyvals)[keyvals == 'yellow'] <- 'high'
names(keyvals)[keyvals == 'gray'] <- 'n.s.'
names(keyvals)[keyvals == 'purple'] <- 'low'

siggenes <- na.omit(ar$X[ar$p_val_fdr < 0.05])

volcano <- EnhancedVolcano(ar, 
                lab = ar$X, 
                x = 'avg_log2FC', 
                y = 'p_val_fdr',
                ylim = c(0, -log10(min(ar$p_val_fdr, na.rm = TRUE))),
                pCutoff = 0.05,
                FCcutoff = 1,
                colCustom = keyvals,
                drawConnectors = TRUE,
                max.overlaps = Inf)
ggsave("volcano.png", plot = volcano, device = "png", units = "cm", width = 20, height = 20)



#---- working
input <- read.csv("./results/by_gene_marker/ar/DEG-by-sex_wilcox.csv")




