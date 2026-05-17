library(Seurat)
library(tidyverse)
library(org.Dr.eg.db)

# ── 1. Load data — ENSEMBL_ID column becomes rownames ────────────────────────
dat <- read.csv("~/Downloads/GSE137525_counts.csv", row.names = 1)

# ── 2. Map Ensembl IDs → gene symbols using local annotation db ──────────────
id2sym <- AnnotationDbi::select(
  org.Dr.eg.db,
  keys    = rownames(dat),
  columns = c("ENSEMBL", "SYMBOL"),
  keytype = "ENSEMBL"
) %>%
  filter(!is.na(SYMBOL), SYMBOL != "") %>%
  distinct(ENSEMBL, .keep_all = TRUE)

cat("Genes mapped:", nrow(id2sym), "\n")

# ── 3. Subset dat to mapped genes and rename rows to symbols ─────────────────
dat_sym <- dat[id2sym$ENSEMBL, ]
dat_sym$symbol <- id2sym$SYMBOL

# For duplicated symbols, keep the row with highest total counts
dat_sym <- dat_sym %>%
  mutate(total_counts = rowSums(across(where(is.numeric)))) %>%
  group_by(symbol) %>%
  slice_max(total_counts, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  dplyr::select(-total_counts) %>%
  column_to_rownames("symbol")

cat("Genes after deduplication:", nrow(dat_sym), "\n") # rename to gene symbols

# ── 4. Create Seurat object ───────────────────────────────────────────────────
seu <- CreateSeuratObject(
  counts       = dat_sym,
  project      = "GSE137525_Lange2020",
  min.cells    = 3,
  min.features = 200
)

# ── 5. Add metadata: prefix before first underscore ──────────────────────────

seu =  NormalizeData(seu, normalization.method = "LogNormalize", scale.factor = 10000)
seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(seu)
seu <- ScaleData(seu, features = all.genes)
seu <- RunPCA(seu, features = VariableFeatures(object = seu))
ElbowPlot(seu)
# 8
seu <- FindNeighbors(seu, dims = 1:8)
seu <- FindClusters(seu, resolution = 0.1)
seu <- RunUMAP(seu, dims = 1:8)
DimPlot(seu, reduction = "umap")
seu = FindSubCluster(seu, 3, graph.name= 'RNA_snn')
DimPlot(seu, group.by = 'sub.cluster')

FeaturePlot(seu, 'tubb5')
FeaturePlot(seu, 'aplnra')
FeaturePlot(seu, 'sybu')

seu$annot = ifelse(seu$seurat_clusters==0, 'RGC', NA)
seu$annot = ifelse(seu$seurat_clusters==1, 'NB.1', seu$annot)
seu$annot = ifelse(seu$seurat_clusters==2 , 'NB.2', seu$annot)
seu$annot = ifelse(seu$sub.cluster=='3_0' , 'MN', seu$annot)
seu$annot = ifelse(seu$sub.cluster=='3_1' , 'OPC', seu$annot)


seu$sample_group <- gsub("[^A-Za-z]", "", sub("_.*", "", colnames(seu)))
table(seu$sample_group)

DimPlot(seu, reduction = "umap", group.by = 'annot')

marks_dr = FindAllMarkers(seu, group.by = 'annot')



# ── 6. Load clownfish object ──────────────────────────────────────────────────
obj <- readRDS("~/Desktop/optimal_clustering_rna_only.rds")

#### FUCK IT FUCK ENSEMBL WE ARE ONLY DOING GENES THAT HAVE THE SAME NAMES###
## ── 7. Build zebrafish ↔ clownfish ortholog table via BiomaRt ─────────────────
## Uses AnnotationDbi symbols (SYMBOL column) as input, not Ensembl IDs
## Mirrors used for robustness — try next if one fails
#mart_zf <- tryCatch(
#  useEnsembl("ensembl", dataset = "drerio_gene_ensembl", mirror = "useast"),
#  error = function(e) useEnsembl("ensembl", dataset = "drerio_gene_ensembl", mirror = "uswest")
#)

## Query 1: Ensembl gene ID → ZFIN symbol
#zfin_map <- getBM(
#  attributes = c("ensembl_gene_id", "zfin_id_symbol"),
#  filters    = "ensembl_gene_id",
#  values     = id2sym$ENSEMBL,
#  mart       = mart_zf
#) %>%
#  filter(zfin_id_symbol != "") %>%
#  distinct(ensembl_gene_id, .keep_all = TRUE)

## Query 2: Ensembl gene ID → clownfish ortholog
#cf_map <- getBM(
#  attributes = c(
 #   "ensembl_gene_id",
#    "aocellaris_homolog_associated_gene_name",
#    "aocellaris_homolog_orthology_confidence"
#  ),
#  filters = "ensembl_gene_id",
#  values  = id2sym$ENSEMBL,
#  mart    = mart_zf
#) %>%
#  filter(
#    aocellaris_homolog_associated_gene_name != "",
#    aocellaris_homolog_orthology_confidence == 1
#  ) %>%
#  distinct(ensembl_gene_id, .keep_all = TRUE)

## Join both via Ensembl ID, then bring in the SYMBOL from id2sym
#orthologs_filt <- zfin_map %>%
#  inner_join(cf_map, by = "ensembl_gene_id") %>%
#  inner_join(id2sym %>% dplyr::select(ENSEMBL, SYMBOL),
#             by = c("ensembl_gene_id" = "ENSEMBL")) %>%
#  filter(SYMBOL %in% rownames(seu)) %>%
#  distinct(SYMBOL, .keep_all = TRUE) %>%
#  distinct(aocellaris_homolog_associated_gene_name, .keep_all = TRUE)

#cat("Ortholog pairs retained:", nrow(orthologs_filt), "\n")

## ── 8. Subset zebrafish to orthologs and rename features to clownfish symbols ─
#seu_sub   <- seu[orthologs_filt$SYMBOL, ]

## If attribute name errors, run:
## listAttributes(mart_zf) %>% filter(str_detect(name, "ocellaris"))

#orthologs_filt <- orthologs %>%
#  filter(
  #  zfin_id_symbol != "",
  #  aocellaris_homolog_associated_gene_name != "",
#    aocellaris_homolog_orthology_confidence == 1
#  ) %>%
#  distinct(zfin_id_symbol, .keep_all = TRUE) %>%
#  distinct(aocellaris_homolog_associated_gene_name, .keep_all = TRUE)

#cat("Ortholog pairs retained:", nrow(orthologs_filt), "\n")

## ── 8. Subset zebrafish to orthologs and rename features to clownfish symbols ─
#seu_sub   <- seu[orthologs_filt$zfin_id_symbol, ]

#counts_zf <- GetAssayData(seu_sub, layer = "counts")
#rownames(counts_zf) <- orthologs_filt$aocellaris_homolog_associated_gene_name[
#  match(rownames(counts_zf), orthologs_filt$zfin_id_symbol)
#]

#seu_zf_renamed <- CreateSeuratObject(
#  counts    = counts_zf,
#  meta.data = seu_sub@meta.data
#)

# ── 9. Subset clownfish to shared orthologous genes ───────────────────────────
#cf_keep <- intersect(rownames(seu_zf_renamed), rownames(obj))
#obj_sub         <- obj[cf_keep, ]
#seu_zf_renamed  <- seu_zf_renamed[cf_keep, ]

#cat("Shared features after filtering:", length(cf_keep), "\n")

## ── 10. Normalize and reduce zebrafish reference ──────────────────────────────
#seu_zf_renamed <- NormalizeData(seu_zf_renamed)
#seu_zf_renamed <- FindVariableFeatures(seu_zf_renamed, nfeatures = 3000)
#seu_zf_renamed <- ScaleData(seu_zf_renamed)
#seu_zf_renamed <- RunPCA(seu_zf_renamed)

# ── 11. Get marker genes per cell type from zebrafish reference ───────────────
# marks_dr was already computed above with FindAllMarkers(seu, group.by = 'annot')

marks_filtered <- marks_dr %>%
  filter(p_val_adj < 0.05, avg_log2FC > 0.5) %>%
  group_by(cluster) %>%
  slice_max(avg_log2FC, n = 50) %>%          # top 50 markers per cell type
  ungroup()

cell_types <- unique(marks_filtered$cluster)

# Build named list of gene vectors; keep only genes present in the query object
gene_lists <- lapply(cell_types, function(ct) {
  genes <- marks_filtered %>%
    filter(cluster == ct) %>%
    pull(gene)
  intersect(genes, rownames(obj))
})
names(gene_lists) <- cell_types

# Report overlap
cat("Gene overlap per cell type:\n")
lapply(names(gene_lists), function(ct) cat(ct, ":", length(gene_lists[[ct]]), "genes\n"))

# ── 12. Score with AUCell ─────────────────────────────────────────────────────
library(AUCell)

# Build rankings from the query object's count matrix
counts_mat <- GetAssayData(obj, layer = "counts") %>% as.matrix()
cell_rankings <- AUCell_buildRankings(counts_mat, plotStats = FALSE)

# Score each cell type module
auc_scores <- AUCell_calcAUC(gene_lists, cell_rankings)

# Extract AUC matrix (cell types × cells) and transpose to cells × cell types
auc_mat <- t(auc_scores@assays@data[['AUC']]%>%as.matrix())  # now cells × cell types
colnames(auc_mat) <- paste0("AUC_", colnames(auc_mat))

# Add to metadata
obj <- AddMetaData(obj, metadata = as.data.frame(auc_mat))

# ── 13. Assign predicted_annot as the cell type with highest z-scored AUC ────
auc_mat_z <- scale(auc_mat)
colnames(auc_mat_z) <- paste0("AUCz_", cell_types)
obj <- AddMetaData(obj, metadata = as.data.frame(auc_mat_z))

z_threshold <- 0.5
winner_scores <- apply(auc_mat_z, 1, max)
winner_labels <- cell_types[apply(auc_mat_z, 1, which.max)]
obj$predicted_annot <- ifelse(winner_scores >= z_threshold, winner_labels%>%unfactor(), "Unassigned")

cat("Predicted annotation distribution after z-scoring:\n")
print(table(obj$predicted_annot) / ncol(obj))

# ── 14. Visualize ─────────────────────────────────────────────────────────────
DimPlot(obj, group.by = "predicted_annot", label = TRUE)
DimPlot(obj, group.by = "final_clusters", label = TRUE)

auc_cols <- paste0("AUC_", cell_types)
FeaturePlot(obj, features = auc_cols, ncol = 3)

#### NB1 ###
obj@meta.data$Status = factor(obj@meta.data$Status, levels = c('NRM', 'M', 'D', 'E', 'NF', 'F'))

dag = obj@meta.data %>%
  group_by(individual, Status, final_clusters) %>%
  summarize(nb1 = mean(AUCz_NB.1))

ggplot(subset(dag, final_clusters == 6), aes(x = Status, y = nb1)) +
  geom_boxplot() +
  geom_point()

dag2 = obj@meta.data %>%
  group_by(individual, Status, final_clusters) %>%
  summarize(nb1 = sum(predicted_annot == 'NB.1'), tot = n()) %>%
  mutate(prop = nb1 / tot)

ggplot(subset(dag2, final_clusters == 6), aes(x = Status, y = prop)) +
  stat_summary(geom = 'bar', fun = 'mean') +
  geom_point()

dag3 = obj@meta.data %>%
  group_by(individual, Status, final_clusters) %>%
  summarize(nb2 = sum(predicted_annot == 'NB.2'), tot = n()) %>%
  mutate(prop = nb2 / tot)

ggplot(subset(dag3, final_clusters == 6), aes(x = Status, y = prop)) +
  stat_summary(geom = 'bar', fun = 'mean') +
  geom_point()

dag4 = obj@meta.data %>%
  group_by(individual, Status, final_clusters) %>%
  summarize(mn = sum(predicted_annot == 'MN'), tot = n()) %>%
  mutate(prop = mn / tot)

ggplot(subset(dag4, final_clusters == 6), aes(x = Status, y = prop)) +
  stat_summary(geom = 'bar', fun = 'mean') +
  geom_point()

dag5 = obj@meta.data %>%
  group_by(individual, Status, final_clusters) %>%
  summarize(nb = sum(predicted_annot %in% c('NB.1', 'NB.2')), tot = n()) %>%
  mutate(prop = nb / tot)

ggplot(subset(dag5, final_clusters == 6), aes(x = Status, y = prop)) +
  stat_summary(geom = 'bar', fun = 'mean') +
  geom_point()

typewise_props = data.frame()
for(i in c('MN',
           'NB.1',
           'NB.2')){
  
  dags = obj@meta.data %>%
  group_by(individual, Status, final_clusters) %>%
  summarize(cells = sum(predicted_annot %in% i), tot = n())%>%
    subset(final_clusters ==6)%>%
    subset(Status != 'NRM')
  
  mod_dat = cbind(dags$cells, dags$tot-dags$cells)
  mod = glm(mod_dat~Status, data = dags, 
            family = binomial('logit'))
  
  av_ = anova(mod, test= 'Chisq')
  av_p = av_$`Pr(>Chi)`[2]

  newd = data.frame(clust = i,
                    p = av_p)
  typewise_props = rbind(typewise_props, newd)
}

  dags = obj@meta.data %>%
  group_by(individual, Status, final_clusters) %>%
  summarize(cells = sum(predicted_annot %in% c('NB.1','NB.2')), tot = n())%>%
    subset(final_clusters ==6)%>%
    subset(Status != 'NRM')
  
  mod_dat = cbind(dags$cells, dags$tot-dags$cells)
  mod = glm(mod_dat~Status, data = dags, 
            family = binomial('logit'))
  
  av_ = anova(mod, test= 'Chisq')
  av_p = av_$`Pr(>Chi)`[2]
  