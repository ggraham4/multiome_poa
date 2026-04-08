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

# ── 11. Find transfer anchors ─────────────────────────────────────────────────
#obj_all = obj

`%notin%` = Negate(`%in%`)
obj= subset(obj,
                                final_clusters %notin% c(
                                  2, # olig
                                  13, # opc
                                  20, # leuko
                                  26, # DG
                                  11, # mg
                                  15 # eg
                                ) )


transfer_anchors <- FindTransferAnchors(
  reference            = seu,
  query                = obj,
  dims                 = 1:30,
  reference.reduction  = "pca",
  normalization.method = "LogNormalize",
  features             = intersect(VariableFeatures(seu),
                                   VariableFeatures(obj))  # shared HVGs only
)

# ── 12. Transfer annot and seurat_clusters labels ─────────────────────────────
predictions <- TransferData(
  anchorset = transfer_anchors,
  refdata   = list(
    annot           = seu$annot,
    seurat_clusters = as.character(Idents(seu))
  ),
  dims = 1:30
)

# ── 13. Add predictions to clownfish object ───────────────────────────────────
obj <- AddMetaData(obj, metadata = predictions[["annot"]][["predicted.id"]],
                   col.name = "predicted_annot")
#obj <- AddMetaData(obj, metadata = predictions$seurat_clusters,
#                   col.name = "predicted_seurat_clusters")

# ── 14. Visualize ─────────────────────────────────────────────────────────────
DimPlot(obj, group.by = "predicted_annot", label = TRUE) 
DimPlot(obj, group.by = 'final_clusters', label = T)

### Group and summarize
summary_stats = obj@meta.data %>%
  group_by(final_clusters, predicted_annot) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(final_clusters) %>%
  mutate(percent = n / sum(n) * 100) %>%
  arrange(final_clusters, desc(percent))

ggplot(summary_stats, aes(x = final_clusters,y = percent ,fill = predicted_annot))+
  geom_bar(position = 'stack', stat = 'identity')+
  scale_fill_manual(values = c('red', 'green', 'blue','black', 'purple'))

ggplot(summary_stats, aes(x =final_clusters ,y = percent ,fill = predicted_annot))+
  geom_bar(position = 'stack', stat = 'identity')+
  scale_fill_manual(values = c('red', 'green', 'green','black', 'purple'))

cluster_order <- order_df %>%
  arrange(desc(nb_sum)) %>%
  pull(final_clusters)

summary_stats$final_clusters <- factor(
  summary_stats$final_clusters,
  levels = cluster_order
)

ggplot(summary_stats, aes(x = final_clusters, y = percent, fill = predicted_annot)) +
  geom_bar(position = 'stack', stat = 'identity') +
  scale_fill_manual(values = c('red', 'green', 'green','black', 'purple'))


# whats going on in 6
sub_6 = subset(obj, final_clusters==6)
marks_sub_6 =FindAllMarkers(sub_6, group.by = 'predicted_annot', only.pos = T)



## sex differences
# I think what I can do easier is just ask if the proportion of MNs in a cluster differs

out_dat = data.frame()
clusters = unique(obj$final_clusters)
for(i in clusters){
  model_dat = obj@meta.data %>%
    group_by(individual, Status, final_clusters) %>%
    summarize(
      n_cells = n(),
      n_mn = sum(predicted_annot == 'MN')
    ) %>%
    subset(Status %in%c('M',
                        'D',
                        'F')& final_clusters == i) %>%   # <-- removed predicted_annot filter
    mutate(successes = n_mn,
           failures = n_cells - n_mn)
  
  model_mat = cbind(model_dat$successes, model_dat$failures)
  model = glm(model_mat ~ Status, data = model_dat, family = binomial('logit'))
  av_ = anova(model, test = 'Chisq')
  p = av_$`Pr(>Chi)`[2]
  newd = data.frame(cluster = i,
                    av_p = p)
  out_dat = rbind(newd, out_dat)
}

# cluster 0
model_dat = obj@meta.data %>%
  group_by(individual, Status, final_clusters) %>%   # <-- same fix here
  summarize(
    n_cells = n(),
    n_mn = sum(predicted_annot == 'MN')
  ) %>%
  subset(Status != 'NRM') %>%
  mutate(successes = n_mn,
         failures = n_cells - n_mn)

model_dat$Status = factor(model_dat$Status, levels = c('M',
                                                       'D',
                                                       'E',
                                                       'NF',
                                                       'F'))

ggplot(model_dat%>%subset(final_clusters ==0), aes(x = Status, y = (n_cells - n_mn) / n_cells)) +
  geom_boxplot()+
    geom_point()+
  labs(y = 'Immature Cells', title= 'Cluster 0')

ggplot(model_dat%>%subset(final_clusters ==7), aes(x = Status, y = (n_cells - n_mn) / n_cells)) +
  geom_boxplot()+
    geom_point()+
  labs(y = 'Immature Cells', title= 'Cluster 7')

ggplot(model_dat%>%subset(final_clusters ==10), aes(x = Status, y = (n_cells - n_mn) / n_cells))+
  geom_boxplot()+
    geom_point()+
  labs(y = 'Immature Cells', title= 'Cluster 10')

ggplot(model_dat%>%subset(final_clusters ==6), aes(x = Status, y = (n_cells - n_mn) / n_cells)) +
  geom_boxplot()+
    geom_point()+
  labs(y = 'Immature Cells', title= 'Cluster 6')

ggplot(model_dat%>%subset(final_clusters ==9), aes(x = Status, y = (n_cells - n_mn) / n_cells))+
  geom_boxplot()+
    geom_point()+
  labs(y = 'Immature Cells', title= 'Cluster 9')

ggplot(model_dat%>%subset(final_clusters ==25), aes(x = Status, y = (n_cells - n_mn) / n_cells))+
  geom_boxplot()+
    geom_point()+
  labs(y = 'Immature Cells', title= 'Cluster 25')

library(CytoTRACE)
sub_6_cyto = CytoTRACE(sub_6@assays$RNA$data%>%as.matrix())
sub_6$cyto =sub_6_cyto$CytoTRACE

mdat = sub_6@meta.data %>%
  group_by(individual, Status) %>%  
  summarize(
    mean_cyto =mean(cyto
  )) 

mdat$Status = factor(mdat$Status, levels = c('NRM', 'M',
                                                       'D',
                                                       'E',
                                                       'NF',
                                                       'F'))


  ggplot(mdat, aes(x = Status, y = mean_cyto))+
  geom_boxplot()+
    geom_point()
  
  mdat2 = sub_6@meta.data %>%
  group_by(individual, Status, predicted_annot) %>%  
  summarize(
    mean_cyto =mean(cyto
  )) 
    ggplot(mdat2, aes(x = predicted_annot, y = mean_cyto, color =Status))+
  geom_boxplot()
# ok s I think this shows that this is kinda meaningless cause wdym RGC is the lowest
    
    sum(obj@meta.data$predicted_annot=='RGC')/nrow(obj@meta.data)
        sum(obj@meta.data$predicted_annot=='RGC')/nrow(obj@meta.data)

                sum(obj@meta.data$predicted_annot%in%c('NB.1','NB.2'))/nrow(obj@meta.data)
#22%
                    