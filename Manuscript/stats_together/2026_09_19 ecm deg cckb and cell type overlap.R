library(Seurat)
library(patchwork)
library(tidyverse)
library(ggplot2)
library(Polychrome)
library(emmeans)
library(ggsignif)
library(CytoTRACE)
library(clusterProfiler)
library(AUCell)
library(lme4)
library(multcomp)
  
obj  = readRDS("~/Desktop/optimal_clustering_rna_only.rds")


sub_6 = FindSubCluster(obj, 6, graph.name = "harmony.wsnn")
Idents(sub_6) <- "sub.cluster"
sub_6 = subset(sub_6, final_clusters == 6 & Status %in% c('M','D','F'))

sub_6$Status = factor(
  sub_6$Status,
  levels = c("NRM", "M", "D", "E", "NF", "F")
)

cyto = CytoTRACE(sub_6@assays$RNA$data %>% as.matrix())
sub_6$cyto = cyto$CytoTRACE

degs_plasticity = c(
  "LOC111588913",
  "cntn4",
  "LOC111567620",
  "pcdh10b",
  "sdc2",
  "LOC111585095",
  "bcan",
  "LOC111568896"
)

sub_6$gene_pos = colSums(
  sub_6@assays$RNA$data[degs_plasticity, ]
) > 0

sub_6$cckb = sub_6@assays$RNA$data['cckb', ] > 0

DotPlot(sub_6, 
        c('gene_pos',
          'cckb'))+
  coord_flip()


table(ecm_deg = sub_6$gene_pos,cck =  sub_6$cckb)
749/(234+729)
131/(131+18)

mat = table(ecm_deg = sub_6$gene_pos,cck =  sub_6$cckb)%>%
  as.matrix()

# significantly higher
fisher.test(mat)


sub_6$drd3 = sub_6@assays$RNA$data['drd3', ] > 0

table(ecm_deg = sub_6$gene_pos,drd3 =  sub_6$drd3)

# not significant
fisher.test(table(ecm_deg = sub_6$gene_pos,drd3 =  sub_6$drd3)%>%
  as.matrix())

#significantly higher
sub_6$tacr3a = sub_6@assays$RNA$data['tacr3a', ] > 0
table(ecm_deg = sub_6$gene_pos,tacr3a =  sub_6$tacr3a)
fisher.test(table(sub_6$gene_pos, sub_6$tacr3a)%>%
  as.matrix())

# not significant
sub_6$npy7r = sub_6@assays$RNA$data['npy7r', ] > 0
table(ecm_deg = sub_6$gene_pos,npy7r =  sub_6$npy7r)
fisher.test(table(sub_6$gene_pos, sub_6$npy7r)%>%
  as.matrix())

# significantly higher
sub_6$gnrh1 = sub_6@assays$RNA$data['LOC111571064', ] > 0
table(ecm_deg = sub_6$gene_pos,gnrh1 =  sub_6$gnrh1)
fisher.test(table(sub_6$gene_pos, sub_6$gnrh1)%>%
  as.matrix())

# significantly higher #
sub_6$pgr = sub_6@assays$RNA$data['pgr', ] > 0
table(ecm_deg = sub_6$gene_pos,pgr =  sub_6$pgr)
fisher.test(table(sub_6$gene_pos, sub_6$pgr)%>%
  as.matrix())




# significantly higher#
sub_6$esr2b = sub_6@assays$RNA$data['esr2b', ] > 0
table(ecm_deg = sub_6$gene_pos,esr2b =  sub_6$esr2b)
fisher.test(table(sub_6$gene_pos, sub_6$esr2b)%>%
  as.matrix())

#significantly higher
sub_6$ar = sub_6@assays$RNA$data['ar', ] > 0
table(ecm_deg = sub_6$gene_pos,ar =  sub_6$ar)
fisher.test(table(sub_6$gene_pos, sub_6$ar)%>%
  as.matrix())

# significantly higher
sub_6$arlike = sub_6@assays$RNA$data['LOC111568069', ] > 0
table(ecm_deg = sub_6$gene_pos,ar =  sub_6$arlike)
fisher.test(table(sub_6$gene_pos, sub_6$arlike)%>%
  as.matrix())

### are they the same cells?




library(tidyverse)
library(ggplot2)

# ============================================================
# Genes to test
# ============================================================

genes <- c(
  "cckb",
  "drd3",
  "tacr3a",
  "npy7r",
  "gnrh1",
  "pgr",
  "esr2b",
  "ar",
  "arlike"
)

# Actual RNA gene names
gene_names <- c(
  cckb   = "cckb",
  drd3   = "drd3",
  tacr3a = "tacr3a",
  npy7r  = "npy7r",
  gnrh1  = "LOC111571064",
  pgr    = "pgr",
  esr2b  = "esr2b",
  ar     = "ar",
  arlike = "LOC111568069"
)


# ============================================================
# Calculate proportions and Fisher's exact test
# ============================================================

plot_data <- map_dfr(genes, function(g) {
  
  # ECM+ cells
  ecm <- sub_6$gene_pos
  
  # Gene+ cells
  gene_pos <- sub_6@assays$RNA$data[gene_names[g], ] > 0
  
  # 2 x 2 contingency table
  tab <- table(
    ECM = ecm,
    Gene = gene_pos
  )
  
  # One-sided Fisher test:
  # Is ECM+ enriched among Gene+ cells?
  fisher <- fisher.test(
    tab,
    alternative = "greater"
  )
  
  # Proportion of ALL cells that are ECM+
  prop_ecm_all <- mean(ecm)
  
  # Proportion of Gene+ cells that are ECM+
  prop_ecm_gene <- sum(ecm & gene_pos) / sum(gene_pos)
  
  tibble(
    gene = g,
    group = c("All cells", "Gene+ cells"),
    proportion = c(
      prop_ecm_all,
      prop_ecm_gene
    ),
    p_value = fisher$p.value,
    odds_ratio = unname(fisher$estimate),
    n_gene_pos = sum(gene_pos),
    n_ecm_gene_pos = sum(ecm & gene_pos)
  )
})


# ============================================================
# Significance labels
# ============================================================

sig_data <- plot_data %>%
  group_by(gene) %>%
  summarise(
    p_value = first(p_value),
    max_prop = max(proportion, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    significance = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01  ~ "**",
      p_value < 0.05  ~ "*",
      TRUE ~ ""
    )
  )


# ============================================================
# Plot
# ============================================================

p <- ggplot(
  plot_data,
  aes(
    x = group,
    y = proportion,
    fill = group
  )
) +
  
  geom_col(
    width = 0.65,
    color = "black",
    linewidth = 0.3
  ) +
  
  # Significance stars
  geom_text(
    data = sig_data,
    aes(
      x = "Gene+ cells",
      y = max_prop + 0.05,
      label = significance
    ),
    inherit.aes = FALSE,
    size = 6
  ) +
  
  facet_wrap(
    ~ gene,
    nrow = 1
  ) +
  
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1),
    limits = c(0, 1.05),
    expand = expansion(mult = c(0, 0.02))
  ) +
  
  scale_fill_manual(
    values = c(
      "All cells" = "grey70",
      "Gene+ cells" = "steelblue"
    )
  ) +
  
  labs(
    x = NULL,
    y = "Proportion of cells that are ECM+"
  ) +
  
  theme_classic(base_size = 13) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(
      face = "bold",
      size = 12
    ),
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    ),
    legend.position = "none",
    panel.spacing = unit(1, "lines")
  )

p


# ============================================================
# Print Fisher test results
# ============================================================

results <- plot_data %>%
  distinct(
    gene,
    p_value,
    odds_ratio,
    n_gene_pos,
    n_ecm_gene_pos
  ) %>%
  arrange(p_value)

results

##### matrix #####
# ============================================================
# Pairwise overlap matrix
# Each cell = proportion of Gene A+ cells that are also Gene B+
# ============================================================

genes <- c(
  "cckb",
  "drd3",
  "tacr3a",
  "npy7r",
  "gnrh1",
  "pgr",
  "esr2b",
  "ar",
  "arlike"
)

# Actual RNA gene names
gene_names <- c(
  cckb   = "cckb",
  drd3   = "drd3",
  tacr3a = "tacr3a",
  npy7r  = "npy7r",
  gnrh1  = "LOC111571064",
  pgr    = "pgr",
  esr2b  = "esr2b",
  ar     = "ar",
  arlike = "LOC111568069"
)

# ------------------------------------------------------------
# Make binary matrix: cells x genes
# ------------------------------------------------------------

gene_matrix <- sapply(genes, function(g) {
  sub_6@assays$RNA$data[gene_names[g], ] > 0
})

colnames(gene_matrix) <- genes


# ------------------------------------------------------------
# Calculate pairwise overlap
# Rows = Gene A+
# Columns = Gene B+
#
# Value = proportion of Gene A+ cells that are also Gene B+
# ------------------------------------------------------------

overlap_matrix <- matrix(
  NA,
  nrow = length(genes),
  ncol = length(genes),
  dimnames = list(genes, genes)
)

for (i in seq_along(genes)) {
  
  for (j in seq_along(genes)) {
    
    gene_A <- gene_matrix[, i]
    gene_B <- gene_matrix[, j]
    
    overlap_matrix[i, j] <- sum(gene_A & gene_B) / sum(gene_A)
  }
}


# View matrix
round(overlap_matrix, 3)


# ------------------------------------------------------------
# Heatmap
# ------------------------------------------------------------

overlap_df <- as.data.frame(overlap_matrix) %>%
  rownames_to_column("Gene_A") %>%
  pivot_longer(
    cols = -Gene_A,
    names_to = "Gene_B",
    values_to = "Proportion"
  )

overlap_df$Gene_A <- factor(
  overlap_df$Gene_A,
  levels = rev(genes)
)

overlap_df$Gene_B <- factor(
  overlap_df$Gene_B,
  levels = genes
)


ggplot(
  overlap_df,
  aes(
    x = Gene_B,
    y = Gene_A,
    fill = Proportion
  )
) +
  geom_tile(
    color = "white",
    linewidth = 0.3
  ) +
  geom_text(
    aes(label = sprintf("%.2f", Proportion)),
    size = 3
  ) +
  scale_fill_gradient(
    low = "white",
    high = "steelblue",
    limits = c(0, 1),
    labels = scales::percent_format(accuracy = 1)
  ) +
  labs(
    x = "Gene B",
    y = "Gene A",
    fill = "Proportion"
  ) +
  coord_fixed() +
  theme_classic(base_size = 13) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  )+
  labs(title = 'Proportion of Gene A+ cells that are also Gene B+')



