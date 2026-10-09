#### Reproduces every statistic in Appendix S1 ####
# Run from the project root (the folder containing "Measures/", "DEG Outputs/", "Manuscript/", ...).
# Output: "Manuscript/Supplementary Tables/Appendix S1. Statistics.xlsx" (one sheet per figure).
#
# Sections 1-16 are the analyses; section 17 assembles the workbook.
# Phase codes: M male, I intermediate, IP intermediate partner, LI late intermediate,
# LIP late intermediate partner, F female, NF new female, NM new male.

#### 0. Setup ####
library(Seurat)
library(tidyverse)
library(emmeans)
library(car)
library(lme4)
library(CytoTRACE)
library(AUCell)
library(openxlsx)
set.seed(0)
options(pillar.sigfig = 5)

# setwd("/path/to/multiome_poa")   # uncomment and edit if not already in the project root

status_to_phase <- c(M = "M", D = "I", S = "IP", E = "LI", EP = "LIP",
                     F = "F", NF = "NF", NM = "NM", NRM = "NRM")
mif <- c("M", "I", "F")
all_phases <- c("M", "I", "IP", "LI", "LIP", "F", "NF", "NM")

add_phase <- function(d, keep = mif) {
  d %>%
    mutate(Phase = unname(status_to_phase[as.character(Status)])) %>%
    filter(Phase %in% keep) %>%
    mutate(Phase = factor(Phase, levels = keep))
}

model_label <- function(fit) {
  fn <- if (inherits(fit, "merMod")) "lmer" else if (inherits(fit, "glm")) "glm" else "lm"
  paste0(fn, "(", deparse1(formula(fit)), if (fn == "glm") ", family = binomial", ")")
}

# Type-III test per term + uncorrected emmeans pairwise p-values on the grouping term
# contrasts = NULL returns every pairwise contrast
stat_table <- function(models, figure, group = "Phase", contrasts = c("M-F", "M-I", "I-F")) {
  imap_dfr(models, function(fit, trait) {
    a <- as.data.frame(Anova(fit, type = 3))
    a <- a[!rownames(a) %in% c("(Intercept)", "Residuals"), , drop = FALSE]
    stat <- intersect(c("F value", "LR Chisq", "Chisq"), names(a))[1]
    pw <- as.data.frame(pairs(emmeans(fit, group), adjust = "none"))
    lab <- gsub(" ", "", as.character(pw$contrast))
    want <- if (is.null(contrasts)) lab else contrasts
    pw_p <- map_dbl(want, function(x) {
      i <- which(lab %in% c(x, paste(rev(strsplit(x, "-")[[1]]), collapse = "-")))
      if (length(i)) pw$p.value[i[1]] else NA_real_
    })
    out <- tibble(Trait = trait, Model = model_label(fit), Figure = figure, Term = rownames(a),
                  !!stat := a[[stat]], df = a$Df, `Type-III p-value` = a[[ncol(a)]])
    if (!inherits(fit, "glm") && inherits(fit, "lm")) out$df_resid <- df.residual(fit)
    for (k in seq_along(want)) out[[want[k]]] <- ifelse(out$Term == group, pw_p[k], NA_real_)
    out
  })
}

mean_se <- function(data, traits, group = "Phase") {
  data %>%
    dplyr::select(all_of(c(group, traits))) %>%
    pivot_longer(-all_of(group), names_to = "trait") %>%
    filter(!is.na(value)) %>%
    group_by(trait, .data[[group]]) %>%
    summarize(mean = mean(value), se = sd(value) / sqrt(n()), n = n(), .groups = "drop")
}

# Per-individual cell counts in each cluster
cluster_counts <- function(o, cluster_col) {
  o@meta.data %>%
    count(individual, Status, cluster = as.character(.data[[cluster_col]]), name = "ncells") %>%
    group_by(individual) %>%
    mutate(total_cells = sum(ncells)) %>%
    ungroup() %>%
    add_phase() %>%
    mutate(prop = ncells / total_cells)
}

prop_glm <- function(d) glm(cbind(ncells, total_cells - ncells) ~ Phase, family = binomial, data = d)

go_terms <- readRDS("Function Scripts/Dependencies/Term2gene_clown_go2.rds") %>%
  left_join(readRDS("Function Scripts/Dependencies/Term2name.rds"), by = "go_id")

# AUCell GO module score, averaged per individual
go_module <- function(term, o) {
  set.seed(1)
  genes <- intersect(go_terms$aocellaris_name[go_terms$go_id == term], rownames(o))
  message(length(genes), " genes found for: ", unique(go_terms$go_name[go_terms$go_id == term]))
  ranks <- AUCell_buildRankings(o@assays$RNA$data, plotStats = FALSE, verbose = FALSE)
  o$score <- as.numeric(getAUC(AUCell_calcAUC(setNames(list(genes), term), ranks, verbose = FALSE))[1, ])
  o@meta.data %>%
    group_by(individual, Status) %>%
    summarize(score = mean(score), .groups = "drop") %>%
    add_phase()
}

measures <- read.csv("Measures/2025_12_26 all_data.csv") %>%
  mutate(across(any_of(c("Behaviors_Day_2", "Time_Day_2", "Log_11KT", "Percent_Testicular",
                         "Percent_Ovarian", "Log10_Volume", "length_final_cm", "mass_final_cm",
                         "Change_Mass", "Change_Length", "Log10_Testicular_Estimate",
                         "Log10_Ovarian_Estimate")), as.numeric)) %>%
  add_phase(keep = all_phases)


#### 1. Fig. 1B-C (Table 1) ####
fig1_models <- list(
  `# parental acts` = lm(Behaviors_Day_2 ~ Phase, data = measures),
  `time in nest`    = lm(Time_Day_2 ~ Phase, data = measures),
  `log (11-KT)`     = lm(Log_11KT ~ Phase, data = measures),
  `% testicular`    = lm(Percent_Testicular ~ Phase, data = measures),
  `% ovarian`       = lm(Percent_Ovarian ~ Phase, data = measures),
  `volume`          = lm(Log10_Volume ~ Phase + length_final_cm, data = measures)
)
fig1_table <- stat_table(fig1_models, "Fig. 1B-C")

# Length-corrected gonad volume (residual of log10 volume on final length)
vol_len <- lm(Log10_Volume ~ length_final_cm, data = measures)
measures$resid_volume <- residuals(vol_len)[rownames(measures)]

fig1_summary <- mean_se(measures, c("Behaviors_Day_2", "Time_Day_2", "Log_11KT",
                                    "Percent_Testicular", "Percent_Ovarian", "resid_volume"))


#### 2. Fig. S1A-F ####
figS1ab_table <- stat_table(list(
  `body mass`   = lm(mass_final_cm ~ Phase, data = measures),
  `body length` = lm(length_final_cm ~ Phase, data = measures)
), "Fig. S1A-B")
figS1ab_summary <- mean_se(measures, c("mass_final_cm", "length_final_cm"))

dat_experiment <- filter(measures, Condition == "Experiment")
figS1cf_table <- stat_table(list(
  `change in mass`      = lm(Change_Mass ~ Phase, data = dat_experiment),
  `change in length`    = lm(Change_Length ~ Phase, data = dat_experiment),
  `log10 testis volume` = lm(Log10_Testicular_Estimate ~ Phase, data = measures),
  `log10 ovary volume`  = lm(Log10_Ovarian_Estimate ~ Phase, data = measures)
), "Fig. S1C-F")
figS1cf_summary <- bind_rows(
  mean_se(dat_experiment, c("Change_Mass", "Change_Length")),
  mean_se(measures, c("Log10_Testicular_Estimate", "Log10_Ovarian_Estimate"))
)


#### 3. snMultiome objects ####
obj <- readRDS("~/Desktop/optimal_clustering_rna_only.rds")

# non-neuronal: 1 RGC, 2 OL, 11 MG, 13 OPC, 15 EC, 20 Leuko, 26 DG
neurons_only <- obj[, !obj$res0.8_50nn_40PC_45LSI %in% c(1, 2, 11, 13, 15, 20, 26)]

obj <- FindSubCluster(obj, 1, "harmony.wsnn", resolution = 0.2, subcluster.name = "sub_res0.2")
sub_1 <- subset(obj, final_clusters == 1)
Idents(sub_1) <- "sub_res0.2"
sub_1$cyto <- CytoTRACE(as.matrix(sub_1@assays$RNA$data))$CytoTRACE

sub_6 <- subset(obj, final_clusters == 6 & Status %in% c("M", "D", "F"))
sub_6$cyto <- CytoTRACE(as.matrix(sub_6@assays$RNA$data))$CytoTRACE
expr6 <- sub_6@assays$RNA$data


#### 4. Fig. 3A-B, S6A-B: cell-type abundance ####
counts_all <- cluster_counts(obj, "res0.8_50nn_40PC_45LSI")
counts_neuron <- cluster_counts(neurons_only, "res0.8_50nn_40PC_45LSI")

figS6a_table <- map(split(counts_all, counts_all$cluster), prop_glm) %>%
  stat_table("Fig. S6A") %>%
  mutate(q.value = p.adjust(`Type-III p-value`, "fdr"))
figS6b_table <- map(split(counts_neuron, counts_neuron$cluster), prop_glm) %>%
  stat_table("Fig. S6B") %>%
  mutate(q.value = p.adjust(`Type-III p-value`, "fdr"))
figS6a_summary <- map_dfr(split(counts_all, counts_all$cluster), mean_se, "prop", .id = "cluster")
figS6b_summary <- map_dfr(split(counts_neuron, counts_neuron$cluster), mean_se, "prop", .id = "cluster")


#### 5. Fig. 3C-D: DEG and DAR enrichment ####
chisq_enrichment <- function(d) {
  chi <- chisq.test(d$n)
  d %>%
    mutate(expected = chi$expected,
           enrichment = n / expected,
           residual = (n - expected) / sqrt(expected),
           p_value = 2 * pnorm(-abs(residual)),
           p_adj = p.adjust(p_value, "BH"),
           signif = p_adj < 0.05 & enrichment > 1)
}

fig3c_deg <- read.csv("DEG Outputs/FINAL degs classified w singular.csv") %>%
  count(cluster) %>%
  chisq_enrichment()
fig3d_dar <- read.csv("Collaboration/all_clusters_DARs_peak_level_classified_with_support.csv") %>%
  count(cluster = cluster_id) %>%
  chisq_enrichment()
fig3c_chi <- chisq.test(fig3c_deg$n)
fig3d_chi <- chisq.test(fig3d_dar$n)


#### 6. Fig. 3E: classifier scores with QC ####
# QC: keep clusters whose leave-one-out M/F validation gives mean P(female) < 0.10 for males and > 0.90 for females
mf_val <- read.csv("Sex Classifier and Linearity/progress_classifier/Including DARs/logistic_mf_validation_07_13_2026.csv")
validation <- mf_val %>%
  group_by(cluster, status) %>%
  summarize(mean_prob = mean(proba), .groups = "drop") %>%
  pivot_wider(names_from = status, values_from = mean_prob, names_prefix = "val_")
qc_pass <- validation$cluster[validation$val_m < 0.10 & validation$val_f > 0.90]

classifier_scores <- read.csv("Manuscript/Supplementary Tables/classifier_scores.csv") %>%
  filter(status == "D")   # intermediates only

fig3e_by_cluster <- classifier_scores %>%
  group_by(cluster) %>%
  summarize(n_features = first(n_degs), mean_score = mean(prediction),
            se_score = sd(prediction) / sqrt(n()), n = n(), .groups = "drop") %>%
  left_join(validation, by = "cluster") %>%
  mutate(qc_pass = cluster %in% qc_pass)

fig3e_overall <- classifier_scores %>%
  filter(cluster %in% qc_pass) %>%
  summarize(mean_score = mean(prediction), se_score = sd(prediction) / sqrt(n()),
            n_rows = n(), n_clusters = n_distinct(cluster))


#### 7. Fig. 4D-E: RGC subcluster abundance ####
counts_rgc <- cluster_counts(sub_1, "sub_res0.2")
fig4de_table <- stat_table(list(
  `1_NSC (1_1)` = prop_glm(filter(counts_rgc, cluster == "1_1")),
  `1_NP (1_0)`  = prop_glm(filter(counts_rgc, cluster == "1_0"))
), "Fig. 4D-E")
fig4de_summary <- map_dfr(split(counts_rgc, counts_rgc$cluster), mean_se, "prop", .id = "cluster")


#### 8. Fig. 4F: neuron differentiation module in 1_NP ####
neuron_diff <- go_module("GO:0030182", subset(sub_1, sub_res0.2 == "1_0"))
fig4f_table <- stat_table(list(`neuron differentiation (GO:0030182), 1_NP` = lm(score ~ Phase, data = neuron_diff)), "Fig. 4F")
fig4f_summary <- mean_se(neuron_diff, "score")


#### 9. Fig. S7B: CytoTRACE by 1_RGC subcluster ####
cyto_sub1 <- sub_1@meta.data %>%
  filter(Status != "NRM") %>%
  group_by(individual, Status, sub_res0.2) %>%
  summarize(mean_cyto = mean(cyto), .groups = "drop") %>%
  mutate(sub_res0.2 = factor(sub_res0.2))
figS7b_table <- stat_table(list(`CytoTRACE by 1_RGC subcluster` = lm(mean_cyto ~ sub_res0.2, data = cyto_sub1)),
                           "Fig. S7B", group = "sub_res0.2", contrasts = NULL)
figS7b_summary <- mean_se(cyto_sub1, "mean_cyto", group = "sub_res0.2")


#### 10. Fig. 5B: CytoTRACE in 6_POA_Mixed ####
cyto_6 <- sub_6@meta.data %>%
  group_by(individual, Status) %>%
  summarize(cyto = mean(cyto), .groups = "drop") %>%
  add_phase()
fig5b_table <- stat_table(list(`CytoTRACE, 6_POA_Mixed` = lm(cyto ~ Phase, data = cyto_6)), "Fig. 5B")
fig5b_summary <- mean_se(cyto_6, "cyto")


#### 11. Fig. 5C-D, 5G: GO modules in 6_POA_Mixed ####
brain_dev <- go_module("GO:0007420", sub_6)
axon_guid <- go_module("GO:0008046", sub_6)
syn_plas  <- go_module("GO:0048167", sub_6)

fig5cd_table <- stat_table(list(
  `brain development (GO:0007420)`               = lm(score ~ Phase, data = brain_dev),
  `axon guidance receptor activity (GO:0008046)` = lm(score ~ Phase, data = axon_guid)
), "Fig. 5C-D")
fig5cd_summary <- bind_rows(`brain development (GO:0007420)` = mean_se(brain_dev, "score"),
                            `axon guidance receptor activity (GO:0008046)` = mean_se(axon_guid, "score"),
                            .id = "module")
syn_plas_table <- stat_table(list(`regulation of synaptic plasticity (GO:0048167)` = lm(score ~ Phase, data = syn_plas)), "Fig. 5G")
syn_plas_summary <- mean_se(syn_plas, "score")


#### 12. Fig. 6A-C, S10A-E: gene+ populations in 6_POA_Mixed ####
genes_6 <- c(drd3 = "drd3", tacr3a = "tacr3a", cckb = "cckb", npy7r = "npy7r", nmbr = "nmbr",
             pgr = "pgr", gnrh1 = "LOC111571064", ar_like = "LOC111568069")

# Proportion of gene+ cells
gene_counts <- imap_dfr(genes_6, function(g, nm) {
  sub_6@meta.data %>%
    mutate(pos = expr6[g, ] > 0) %>%
    group_by(individual, Status) %>%
    summarize(n_pos = sum(pos), n_cells = n(), .groups = "drop") %>%
    mutate(gene = nm)
}) %>%
  add_phase() %>%
  mutate(prop = n_pos / n_cells)
gene_split <- split(gene_counts, factor(gene_counts$gene, levels = names(genes_6)))

fig6_prop_table <- map(gene_split, ~ glm(cbind(n_pos, n_cells - n_pos) ~ Phase,
                                         family = binomial, data = .x)) %>%
  stat_table("Fig. 6A-C, S10A-E")
fig6_prop_summary <- map_dfr(gene_split, mean_se, "prop", .id = "gene")

# CytoTRACE of gene+ cells
cyto_genes <- c("drd3", "tacr3a", "cckb", "pgr", "ar_like")
cyto_pos <- map_dfr(cyto_genes, function(g) {
  sub_6@meta.data[expr6[genes_6[g], ] > 0, ] %>% mutate(gene = g)
}) %>%
  add_phase()
cyto_split <- split(cyto_pos, factor(cyto_pos$gene, levels = cyto_genes))

fig6_cyto_table <- map(cyto_split, ~ lmer(cyto ~ nCount_RNA + Phase + (1 | individual), data = .x)) %>%
  stat_table("Fig. 6A-C, S10A-E")
cyto_pos_ind <- cyto_pos %>%
  group_by(gene, individual, Phase) %>%
  summarize(cyto = mean(cyto), .groups = "drop")
fig6_cyto_summary <- map_dfr(split(cyto_pos_ind, cyto_pos_ind$gene), mean_se, "cyto", .id = "gene")


#### 13. Fig. 6D-F: steroid receptor-associated expression ####
steroid <- read.csv("Manuscript/updatedcluster_6_steroid_receptor_SPECIFIC_PROMOTER_ZSCORES.csv") %>%
  group_by(individual, Phase = group) %>%
  summarize(mean_esr2b = mean(ESR2B_score), mean_ar = mean(AR_score),
            mean_pgr = mean(PGR_score), .groups = "drop") %>%
  filter(Phase %in% mif) %>%
  mutate(Phase = factor(Phase, levels = mif))

fig6df_table <- stat_table(list(
  `ESR2B-associated expression` = lm(mean_esr2b ~ Phase, data = steroid),
  `AR-associated expression`    = lm(mean_ar ~ Phase, data = steroid),
  `PGR-associated expression`   = lm(mean_pgr ~ Phase, data = steroid)
), "Fig. 6D-F")
fig6df_summary <- mean_se(steroid, c("mean_esr2b", "mean_ar", "mean_pgr"))


#### 14. Named DEGs: negative binomial GLMM ####
deg_genes <- list(
  `0` = c(igf2bp1 = "igf2bp1"),
  `1` = c(cyp19a1b = "LOC111577263"),
  `6` = c(drd3 = "drd3", tacr3a = "tacr3a", npy7r = "npy7r", nmbr = "nmbr",
          cckb = "cckb", pgr = "pgr")   # gnrh1 is analysed with lm (section 16)
)
cluster_names <- c(`0` = "0_Spall_Mixed", `1` = "1_RGC", `6` = "6_POA_Mixed")

deg_stats <- function(obj, cluster, genes, clustering = "res0.8_50nn_40PC_45LSI") {
  keep <- obj@meta.data[[clustering]] == cluster & obj$Status %in% c("M", "D", "F")
  counts <- obj@assays$RNA$counts[, keep]
  counts <- as.matrix(counts[Matrix::rowSums(counts) > 0, ])
  meta <- data.frame(subject = factor(obj$individual[keep]),
                     Status = factor(as.character(obj$Status[keep])))

  # size factors and dispersions from all genes in the cluster, as in the full DEG run
  sf <- scran::calculateSumFactors(counts, max.cluster.size = ncol(counts), positive = TRUE)
  disp <- glmGamPoi::glm_gp(counts, col_data = meta["Status"], size_factors = sf,
                            design = ~ Status, on_disk = FALSE)$overdispersion_shrinkage_list$ql_disp_estimate
  meta$log_sf <- log(sf)

  genes <- genes[genes %in% rownames(counts)]
  lv <- levels(meta$Status)
  w <- function(a, b) as.numeric(lv == a) - as.numeric(lv == b)

  imap_dfr(genes, function(g, nm) {
    d <- mutate(meta, outcome = counts[g, ])
    m <- glmer(outcome ~ Status + (1 | subject), data = d, offset = log_sf,
               family = MASS::negative.binomial(theta = 1 / disp[match(g, rownames(counts))]))
    av <- Anova(m, type = 3)
    ct <- as.data.frame(contrast(emmeans(m, ~ Status),
                                 list(`I/M` = w("D", "M"), `F/M` = w("F", "M"), `F/I` = w("F", "D")),
                                 type = "response"))
    out <- tibble(cluster = cluster, Cluster = cluster_names[cluster], Gene = nm, Symbol = g,
                  `Wald Chisq` = av["Status", "Chisq"], df = av["Status", "Df"],
                  `Type-III p-value` = av["Status", "Pr(>Chisq)"])
    for (k in seq_len(nrow(ct))) {
      out[[paste(ct$contrast[k], "fold")]] <- ct$ratio[k]
      out[[paste(ct$contrast[k], "p")]] <- ct$p.value[k]
    }
    out$converged <- is.null(m@optinfo$conv$lme4$messages)
    out$singular <- isSingular(m)
    out
  })
}

# q-values come from the full-cluster DEG runs (FDR across all genes in the cluster)
saved_q <- map_dfr(names(deg_genes), function(cl) {
  read.csv(paste0("DEG Outputs/05_11_2025 Neg Bin w Doms New_clustering/cluster_", cl, ".csv")) %>%
    transmute(cluster = cl, Symbol = gene, q.value = suppressWarnings(as.numeric(av_q.value)))
})

deg_table <- imap_dfr(deg_genes, ~ deg_stats(obj, .y, .x)) %>%
  left_join(saved_q, by = c("cluster", "Symbol"))


#### 15. Named DEGs: mean normalized expression per individual ####
expr_genes <- deg_genes
expr_genes$`6` <- c(expr_genes$`6`, gnrh1 = "LOC111571064")

expr_by_ind <- imap_dfr(expr_genes, function(genes, cl) {
  keep <- obj$res0.8_50nn_40PC_45LSI == cl & obj$Status %in% c("M", "D", "F")
  genes <- genes[genes %in% rownames(obj)]
  ex <- obj@assays$RNA$data[genes, keep, drop = FALSE]
  imap_dfr(genes, function(g, nm) {
    tibble(individual = obj$individual[keep], Status = as.character(obj$Status[keep]), expr = ex[g, ]) %>%
      group_by(individual, Status) %>%
      summarize(expr = mean(expr), .groups = "drop") %>%
      add_phase() %>%
      mutate(Cluster = cluster_names[cl], Gene = nm)
  })
})

expr_summary <- expr_by_ind %>%
  group_by(Cluster, Gene, Phase) %>%
  summarize(mean = mean(expr), se = sd(expr) / sqrt(n()), n = n(), .groups = "drop")


#### 16. gnrh1 (Fig. S10C): lm on per-individual mean expression ####
gnrh1_table <- stat_table(list(`gnrh1 expression, 6_POA_Mixed` = lm(expr ~ Phase, data = filter(expr_by_ind, Gene == "gnrh1"))),
                          "Fig. S10C")


#### 17. Write Appendix S1 ####
chi2 <- "χ²"
trait_labels <- c(Behaviors_Day_2 = "# parental acts", Time_Day_2 = "Time in nest", Log_11KT = "log 11-KT",
                  Percent_Testicular = "Proportion testicular tissue", Percent_Ovarian = "Proportion ovarian tissue",
                  resid_volume = "Length-corrected log10 gonad volume (residual)",
                  mass_final_cm = "Final body mass", length_final_cm = "Final body length",
                  Change_Mass = "Final mass (% of initial)", Change_Length = "Final length (% of initial)",
                  Log10_Testicular_Estimate = "log10 testis volume", Log10_Ovarian_Estimate = "log10 ovary volume",
                  mean_esr2b = "ESR2B-associated expression", mean_ar = "AR-associated expression",
                  mean_pgr = "PGR-associated expression", score = "Module score", cyto = "CytoTRACE",
                  mean_cyto = "CytoTRACE")
gene_labels <- c(ar_like = "ar-like (LOC111568069)", gnrh1 = "gnrh1 (LOC111571064)")
deg_note <- "glmer negative binomial: counts ~ Phase + (1 | individual), offset = log(size factor); folds are model-based ratios."

relab <- function(x, map) { x <- as.character(x); i <- x %in% names(map); x[i] <- map[x[i]]; x }

# --- table builders ---
tests_tbl <- function(d, stat_label) {
  stat <- intersect(c("F value", "LR Chisq", "Chisq", "Wald Chisq"), names(d))[1]
  out <- tibble(Measure = d$Trait, Model = d$Model, Term = d$Term, !!stat_label := d[[stat]])
  if ("df_resid" %in% names(d)) {
    out$`df (num)` <- d$df
    out$`df (resid)` <- d$df_resid
  } else {
    out$df <- d$df
  }
  out$`p-value (Type III)` <- d$`Type-III p-value`
  if ("q.value" %in% names(d)) out$`q-value (FDR)` <- d$q.value
  for (k in c("M-F", "M-I", "I-F")) if (k %in% names(d)) out[[paste0("p (", sub("-", " vs ", k), ")")]] <- d[[k]]
  out
}

means_long <- function(s, labels = trait_labels, group = "Phase", group_name = "Phase") {
  s %>%
    mutate(Measure = relab(trait, labels)) %>%
    arrange(Measure, .data[[group]]) %>%
    transmute(Measure, !!group_name := as.character(.data[[group]]), Mean = mean, SE = se, n)
}

wide_mif <- function(s, id, idname) {
  s %>%
    filter(Phase %in% mif) %>%
    mutate(Phase = factor(Phase, levels = mif)) %>%
    dplyr::select(all_of(id), Phase, mean, se, n) %>%
    pivot_wider(names_from = Phase, values_from = c(mean, se, n), names_glue = "{.value} {Phase}") %>%
    dplyr::select(all_of(id), all_of(unlist(lapply(mif, function(p) paste(c("mean", "se", "n"), p))))) %>%
    rename_with(~ sub("^se ", "SE ", sub("^mean ", "Mean ", .x))) %>%
    rename(!!idname := all_of(id))
}

abund_tbl <- function(tab, summ) {
  t <- tests_tbl(tab, paste("LR", chi2)) %>% dplyr::select(-Model, -Term) %>% rename(Cluster = Measure)
  w <- wide_mif(summ, "cluster", "Cluster") %>%
    rename_with(~ sub("^SE ", "SE prop. ", sub("^Mean ", "Mean prop. ", .x)))
  left_join(t, w, by = "Cluster") %>% arrange(as.numeric(Cluster))
}

enrich_tbl <- function(d) {
  tibble(Cluster = as.character(d$cluster), Observed = d$n, Expected = d$expected,
         `Enrichment (obs/exp)` = d$enrichment, `Std. residual` = d$residual,
         p = d$p_value, `q (BH)` = d$p_adj, `Significantly enriched` = d$signif)
}

deg_tbl <- function(rows) {
  tibble(Cluster = rows$Cluster, Gene = rows$Gene, Symbol = rows$Symbol,
         !!paste("Wald", chi2) := rows$`Wald Chisq`, df = rows$df,
         `p-value (Type III)` = rows$`Type-III p-value`, `q-value (FDR, full DEG run)` = rows$q.value,
         `Fold I/M` = rows$`I/M fold`, `p (I vs M)` = rows$`I/M p`,
         `Fold F/M` = rows$`F/M fold`, `p (F vs M)` = rows$`F/M p`,
         `Fold F/I` = rows$`F/I fold`, `p (F vs I)` = rows$`F/I p`,
         Converged = rows$converged, `Singular fit` = rows$singular)
}

expr_tbl <- function(genes, cluster) {
  expr_summary %>%
    filter(Cluster == cluster, Gene %in% genes) %>%
    arrange(Gene, Phase) %>%
    transmute(Cluster, Gene, Phase = as.character(Phase), `Mean expression` = mean, SE = se, n)
}

# --- workbook writer: stacked, titled tables on each sheet ---
wb <- createWorkbook()
modifyBaseFont(wb, fontName = "Arial", fontSize = 10)
st_sheet <- createStyle(fontName = "Arial", fontSize = 14, textDecoration = "bold")
st_title <- createStyle(fontName = "Arial", fontSize = 12, textDecoration = "bold")
st_bold  <- createStyle(fontName = "Arial", fontSize = 10, textDecoration = "bold")
st_note  <- createStyle(fontName = "Arial", fontSize = 9, textDecoration = "italic", fontColour = "#555555")
st_hdr   <- createStyle(fontName = "Arial", fontSize = 10, textDecoration = "bold", fontColour = "#FFFFFF",
                        fgFill = "#3B5B7A", halign = "center", valign = "center", wrapText = TRUE)
st_body  <- createStyle(border = "bottom", borderColour = "#BFBFBF")
st_hl    <- createStyle(fgFill = "#FFF2CC")
st_dec   <- createStyle(numFmt = "0.0000")
st_int   <- createStyle(numFmt = "0")
st_sci   <- createStyle(numFmt = "0.00E+00")
next_row <- list()

new_sheet <- function(name, title, note = NULL) {
  addWorksheet(wb, name)
  writeData(wb, name, title, startRow = 1); addStyle(wb, name, st_sheet, rows = 1, cols = 1)
  r <- 2
  if (!is.null(note)) { writeData(wb, name, note, startRow = 2); addStyle(wb, name, st_note, rows = 2, cols = 1); r <- 3 }
  next_row[[name]] <<- r + 1
}

add_text <- function(name, txt, style = st_bold) {
  r <- next_row[[name]]
  writeData(wb, name, txt, startRow = r); addStyle(wb, name, style, rows = r, cols = 1)
  next_row[[name]] <<- r + 1
}

add_table <- function(name, title, df, note = NULL, highlight = NULL) {
  df <- as.data.frame(df, check.names = FALSE)
  r <- next_row[[name]]
  writeData(wb, name, title, startRow = r); addStyle(wb, name, st_title, rows = r, cols = 1); r <- r + 1
  if (!is.null(note)) { writeData(wb, name, note, startRow = r); addStyle(wb, name, st_note, rows = r, cols = 1); r <- r + 1 }
  writeData(wb, name, df, startRow = r, headerStyle = st_hdr, keepNA = FALSE)
  rows <- r + seq_len(nrow(df))
  addStyle(wb, name, st_body, rows = rows, cols = seq_along(df), gridExpand = TRUE, stack = TRUE)
  for (j in seq_along(df)) {
    x <- df[[j]]
    if (!is.numeric(x)) next
    nm <- names(df)[j]
    if (grepl("^(n|df|Observed)( |$)|^df \\(", nm)) {
      addStyle(wb, name, st_int, rows = rows, cols = j, stack = TRUE)
    } else {
      addStyle(wb, name, st_dec, rows = rows, cols = j, stack = TRUE)
      if (grepl("^(p|q)( |-|$)", nm)) {
        small <- which(!is.na(x) & x != 0 & abs(x) < 1e-4)
        if (length(small)) addStyle(wb, name, st_sci, rows = rows[small], cols = j, stack = TRUE)
      }
    }
  }
  if (!is.null(highlight) && any(highlight)) {
    addStyle(wb, name, st_hl, rows = rows[which(highlight)], cols = seq_along(df), gridExpand = TRUE, stack = TRUE)
  }
  setColWidths(wb, name, cols = seq_along(df), widths = "auto")
  next_row[[name]] <<- r + nrow(df) + 3
}

# --- README ---
addWorksheet(wb, "README")
readme <- c("Appendix S1. Statistical results", "",
            "Each sheet corresponds to a main or supplementary figure. Values are reported at full precision; the manuscript text rounds them.", "",
            "Phase codes",
            "M = male; I = intermediate; IP = intermediate partner; LI = late intermediate; LIP = late intermediate partner; F = female; NF = new female; NM = new male.", "",
            "Conventions",
            "Omnibus tests are Type-III tests (car::Anova). lm models report F; binomial glm models report likelihood-ratio chi-squared; mixed models (lmer/glmer) report Wald chi-squared.",
            "Pairwise p-values (M vs F, M vs I, I vs F) are uncorrected emmeans contrasts on the phase term.",
            "q-values are Benjamini-Hochberg FDR across clusters (cell abundance, enrichment) or across all genes tested in the cluster (DEGs).",
            "Cell-abundance proportions, module scores, CytoTRACE scores and gene expression are averaged per individual before summarizing; n = number of individuals.",
            "snMultiome analyses compare M, I and F only. Cluster labels are given in Appendix S2.",
            "Rows highlighted in yellow are the specific results shown in the main-text figure panel.", "", "Sheets")
sheet_index <- tibble(
  Sheet = c("Fig1", "FigS1", "Fig3A-B_S6", "Fig3C-D", "Fig3E", "Fig4", "FigS7B", "Fig5", "Fig6A-C_S10", "Fig6D-F", "Other_DEGs"),
  Contents = c("Fig. 1B-C: behavior, 11-KT, gonadal histology and volume (Table 1)",
               "Fig. S1A-F: body size, growth, gonad volume estimates",
               "Fig. 3A-B, S6A-B: cell-type abundance (all cells; neurons only)",
               "Fig. 3C-D: enrichment of DEGs and DARs by cluster",
               "Fig. 3E: sex-change classifier scores in intermediates, with QC",
               "Fig. 4A, D-F: cyp19a1b in 1_RGC; RGC subcluster abundance; neuron differentiation module",
               "Fig. S7B: CytoTRACE by 1_RGC subcluster",
               "Fig. 5B-D, G: CytoTRACE and GO module scores in 6_POA_Mixed",
               "Fig. 6A-C, S10A-E: candidate neuroendocrine populations in 6_POA_Mixed",
               "Fig. 6D-F: steroid receptor-associated expression in 6_POA_Mixed",
               "Other DEGs cited in the text (igf2bp1, 0_Spall_Mixed)"))
writeData(wb, "README", readme, startRow = 1)
addStyle(wb, "README", st_sheet, rows = 1, cols = 1)
addStyle(wb, "README", st_title, rows = which(readme %in% c("Phase codes", "Conventions", "Sheets")), cols = 1)
writeData(wb, "README", sheet_index, startRow = length(readme) + 1, colNames = FALSE)
addStyle(wb, "README", st_bold, rows = length(readme) + seq_len(nrow(sheet_index)), cols = 1)
setColWidths(wb, "README", cols = 1:2, widths = c(16, 100))

# --- Fig1 ---
new_sheet("Fig1", "Fig. 1B-C. Behavior, 11-KT, gonadal histology and gonad volume (Table 1)",
          "All eight phases. Gonad volume model includes final body length as a covariate.")
add_table("Fig1", "Omnibus tests", tests_tbl(fig1_table, "F"))
add_table("Fig1", "Means by phase", means_long(fig1_summary),
          "Gonad volume is the residual of log10 volume regressed on final body length.")

# --- FigS1 ---
new_sheet("FigS1", "Fig. S1A-F. Body size, growth and gonad volume estimates")
add_table("FigS1", "Omnibus tests", bind_rows(tests_tbl(figS1ab_table, "F"), tests_tbl(figS1cf_table, "F")),
          "Change in mass and length: experimental fish only (paired males), so M and F are not included.")
add_table("FigS1", "Means by phase", bind_rows(means_long(figS1ab_summary), means_long(figS1cf_summary)))

# --- Fig3A-B / S6 ---
new_sheet("Fig3A-B_S6", "Fig. 3A-B, S6A-B. Cell-type abundance across sex change",
          "Binomial GLM on per-individual cell counts per cluster (cbind(cluster cells, other cells) ~ Phase).")
t6a <- abund_tbl(figS6a_table, figS6a_summary)
add_table("Fig3A-B_S6", "Fig. S6A: all cells (Fig. 3A = cluster 1)", t6a, highlight = t6a$Cluster == "1")
t6b <- abund_tbl(figS6b_table, figS6b_summary)
add_table("Fig3A-B_S6", "Fig. S6B: neurons only (Fig. 3B = cluster 6)", t6b, highlight = t6b$Cluster == "6")

# --- Fig3C-D ---
new_sheet("Fig3C-D", "Fig. 3C-D. Enrichment of DEGs and DARs across clusters",
          "Chi-squared goodness-of-fit against equal expected counts per cluster; per-cluster enrichment = observed/expected, two-sided p from standardized residuals, BH-adjusted.")
for (x in list(list(d = fig3c_deg, chi = fig3c_chi, lab = "Fig. 3C: DEGs"),
               list(d = fig3d_dar, chi = fig3d_chi, lab = "Fig. 3D: DARs"))) {
  add_text("Fig3C-D", sprintf("%s: total = %d; %s(%d) = %.3f, p = %.3g", x$lab, sum(x$d$n), chi2,
                              as.integer(x$chi$parameter), x$chi$statistic, x$chi$p.value))
  e <- enrich_tbl(x$d)
  add_table("Fig3C-D", paste(x$lab, "by cluster"), e, highlight = e$`Significantly enriched`)
}

# --- Fig3E ---
new_sheet("Fig3E", "Fig. 3E. Sex-change classifier scores in intermediates",
          "Logistic regression (C = 10) trained on M vs F per cluster using DEG expression and DAR accessibility; applied to intermediates (0 = male-like, 1 = female-like).")
add_text("Fig3E", sprintf("QC: clusters retained if leave-one-out validation gave mean P(female) < 0.10 for males and > 0.90 for females (%d of %d clusters).",
                          length(qc_pass), nrow(fig3e_by_cluster)), createStyle(fontName = "Arial", fontSize = 10))
add_text("Fig3E", sprintf("Mean across QC-passing clusters (intermediates, individual x cluster): %.4f +/- %.4f SE (n = %d).",
                          fig3e_overall$mean_score, fig3e_overall$se_score, fig3e_overall$n_rows))
next_row[["Fig3E"]] <- next_row[["Fig3E"]] + 1
e3 <- fig3e_by_cluster %>%
  arrange(cluster) %>%
  transmute(Cluster = as.character(cluster), `n features (DEGs + DARs)` = n_features,
            `LOO validation: mean P(female), males` = val_m, `LOO validation: mean P(female), females` = val_f,
            `Passed QC` = qc_pass, `Intermediates: mean P(female)` = mean_score, SE = se_score, n = n)
add_table("Fig3E", "Per-cluster scores", e3, highlight = e3$`Passed QC`)

# --- Fig4 ---
new_sheet("Fig4", "Fig. 4. 1_RGC: aromatase, subcluster abundance and neuron differentiation")
add_table("Fig4", "Fig. 4A: cyp19a1b differential expression (1_RGC)", deg_tbl(filter(deg_table, Gene == "cyp19a1b")), deg_note)
add_table("Fig4", "Fig. 4A: cyp19a1b mean normalized expression", expr_tbl("cyp19a1b", "1_RGC"))
add_table("Fig4", "Fig. 4D-E: subcluster abundance within 1_RGC",
          tests_tbl(fig4de_table, paste("LR", chi2)) %>%
            mutate(Measure = dplyr::recode(Measure, `1_NSC (1_1)` = "1_NSC (subcluster 1_1)", `1_NP (1_0)` = "1_NP (subcluster 1_0)")) %>%
            dplyr::select(-Term),
          "Binomial GLM on per-individual subcluster counts relative to all 1_RGC cells.")
add_table("Fig4", "Fig. 4D-E: mean subcluster proportion",
          fig4de_summary %>% filter(cluster %in% c("1_0", "1_1")) %>%
            mutate(cluster = dplyr::recode(cluster, `1_0` = "1_NP (1_0)", `1_1` = "1_NSC (1_1)")) %>%
            wide_mif("cluster", "Subcluster"))
add_table("Fig4", "Fig. 4F: neuron differentiation module (GO:0030182) in 1_NP", dplyr::select(tests_tbl(fig4f_table, "F"), -Term),
          "AUCell module score averaged per individual.")
add_table("Fig4", "Fig. 4F: mean module score", means_long(fig4f_summary, c(score = "GO:0030182 module score")))

# --- FigS7B ---
new_sheet("FigS7B", "Fig. S7B. CytoTRACE score by 1_RGC subcluster",
          "lm on per-individual mean CytoTRACE score; 1_0 = 1_NP, 1_1 = 1_NSC.")
add_table("FigS7B", "Omnibus test", tibble(Model = figS7b_table$Model, F = figS7b_table$`F value`,
                                           `df (num)` = figS7b_table$df, `df (resid)` = figS7b_table$df_resid,
                                           `p-value (Type III)` = figS7b_table$`Type-III p-value`))
pw <- grep("^1_\\d-1_\\d$", names(figS7b_table), value = TRUE)
add_table("FigS7B", "Pairwise contrasts (uncorrected)",
          tibble(Contrast = sub("-", " vs ", pw), p = unlist(figS7b_table[1, pw])))
add_table("FigS7B", "Mean CytoTRACE by subcluster",
          means_long(figS7b_summary, group = "sub_res0.2", group_name = "Subcluster"))

# --- Fig5 ---
new_sheet("Fig5", "Fig. 5. 6_POA_Mixed: cellular immaturity and GO module scores",
          "lm on per-individual means. Higher CytoTRACE = less mature. Module scores are AUCell.")
add_table("Fig5", "Fig. 5B: CytoTRACE", dplyr::select(tests_tbl(fig5b_table, "F"), -Term))
add_table("Fig5", "Fig. 5B: mean CytoTRACE", means_long(fig5b_summary))
add_table("Fig5", "Fig. 5C-D: brain development and axon guidance receptor activity modules",
          dplyr::select(tests_tbl(fig5cd_table, "F"), -Term))
add_table("Fig5", "Fig. 5C-D: mean module scores", means_long(mutate(fig5cd_summary, trait = module)))
add_table("Fig5", "Fig. 5G: regulation of synaptic plasticity module (GO:0048167)",
          dplyr::select(tests_tbl(syn_plas_table, "F"), -Term))
add_table("Fig5", "Fig. 5G: mean module score", means_long(syn_plas_summary, c(score = "GO:0048167 module score")))

# --- Fig6A-C / S10 ---
new_sheet("Fig6A-C_S10", "Fig. 6A-C, S10A-E. Candidate gonadotroph-regulating populations in 6_POA_Mixed")
add_table("Fig6A-C_S10", "Proportion of gene+ cells",
          tests_tbl(fig6_prop_table, paste("LR", chi2)) %>% mutate(Measure = relab(Measure, gene_labels)) %>% dplyr::select(-Term),
          "Binomial GLM on per-individual counts of gene+ (normalized expression > 0) cells.")
add_table("Fig6A-C_S10", "Mean proportion of gene+ cells",
          wide_mif(mutate(fig6_prop_summary, gene = relab(gene, gene_labels)), "gene", "Gene"))
add_table("Fig6A-C_S10", "CytoTRACE of gene+ cells",
          tests_tbl(fig6_cyto_table, paste("Wald", chi2)) %>% mutate(Measure = relab(Measure, gene_labels)),
          "lmer: CytoTRACE ~ nCount_RNA + Phase + (1 | individual), fit on gene+ cells.")
add_table("Fig6A-C_S10", "Mean CytoTRACE of gene+ cells (per-individual means)",
          wide_mif(mutate(fig6_cyto_summary, gene = relab(gene, gene_labels)), "gene", "Gene"))
add_table("Fig6A-C_S10", "Differential expression (6_POA_Mixed)", deg_tbl(filter(deg_table, Cluster == "6_POA_Mixed")), deg_note)
add_table("Fig6A-C_S10", "gnrh1 expression (Fig. S10C)", dplyr::select(tests_tbl(gnrh1_table, "F"), -Term),
          "gnrh1 analysed with lm on per-individual mean expression (see Methods).")
add_table("Fig6A-C_S10", "Mean normalized expression",
          expr_tbl(c(names(deg_genes$`6`), "gnrh1"), "6_POA_Mixed"))

# --- Fig6D-F ---
new_sheet("Fig6D-F", "Fig. 6D-F. Steroid receptor-associated expression in 6_POA_Mixed",
          "lm on per-individual mean z-scored expression of genes with ESR2B, AR or PGR motifs in their promoters.")
add_table("Fig6D-F", "Omnibus tests", dplyr::select(tests_tbl(fig6df_table, "F"), -Term))
add_table("Fig6D-F", "Means by phase", means_long(fig6df_summary))

# --- Other DEGs ---
new_sheet("Other_DEGs", "Other DEGs cited in the text")
add_table("Other_DEGs", "igf2bp1 (0_Spall_Mixed)", deg_tbl(filter(deg_table, Gene == "igf2bp1")), deg_note)
add_table("Other_DEGs", "Mean normalized expression", expr_tbl("igf2bp1", "0_Spall_Mixed"))

saveWorkbook(wb, "Manuscript/Supplementary Tables/Appendix S1. Statistics.xlsx", overwrite = TRUE)
message("Wrote Manuscript/Supplementary Tables/Appendix S1. Statistics.xlsx")
