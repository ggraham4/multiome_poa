##########################################################################################
#  OMNIBUS MANUSCRIPT PLOTS
#  Every data panel of the manuscript in the order it is first cited in the text.
#
#  HOW TO USE
#   * Run from the project root (the folder that contains Manuscript/, Measures/, ...).
#   * Section 1 holds the plotting functions. Every bar plot, box plot and physiology plot
#     below is drawn by one of them, so they all share the same widths, colors, jitter,
#     fonts and significance brackets. Edit the look in ONE place (section 1).
#   * Section 2 loads the data. These CSVs are written by
#     "2026_10_08 manuscript_statistics_appendix_S1.R", so plots and statistics use
#     identical numbers. Brackets show the pairwise p-values from that script
#     (emmeans, no p-adjustment: M vs I, I vs F, M vs F).
#   * Section 3 is the figures. Each panel is a few lines.
#   * Section 4 saves every panel as an svg, one folder per figure, in "Manuscript/consistent plots/"
#     (e.g. "Fig 1/Fig1C_testicular_pct.svg"). Run export_panels() to write them.
#
#  PANELS THAT ARE NOT DRAWN HERE (UMAPs, heatmaps, GO dotplots, networks, drawings)
#   They are listed in order in section 3 with the script that made them.
#
#  CONTENTS (text order)
#   Fig 1A | S1A,B | 1B | 1C | S1C-F | 1D | S2 | 2A-D | S3-S5 | S6A | 3A | S6B | 3B | 3C | 3D |
#   3E | 4A | 4B | S7A | 4C | S7B | S8 | 4D | 4E | 4F | 5A | 5B | 5C | 5D | S9 | 5E,F | 5G |
#   6A-C | S10A-E | 6D-F | 6G-I | S11 | 7
##########################################################################################


# ========================================================================================
# 1. LOOK AND PLOTTING FUNCTIONS
# ========================================================================================
library(tidyverse)
library(emmeans)
library(multcomp)     # compact-letter display (cld) for the physiology plots
library(ggsignif)     # significance brackets

# ---- constants (same as the existing plotting scripts) --------------------------------
phase_levels <- c("M", "I", "LI", "NF", "F")
phase_colors <- c(M = "#1965B0", I = "#4EB265", LI = "#F7F056", NF = "#7BAFDE", F = "#DC050C")
class_colors <- c("Early Downregulated" = "red",         "Early Upregulated" = "#006400",
                  "Late Downregulated" = "blue",         "Late Upregulated" = "#000000",
                  "Progressively Downregulated" = "purple", "Progressively Upregulated" = "gray",
                  "Transiently Downregulated" = "brown", "Transiently Upregulated" = "orange")
out_dir <- "Manuscript/consistent plots"

# saves one panel as svg into a folder named after its figure
save_panel <- function(plot, figure, file, width = 2, height = 2) {
  folder <- file.path(out_dir, figure)
  dir.create(folder, showWarnings = FALSE, recursive = TRUE)
  ggsave(plot = plot, file = paste0(file, ".svg"), device = "svg", units = "in",
         width = width, height = height, path = folder)
}

# ---- significance labels --------------------------------------------------------------
# stars when p < 0.05, otherwise the p-value in small text
p_value_annotation <- function(p) {
  if (p < 0.001)      list(label = "***", size = 7, adjust = 0.6)
  else if (p < 0.01)  list(label = "**",  size = 7, adjust = 0.6)
  else if (p < 0.05)  list(label = "*",   size = 7, adjust = 0.6)
  else                list(label = paste0("p = ", sprintf("%.3f", p)), size = 2, adjust = 0)
}

# ---- one vertical layout for every plot ----------------------------------------------------
# All positions are fractions of the axis height (0 = bottom of the axis, 1 = top), so two panels placed side by
# side line up: bracket tiers and letters sit at the same height, and the empty space above and below the data is the
# same. Change a number here and every plot follows.
layout <- list(
  bottom_gap = 0.06,   # empty space under the lowest point (box plots, dot plots). Bar plots start at 0 instead
  data_top   = 0.60,   # highest data value (bar + SE, or highest point): bar plots and box plots
  phys_top   = 0.82,   # highest data value in the physiology plots (no brackets, so the data can fill more of the panel)
  lower_tier = 0.72,   # M-I and I-F brackets
  upper_tier = 0.84,   # M-F bracket
  letters    = 0.93    # compact letters in the physiology plots
)

# y-axis limits and a helper at(fraction) -> y value, from the lowest and highest data values
axis_layout <- function(dmin, dmax, from_zero = FALSE, data_top = layout$data_top) {
  if (from_zero) {
    ylo  <- 0
    span <- dmax / data_top
  } else {
    span <- (dmax - dmin) / (data_top - layout$bottom_gap)
    ylo  <- dmin - layout$bottom_gap * span
  }
  list(limits = c(ylo, ylo + span), at = function(fraction) ylo + fraction * span)
}

# adds the three brackets: M-I and I-F on the lower tier, M-F on the upper tier
add_brackets <- function(plot, pv, ax) {
  one <- function(xmin, xmax, y, p) {
    a <- p_value_annotation(p)
    geom_signif(xmin = xmin, xmax = xmax, y_position = y, annotation = a$label, color = "black",
                tip_length = c(0, 0), textsize = a$size, vjust = a$adjust)
  }
  plot + one(1, 1.9, ax$at(layout$lower_tier), pv$m_i) +
         one(2.1, 5, ax$at(layout$lower_tier), pv$i_f) +
         one(1, 5, ax$at(layout$upper_tier), pv$m_f)
}

# ---- bar plot: mean +/- SE of one value per individual (proportions) -------------------
# df needs columns Phase and the value column. pv = list(m_i, i_f, m_f). Bars start at 0, so only the top has empty space.
bar_plot <- function(df, value, ylab, pv, title = NULL, percent = TRUE) {
  df$value <- df[[value]]
  df$Phase <- factor(df$Phase, levels = phase_levels)
  summ <- df %>% group_by(Phase) %>%
    summarise(mean = mean(value, na.rm = TRUE), se = sd(value, na.rm = TRUE) / sqrt(n()), .groups = "drop")
  ax <- axis_layout(0, max(df$value, summ$mean + summ$se, na.rm = TRUE), from_zero = TRUE)
  p <- ggplot(df, aes(x = Phase, y = value)) +
    geom_bar(data = summ, aes(x = Phase, y = mean, fill = Phase), stat = "identity", width = 0.6, inherit.aes = FALSE) +
    geom_errorbar(data = summ, aes(x = Phase, ymin = mean - se, ymax = mean + se), width = 0.4, linewidth = 0.5, inherit.aes = FALSE) +
    geom_point(size = 1, position = position_jitter(width = 0.1, seed = 42)) +
    scale_x_discrete(drop = FALSE) +
    scale_fill_manual(values = phase_colors) +
    theme_classic() + theme(legend.position = "none") +
    labs(x = "Phase", y = ylab, title = title)
  p <- p + if (percent) scale_y_continuous(labels = scales::percent, limits = ax$limits, expand = c(0, 0))
           else         scale_y_continuous(limits = ax$limits, expand = c(0, 0))
  add_brackets(p, pv, ax)
}

# ---- box plot: one value per individual (expression, CytoTRACE, module scores) ---------
# NF is a single fish, so it is drawn as a point without a box.
box_plot <- function(df, value, ylab, pv, title = NULL) {
  df$value <- df[[value]]
  df$Phase <- factor(df$Phase, levels = phase_levels)
  ax <- axis_layout(min(df$value, na.rm = TRUE), max(df$value, na.rm = TRUE))
  p <- ggplot(df, aes(x = Phase, y = value)) +
    geom_boxplot(data = subset(df, Phase != "NF"), aes(x = Phase, y = value, fill = Phase), outlier.shape = NA, width = 0.6) +
    geom_point(size = 1, position = position_jitter(width = 0.1, seed = 42)) +
    scale_x_discrete(drop = FALSE) +
    scale_y_continuous(limits = ax$limits, expand = c(0, 0)) +
    scale_fill_manual(values = phase_colors) +
    theme_classic() + theme(legend.position = "none") +
    labs(x = "Phase", y = ylab, title = title)
  add_brackets(p, pv, ax)
}

# ---- physiology plots (Fig 1, S1): open circles, mean crossbar, SE bars, letters ------------
# Letters = compact letter display of all pairwise comparisons (emmeans, no adjustment).
phase_letters <- function(model) {
  letters_df <- as.data.frame(cld(emmeans(model, "Phase"), Letters = letters, adjust = "none", alpha = 0.05))
  letters_df <- letters_df[!is.na(letters_df$emmean), ]          # drop phases that are not in this model
  letters_df$Phase <- as.character(letters_df$Phase)             # a factor here would add empty columns to the x axis
  letters_df
}

physiology_theme <- theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_rect(linewidth = 1, color = "black"),
        legend.key = element_blank(), legend.position = "none")

# Fig 1B and S1A,B: each phase split into dominant / subordinate partner, arrows give the fold difference
dominance_plot <- function(measures, y, ylab, model) {
  measures$y <- measures[[y]]
  means <- measures %>% group_by(Phase, Dominance) %>%
    summarize(mean = mean(y, na.rm = TRUE), se = sd(y, na.rm = TRUE) / sqrt(sum(!is.na(y))), .groups = "drop")
  ax <- axis_layout(min(measures$y, means$mean - means$se, na.rm = TRUE),
                    max(measures$y, means$mean + means$se, na.rm = TRUE), data_top = layout$phys_top)
  m <- function(ph) mean(measures$y[measures$Phase == ph], na.rm = TRUE)
  pairs_x <- list(c("F", "M", 1.5), c("I", "IP", 3.5), c("LI", "LIP", 5.5), c("NF", "NM", 7.5))
  p <- ggplot(measures, aes(y = y, x = Phase, color = Dominance, shape = Dominance)) +
    stat_summary(data = subset(measures, !Phase %in% c("NM", "NF")), geom = "crossbar", fun = "mean") +
    geom_errorbar(data = means, aes(x = Phase, y = mean, ymin = mean - se, ymax = mean + se),
                  color = "black", width = 0.4, linewidth = 0.5) +
    geom_vline(xintercept = c(2.5, 4.5, 6.5)) +
    geom_point(size = 2, position = position_jitterdodge(0.25)) +
    geom_text(data = phase_letters(model), aes(x = Phase, y = ax$at(layout$letters), label = .group),
              size = 3, inherit.aes = FALSE) +
    scale_y_continuous(limits = ax$limits, expand = c(0, 0)) +
    labs(x = "Phase", y = ylab) + physiology_theme
  for (pr in pairs_x) {                       # fold-difference arrow between the two partners
    a <- m(pr[1]); b <- m(pr[2]); x <- as.numeric(pr[3])
    p <- p + annotate("segment", x = x, xend = x, y = a, yend = b, arrow = arrow(ends = "both", length = unit(0.2, "cm"))) +
      annotate("text", x = x + 0.5, y = mean(c(a, b)), label = sprintf("%.2fx", max(a, b) / min(a, b)), size = 2.5)
  }
  p
}

# Fig 1C and S1C-F: one open circle per fish, mean crossbar, SE bars, letters
phase_dot_plot <- function(df, y, ylab, model, x_order = NULL, percent = FALSE, minimal = FALSE) {
  df$y <- df[[y]]
  df <- df[!is.na(df$y), ]                 # phases with no measurement (e.g. no testis in F) get no x-axis slot
  df$Phase <- as.character(df$Phase)
  means <- df %>% group_by(Phase) %>%
    summarize(mean = mean(y, na.rm = TRUE), se = sd(y, na.rm = TRUE) / sqrt(sum(!is.na(y))), .groups = "drop")
  ax <- axis_layout(min(df$y, means$mean - means$se, na.rm = TRUE),
                    max(df$y, means$mean + means$se, na.rm = TRUE), data_top = layout$phys_top)
  p <- ggplot(df, aes(y = y, x = Phase)) +
    stat_summary(geom = "crossbar", fun = "mean") +
    geom_errorbar(data = means, aes(x = Phase, y = mean, ymin = mean - se, ymax = mean + se), width = 0.4, linewidth = 0.5) +
    geom_point(size = 2, position = position_jitterdodge(0.25), shape = 1) +
    geom_text(data = phase_letters(model), aes(x = Phase, y = ax$at(layout$letters), label = .group), size = 3, inherit.aes = FALSE) +
    labs(x = "Phase", y = ylab)
  p <- p + if (minimal) theme_minimal() else physiology_theme
  p <- p + if (percent) scale_y_continuous(labels = scales::percent, limits = ax$limits, expand = c(0, 0))
           else         scale_y_continuous(limits = ax$limits, expand = c(0, 0))
  if (!is.null(x_order)) p <- p + scale_x_discrete(limits = x_order)
  p
}

# classification stacked bars (Fig 3C,D): number of DEGs / DARs per cluster
classification_plot <- function(df, ylab) {
  counts <- df %>% group_by(cluster, classification) %>% summarize(n = n(), .groups = "drop")
  counts$cluster <- factor(as.numeric(as.character(counts$cluster)), levels = c(-1, 0:26))   # -1 = blank spacer
  ggplot(counts, aes(x = cluster, y = n, fill = classification)) +
    geom_bar(stat = "identity") + theme_minimal() +
    scale_fill_manual(values = class_colors) +
    scale_x_discrete(drop = FALSE, labels = c("", 0:26)) +
    labs(x = "Cluster", y = ylab, fill = "Classification") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
}


# ========================================================================================
# 2. DATA
# ========================================================================================
src <- "Manuscript/Supplementary Tables/Appendix S1 source data"
read_src <- function(name) read.csv(file.path(src, paste0(name, ".csv")), check.names = FALSE)

# p-values from the statistics script. source = table name, trait = measure name in Appendix S1
all_pvalues <- read_src("pairwise_pvalues_for_plots")
get_pvals <- function(source, trait) {
  d <- all_pvalues[all_pvalues$source == source & all_pvalues$Trait == trait, ]
  pick <- function(a, b) d$p[d$contrast %in% c(paste0(a, "-", b), paste0(b, "-", a))]
  list(m_i = pick("M", "I"), i_f = pick("I", "F"), m_f = pick("M", "F"),
       omnibus = unique(d$omnibus_p), q = unique(d$q))
}

# physiology table (all eight phases)
measures <- read_src("fig1_S1_measures_all_phases")
measures$Phase <- factor(measures$Phase, levels = c("F", "M", "I", "IP", "LI", "LIP", "NF", "NM"))
measures$Dominance <- ifelse(measures$Phase %in% c("F", "I", "LI", "NF"), "Dominant", "Subordinate")

# per-individual tables (M, I, LI, NF, F; LI and NF are shown but never tested)
counts_cells   <- read_src("fig3A_S6A_cell_counts_per_individual_cluster")
counts_neurons <- read_src("fig3B_S6B_neuron_counts_per_individual_cluster")
counts_rgc     <- read_src("fig4DE_rgc_subcluster_counts_per_individual")
expression     <- read_src("named_DEGs_mean_expression_per_individual")
gene_pos       <- read_src("fig6ABC_S10_gene_positive_counts_per_individual")
cyto_gene_pos  <- read_src("fig6ABC_S10_cytotrace_gene_positive_per_individual")

cluster_names <- readxl::read_excel("Manuscript/cluster_names_and_regions.xlsx")
cluster_label <- setNames(cluster_names$`Cluster Name`, as.character(cluster_names$`Cluster `))


# ========================================================================================
# 3. FIGURES  (text order)
# ========================================================================================

# ---- Fig 1A: schematic (drawn in PowerPoint; the sex-state curves come from 2026_04_29 updates v1.2 fig 1.R) ----

# ---- Fig S1A,B: final body mass and length (dominant / subordinate pairs) -------------------------------
fig_S1A <- dominance_plot(measures, "mass_final_cm",   "Mass (g)",    lm(mass_final_cm ~ Phase,   data = measures))
fig_S1B <- dominance_plot(measures, "length_final_cm", "Length (cm)", lm(length_final_cm ~ Phase, data = measures))

# ---- Fig 1B: parental behaviors, time in nest, 11-KT -----------------------------------------------------
fig_1B_behaviors <- dominance_plot(measures, "Behaviors_Day_2", "Parental Behaviors", lm(Behaviors_Day_2 ~ Phase, data = measures))
fig_1B_time      <- dominance_plot(measures, "Time_Day_2",      "Time in nest (s)",   lm(Time_Day_2 ~ Phase,      data = measures))
fig_1B_kt        <- dominance_plot(measures, "Log_11KT",        "Log10 11KT (pg/ml)", lm(Log_11KT ~ Phase,        data = measures))

# ---- Fig 1C: % testicular, % ovarian, length-corrected gonad volume (females and NF not shown) -----------
x_gonad <- c("M", "I", "IP", "LI", "LIP", "NM")
gonad   <- subset(measures, !Phase %in% c("F", "NF"))
fig_1C_testicular <- phase_dot_plot(gonad, "Percent_Testicular", "Testicular %", lm(Percent_Testicular ~ Phase, data = measures), x_gonad, percent = TRUE)
fig_1C_ovarian    <- phase_dot_plot(gonad, "Percent_Ovarian",    "Ovarian %",    lm(Percent_Ovarian ~ Phase,    data = measures), x_gonad, percent = TRUE)

volume_model <- lm(Log10_Volume ~ Phase + length_final_cm, data = measures)          # model used for the test and letters
gonad_vol    <- subset(measures, !Phase %in% c("NF", "F") & !is.na(Log10_Volume))
gonad_vol$corrected_volume <- residuals(lm(Log10_Volume ~ length_final_cm, data = gonad_vol))   # volume corrected for body length
fig_1C_volume <- phase_dot_plot(gonad_vol, "corrected_volume", "Corrected Volume", volume_model, x_gonad)

# ---- Fig S1C-F: change in mass and length (experimental fish only), gonad volume estimates -------------
experiment <- subset(measures, Condition == "Experiment")
experiment$Change_Mass_prop   <- as.numeric(experiment$Change_Mass) / 100
experiment$Change_Length_prop <- as.numeric(experiment$Change_Length) / 100
fig_S1C <- phase_dot_plot(experiment, "Change_Mass_prop",   "Change Mass %",   lm(Change_Mass ~ Phase,   data = experiment), percent = TRUE, minimal = TRUE)
fig_S1D <- phase_dot_plot(experiment, "Change_Length_prop", "Change Length %", lm(Change_Length ~ Phase, data = experiment), percent = TRUE, minimal = TRUE)
fig_S1E <- phase_dot_plot(measures, "Log10_Testicular_Estimate", "Log10 Testicular Volume (px)", lm(Log10_Testicular_Estimate ~ Phase, data = measures), minimal = TRUE)
fig_S1F <- phase_dot_plot(measures, "Log10_Ovarian_Estimate",    "Log10 Ovarian Volume (px)",    lm(Log10_Ovarian_Estimate ~ Phase,    data = measures), minimal = TRUE)

# ---- Fig 1D (H&E images), S2 (dissection drawings) -- images, not plotted here ----------------------------

# ---- Figs 2A-D, S3, S4, S5: UMAP, marker dotplot, Visium HD / Cell2location maps, label transfer ----------
#   2A  UMAP:                 DimPlot of the multiome object (see "2026_04_29 sox2 prop and clust 6.R" for the UMAP code)
#   2B  marker dotplot:       "2026_06_22 dotplot.R"
#   2C,D, S3, S4 heatmaps and section maps: "2026_05_19 mapping multiome to visium.r", "2026_06_21 c2l projection on vis sections.R"
#   S5  label transfer:       "2026_08_23 cross projecting parker and multiome.R"

# ---- Fig S6A: cell-type abundance, all clusters with a significant phase effect (q < 0.05) -------------
# One small bar plot per cluster, titled "Cluster <number>" (1.5 x 1.5 in)
proportion_panel <- function(counts, cluster, pv, ylab) {
  d <- counts[counts$cluster == cluster, ]
  d$value <- d$prop
  bar_plot(d, "value", ylab, pv, title = paste0("Cluster ", cluster))
}
significant_clusters <- function(source) {
  q <- unique(all_pvalues[all_pvalues$source == source, c("Trait", "q")])
  q$Trait[!is.na(q$q) & q$q < 0.05]
}
fig_S6A <- lapply(significant_clusters("figS6a"), function(i)
  proportion_panel(counts_cells, i, get_pvals("figS6a", i), "% of Cells"))
names(fig_S6A) <- significant_clusters("figS6a")

# ---- Fig 3A: 1_RGC share of all cells -------------------------------------------------------------------
fig_3A <- proportion_panel(counts_cells, "1", get_pvals("figS6a", "1"), "% of Cells") + labs(title = "1_RGC")

# ---- Fig S6B: neuron clusters with a significant phase effect -----------------------------------------------
fig_S6B <- lapply(significant_clusters("figS6b"), function(i)
  proportion_panel(counts_neurons, i, get_pvals("figS6b", i), "% of Neurons"))
names(fig_S6B) <- significant_clusters("figS6b")

# ---- Fig 3B: 6_POA_Mixed share of neurons ------------------------------------------------------------------
fig_3B <- proportion_panel(counts_neurons, "6", get_pvals("figS6b", "6"), "% of Neurons") + labs(title = "6_POA_Mixed")

# ---- Fig 3C, 3D: number of DEGs and DARs per cluster, colored by temporal classification ----------------
degs <- read.csv("DEG Outputs/FINAL degs classified w singular.csv", fileEncoding = "latin1")
dars <- read.csv("Collaboration/all_clusters_DARs_peak_level_classified_with_support.csv", fileEncoding = "latin1")
fig_3C <- classification_plot(transmute(degs, cluster = cluster, classification = short_label), "DEGs")
fig_3D <- classification_plot(transmute(dars, cluster = cluster_id, classification = classification), "DARs")

# ---- Fig 3E: classifier prediction in intermediates (clusters that pass the leave-one-out QC) -------------
scores <- read_src("fig3E_classifier_scores_intermediates")
validation <- read.csv("Sex Classifier and Linearity/progress_classifier/Including DARs/logistic_mf_validation_07_13_2026.csv") %>%
  group_by(cluster, status) %>% summarize(mean_prob = mean(proba), .groups = "drop") %>%
  pivot_wider(names_from = status, values_from = mean_prob, names_prefix = "mean_prob_")
confident_clusters <- validation$cluster[validation$mean_prob_m < 0.10 & validation$mean_prob_f > 0.90]
classifier <- scores %>% filter(cluster %in% confident_clusters) %>%
  group_by(cluster, neurotransmitter) %>%
  summarize(mean_prediction = mean(prediction), se_prediction = sd(prediction) / sqrt(n()), .groups = "drop")
fig_3E <- ggplot(classifier, aes(x = fct_reorder(as.factor(cluster), mean_prediction, .desc = TRUE), y = mean_prediction)) +
  geom_errorbar(aes(ymin = mean_prediction - se_prediction, ymax = mean_prediction + se_prediction), width = 0.4) +
  geom_point(aes(color = neurotransmitter, shape = neurotransmitter), size = 2) +
  scale_color_manual(values = c("purple", "black", "#5bb450", "maroon")) +
  scale_shape_manual(values = c(15, 16, 17, 18)) +
  labs(x = "Cluster", y = "Prediction") + theme_minimal() + theme(legend.position = "none")

# ---- Fig 4A: brain aromatase (cyp19a1b) in 1_RGC. p-values from the negative-binomial DEG model ---------
cyp19a1b <- subset(expression, Cluster == "1_RGC" & Gene == "cyp19a1b")
fig_4A <- box_plot(cyp19a1b, "expr", "Normalized Expression", get_pvals("named_DEGs", "1_RGC | cyp19a1b"), title = "cyp19a1b | 1_RGC")

# ---- Fig 4B (DEG heatmap), S7A (DAR heatmap), 4C (subcluster UMAP) -- not drawn here ----------------------
#   4B   "2026_05_10 rgc deg matrix.R"          S7A  "2026_07_19 dar matrix1.R"
#   4C   DimPlot of sub_res0.2 (see "2026_09_20 rgc inline stats.R")

# ---- Fig S7B: CytoTRACE by 1_RGC subcluster (letters = all pairwise, no adjustment) -------------------------
cyto_sub <- read_src("figS7B_cytotrace_by_rgc_subcluster_per_individual")
cyto_sub <- subset(cyto_sub, Status != "NRM")
cyto_sub$sub <- factor(cyto_sub[["sub_res0.2"]])
cyto_sub_letters <- as.data.frame(cld(emmeans(lm(mean_cyto ~ sub, data = cyto_sub), "sub"), Letters = letters, adjust = "none", alpha = 0.05))
fig_S7B <- ggplot(cyto_sub, aes(x = sub, y = mean_cyto)) +
  geom_boxplot(outlier.shape = NA) + geom_point(shape = 1, size = 2) +
  geom_text(data = cyto_sub_letters, aes(x = sub, y = 1, label = .group), size = 3, inherit.aes = FALSE) +
  labs(x = "Subcluster", y = "Mean CytoTRACE") + theme_minimal() + ylim(0, 1)

# ---- Fig S7C, S8A-D: GO dotplots and subcluster marker dotplot -- not drawn here --------------------------
#   GO of 1_NSC / 1_NP / 1_gfap / 1_shha: "2026_05_10 go modules rgcs.R"        marker dotplot: "2026_07_17 rgc subcluster marker dotplot.R"

# ---- Fig 4D,E: 1_NSC (subcluster 1_1) and 1_NP (subcluster 1_0) as a share of 1_RGC ---------------------------
fig_4D <- proportion_panel(counts_rgc, "1_1", get_pvals("fig4de", "1_NSC (1_1)"), "% of 1_RGC") + labs(title = "1_NSC")
fig_4E <- proportion_panel(counts_rgc, "1_0", get_pvals("fig4de", "1_NP (1_0)"),  "% of 1_RGC") + labs(title = "1_NP")

# ---- Fig 4F: neuron differentiation module (GO:0030182) in 1_NP ----------------------------------------------
neuron_diff <- read_src("fig4F_neuron_differentiation_module_per_individual")
fig_4F <- box_plot(neuron_diff, "score", "AUCell score", get_pvals("fig4f", "neuron differentiation (GO:0030182), 1_NP"), title = "1_NP")

# ---- Fig 5A: UMAP with 6_POA_Mixed highlighted -- "2026_04_29 sox2 prop and clust 6.R" ------------------------

# ---- Fig 5B: CytoTRACE of 6_POA_Mixed -----------------------------------------------------------------------
cyto_6 <- read_src("fig5B_cytotrace_6_POA_Mixed_per_individual")
fig_5B <- box_plot(cyto_6, "cyto", "CytoTRACE", get_pvals("fig5b", "CytoTRACE, 6_POA_Mixed"), title = "6_POA_Mixed")

# ---- Fig 5C, 5D: brain development and axon guidance receptor activity modules -----------------------------
brain_dev <- read_src("fig5C_brain_development_module_per_individual")
axon_guid <- read_src("fig5D_axon_guidance_module_per_individual")
fig_5C <- box_plot(brain_dev, "score", "AUCell score", get_pvals("fig5cd", "brain development (GO:0007420)"), title = "brain development")
fig_5D <- box_plot(axon_guid, "score", "AUCell score", get_pvals("fig5cd", "axon guidance receptor activity (GO:0008046)"), title = "axon guidance")

# ---- Fig S9 (GO of 6_POA_Mixed DEGs): "2026_05_11 cluster 6 go.R" ; Fig 5E,F (DEG / DAR heatmaps): "2026_05_03 6 npo degs matrix.R", "2026_07_14 dar matrix 6.R" ----

# ---- Fig 5G: regulation of synaptic plasticity module (GO:0048167) ----------------------------------------------
synaptic <- read_src("fig5G_synaptic_plasticity_module_per_individual")
fig_5G <- box_plot(synaptic, "score", "AUCell score", get_pvals("fig5g", "regulation of synaptic plasticity (GO:0048167)"), title = "synaptic plasticity")

# ---- Fig 6A-C and S10A-E: one row per gene: expression | % gene+ cells | CytoTRACE of gene+ cells ----------
gene_row <- function(gene, has_expression = TRUE) {
  row <- list()
  if (has_expression) {
    e <- subset(expression, Cluster == "6_POA_Mixed" & Gene == gene)
    # gnrh1 expression is tested with an lm on per-individual means (the negative-binomial model did not converge)
    pv <- if (gene == "gnrh1") get_pvals("gnrh1_expression", "gnrh1 expression, 6_POA_Mixed")
          else                 get_pvals("named_DEGs", paste("6_POA_Mixed |", gene))
    row$expression <- box_plot(e, "expr", "Expression", pv)       # no title: keeps every Fig 6 / S10 panel the same size
  }
  pos <- gene_pos[gene_pos$gene == gene, ]
  row$percent <- bar_plot(pos, "prop", paste0("Percent ", gene, "+"), get_pvals("fig6_prop", gene))
  if (gene %in% cyto_gene_pos$gene) {
    cy <- cyto_gene_pos[cyto_gene_pos$gene == gene, ]
    row$cytotrace <- box_plot(cy, "cyto", "CytoTRACE", get_pvals("fig6_cyto", gene))
  }
  row
}
fig_6A <- gene_row("drd3")
fig_6B <- gene_row("tacr3a")
fig_6C <- gene_row("cckb")
fig_S10A <- gene_row("npy7r")
fig_S10B <- gene_row("nmbr")
fig_S10C <- gene_row("gnrh1")
fig_S10D <- gene_row("pgr")
fig_S10E <- gene_row("ar_like")

# ---- Fig 6D-F: steroid receptor-associated expression, 6G-I: promoter accessibility ------------------------------
steroid_expression    <- read_src("fig6DF_steroid_expression_scores_per_individual")
steroid_accessibility <- read_src("fig6GI_steroid_accessibility_scores_per_individual")
steroid_expression$Phase    <- steroid_expression$group
steroid_accessibility$Phase <- steroid_accessibility$group
steroid_panel <- function(df, column, trait, source, receptor)
  box_plot(df, column, paste(receptor, "Specific Gene Z-Score"), get_pvals(source, trait))
fig_6D <- steroid_panel(steroid_expression, "ESR2B_score", "ESR2B-associated expression", "fig6df", "ESR2B")
fig_6E <- steroid_panel(steroid_expression, "AR_score",    "AR-associated expression",    "fig6df", "AR")
fig_6F <- steroid_panel(steroid_expression, "PGR_score",   "PGR-associated expression",   "fig6df", "PGR")
fig_6G <- steroid_panel(steroid_accessibility, "ESR2B_accessibility_score", "ESR2B-associated accessibility", "fig6gi", "ESR2B")
fig_6H <- steroid_panel(steroid_accessibility, "AR_accessibility_score",    "AR-associated accessibility",    "fig6gi", "AR")
fig_6I <- steroid_panel(steroid_accessibility, "PGR_accessibility_score",   "PGR-associated accessibility",   "fig6gi", "PGR")

# ---- Fig S11A,B: regulatory networks; Fig 7: model schematic -- drawn outside R ---------------------------------


# ========================================================================================
# 4. EXPORT: every panel as svg, one folder per figure (the DEG / DAR stacked bars 3C,D are not exported)
#    Sizes follow the earlier plotting scripts: 2 x 2 in for bar and box plots, 3 x 2 in for the dominance
#    plots, 2.5 x 2 in for Fig 1C, 1.5 x 1.5 in for the Fig S6 cluster panels, 2.5 x 1.5 in for the classifier.
# ========================================================================================
export_panels <- function() {
  one <- function(figure, file, plot, width = 2, height = 2) save_panel(plot, figure, file, width, height)

  one("Fig S1", "FigS1A_body_mass",           fig_S1A, 3, 2)
  one("Fig S1", "FigS1B_body_length",         fig_S1B, 3, 2)
  one("Fig S1", "FigS1C_change_mass",         fig_S1C, 3, 2)
  one("Fig S1", "FigS1D_change_length",       fig_S1D, 3, 2)
  one("Fig S1", "FigS1E_testis_volume",       fig_S1E, 3, 2)
  one("Fig S1", "FigS1F_ovary_volume",        fig_S1F, 3, 2)

  one("Fig 1", "Fig1B_parental_behaviors",    fig_1B_behaviors, 3, 2)
  one("Fig 1", "Fig1B_time_in_nest",          fig_1B_time,      3, 2)
  one("Fig 1", "Fig1B_11KT",                  fig_1B_kt,        3, 2)
  one("Fig 1", "Fig1C_testicular_pct",        fig_1C_testicular, 2.5, 2)
  one("Fig 1", "Fig1C_ovarian_pct",           fig_1C_ovarian,    2.5, 2)
  one("Fig 1", "Fig1C_corrected_gonad_volume", fig_1C_volume,    2.5, 2)

  for (i in names(fig_S6A)) one("Fig S6", paste0("FigS6A_cluster_", i, "_pct_of_cells"),   fig_S6A[[i]], 1.5, 1.5)
  for (i in names(fig_S6B)) one("Fig S6", paste0("FigS6B_cluster_", i, "_pct_of_neurons"), fig_S6B[[i]], 1.5, 1.5)

  one("Fig 3", "Fig3A_cluster_1_RGC_pct_of_cells",     fig_3A)
  one("Fig 3", "Fig3B_cluster_6_POA_pct_of_neurons",   fig_3B)
  one("Fig 3", "Fig3E_classifier_prediction",          fig_3E, 2.5, 1.5)

  one("Fig 4", "Fig4A_cyp19a1b_expression_1_RGC",      fig_4A)
  one("Fig 4", "Fig4D_1_NSC_pct_of_1_RGC",             fig_4D)
  one("Fig 4", "Fig4E_1_NP_pct_of_1_RGC",              fig_4E)
  one("Fig 4", "Fig4F_neuron_differentiation_module_1_NP", fig_4F)
  one("Fig S7", "FigS7B_cytotrace_by_RGC_subcluster",  fig_S7B)

  one("Fig 5", "Fig5B_cluster_6_cytotrace",            fig_5B)
  one("Fig 5", "Fig5C_brain_development_module",       fig_5C)
  one("Fig 5", "Fig5D_axon_guidance_module",           fig_5D)
  one("Fig 5", "Fig5G_synaptic_plasticity_module",     fig_5G)

  gene_rows <- list("Fig 6" = list(`6A_drd3` = fig_6A, `6B_tacr3a` = fig_6B, `6C_cckb` = fig_6C),
                    "Fig S10" = list(S10A_npy7r = fig_S10A, S10B_nmbr = fig_S10B, S10C_gnrh1 = fig_S10C,
                                     S10D_pgr = fig_S10D, S10E_ar_like = fig_S10E))
  panel_names <- c(expression = "expression", percent = "pct_of_cluster_6", cytotrace = "cytotrace_gene_positive_cells")
  for (figure in names(gene_rows)) for (panel in names(gene_rows[[figure]])) for (k in names(gene_rows[[figure]][[panel]]))
    one(figure, paste0("Fig", panel, "_", panel_names[[k]]), gene_rows[[figure]][[panel]][[k]])

  one("Fig 6", "Fig6D_ESR2B_expression_score",         fig_6D)
  one("Fig 6", "Fig6E_AR_expression_score",            fig_6E)
  one("Fig 6", "Fig6F_PGR_expression_score",           fig_6F)
  one("Fig 6", "Fig6G_ESR2B_accessibility_score",      fig_6G)
  one("Fig 6", "Fig6H_AR_accessibility_score",         fig_6H)
  one("Fig 6", "Fig6I_PGR_accessibility_score",        fig_6I)
}
# export_panels()     # remove the leading # to write all svgs
