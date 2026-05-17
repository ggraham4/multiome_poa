library(Seurat)
library(ggplot2)
library(tidyverse)
library(emmeans)
library(ggsignif)
library(AUCell)
library(lme4)

# All functions take an optional `subcluster` argument.
# If provided, the obj is subset to cells where sub_res0.2 == subcluster
# before any scoring or plotting.
obj  = readRDS("~/Desktop/optimal_clustering_rna_only.rds")

obj = FindSubCluster(obj,
                     1, 'harmony.wsnn', resolution = 0.2, subcluster.name = 'sub_res0.2')

sub_1 = subset(obj, final_clusters==1)

.subset_if_needed <- function(obj, subcluster) {
if (!is.null(subcluster)) {
obj <- subset(obj, sub_res0.2 == subcluster)
}
obj
}

go_module_aucell <- function(term, obj, subcluster = NULL) {
set.seed(1)
obj <- .subset_if_needed(obj, subcluster)

term2gene <- readRDS("Function Scripts/Dependencies/Term2gene_clown_go2.rds")
term2name <- readRDS('/Users/ggraham/Desktop/multiome_poa/Function Scripts/Dependencies/Term2name.rds')

go_terms <- term2gene %>% left_join(term2name, by = 'go_id')

term_genes        <- go_terms$aocellaris_name[go_terms$go_id == term]
term_name         <- unique(go_terms$go_name[go_terms$go_id == term])
term_genes_in_obj <- unique(term_genes[term_genes %in% rownames(obj)])

print(paste0(length(term_genes_in_obj), ' genes found for: ', term_name))

if (length(term_genes_in_obj) < 1) { return(NULL) }

expr_matrix  <- obj@assays$RNA$data
cell_rankings <- AUCell_buildRankings(expr_matrix, plotStats = FALSE, verbose = FALSE)

gene_set  <- setNames(list(term_genes_in_obj), term)
auc_scores <- AUCell_calcAUC(gene_set, cell_rankings, verbose = FALSE)
scores     <- as.numeric(getAUC(auc_scores)[1, ])

return(list(scores = scores, term_name = term_name))
}

plot_go_aucell <- function(term, obj, subcluster = NULL) {
set.seed(1)
obj <- .subset_if_needed(obj, subcluster)

result           <- go_module_aucell(term, obj)
obj$aucell_score <- result$scores

plot_df <- obj@meta.data %>%
group_by(individual, Status) %>%
summarize(mean_mod = mean(aucell_score), .groups = "drop") %>%
mutate(Status = factor(Status, levels = c('NRM', 'M', 'D', 'E', 'NF', 'F')))

ggplot(plot_df, aes(x = Status, y = mean_mod, color = Status)) +
geom_boxplot(aes(fill = Status), alpha = 0.25, outlier.shape = NA) +
geom_point(position = position_jitter(width = 0.15), size = 2) +
theme_classic(base_size = 13) +
labs(
x        = "Status",
y        = "Mean AUC Score",
title    = result$term_name,
subtitle = if (!is.null(subcluster)) paste0(term, '  |  subcluster: ', subcluster) else term
) +
theme(
axis.text.x  = element_text(angle = -45, vjust = 1, hjust = 0),
legend.position = "none"
)
}

model_go_aucell <- function(term, obj, subcluster = NULL) {
set.seed(1)
obj <- .subset_if_needed(obj, subcluster)

result           <- go_module_aucell(term, obj)
obj$aucell_score <- result$scores

mod_df <- obj@meta.data %>%
group_by(individual, Status) %>%
summarize(mean_mod = mean(aucell_score), .groups = "drop") %>%
mutate(Status = factor(Status, levels = c('NRM', 'M', 'D', 'E', 'NF', 'F'))) %>%
subset(Status != 'NRM')

lm(mean_mod ~ Status, data = mod_df)
}

model_go_aucell_mixed <- function(term, obj, subcluster = NULL) {
set.seed(1)
obj <- .subset_if_needed(obj, subcluster)

result           <- go_module_aucell(term, obj)
obj$aucell_score <- result$scores

mod_df <- obj@meta.data %>%
subset(Status != 'NRM')

lmer(aucell_score ~ Status + (1 | individual), data = mod_df)
}

plot_go_aucell_manuscript <- function(term, obj, subcluster = NULL) {
set.seed(1)
obj <- .subset_if_needed(obj, subcluster)

colors <- c('#1965B0', '#4EB265', '#F7F056', '#DC050C', '#DC050C')

format_pval <- function(p) {
if (is.na(p)) return(list(label = 'NA', size = 1, adjust = 0))
if (p >= 0.05) {
list(label = paste0('p = ', round(p, 3)), size = 2, adjust = 0)
} else if (p < 0.001) {
list(label = '***', size = 7, adjust = .6)
} else if (p < 0.01) {
list(label = '**',  size = 7, adjust = .6)
} else {
list(label = '*',   size = 7, adjust = .6)
}
}

result           <- go_module_aucell(term, obj)
obj$aucell_score <- result$scores

plot_df <- obj@meta.data %>%
group_by(individual, Status) %>%
summarize(mean_mod = mean(aucell_score), .groups = "drop") %>%
filter(!Status %in% c('NRM')) %>%
mutate(
Phase = dplyr::recode(Status, 'D' = 'I', 'E' = 'LI'),
Phase = factor(Phase, levels = c('M', 'I', 'LI', 'NF', 'F'))
)

gene_lm    <- lm(mean_mod ~ Phase, data = plot_df)
gene_pairs <- as.data.frame(pairs(emmeans(gene_lm, 'Phase'), adjust = 'none'))

get_p <- function(g1, g2) {
row <- gene_pairs %>% filter(
(contrast == paste(g1, '-', g2)) | (contrast == paste(g2, '-', g1))
)
if (nrow(row) == 0) return(NA)
row$p.value[1]
}

p_M_I <- get_p('M', 'I')
p_M_F <- get_p('M', 'F')
p_I_F <- get_p('I', 'F')

y_max      <- max(plot_df$mean_mod, na.rm = TRUE)
y_lower    <- y_max * 1.10
y_upper    <- y_max * 1.2
y_axis_max <- y_max * 1.3

plot_df_no_nf <- plot_df %>% filter(Phase != 'NF')

ggplot(plot_df, aes(x = Phase, y = mean_mod, fill = Phase)) +
geom_boxplot(data = plot_df_no_nf,
aes(x = Phase, y = mean_mod, fill = Phase),
outlier.shape = NA, width = 0.6) +
geom_point(position = position_jitter(width = 0.15), size = 1, color = 'black') +
scale_fill_manual(values = colors) +
scale_x_discrete(drop = FALSE) +
scale_y_continuous(limits = c(0.95 * min(plot_df$mean_mod), y_axis_max)) +
theme_classic() +
labs(
x        = 'Phase',
y        = 'Mean AUC Score',
title    = result$term_name,
subtitle = if (!is.null(subcluster)) paste0(term, '  |  subcluster: ', subcluster) else term
) +
theme(legend.position = 'none') +
geom_signif(xmin = 1, xmax = 1.9,
y_position = y_lower,
annotation = format_pval(p_M_I)$label,
color = 'black', tip_length = c(0, 0),
textsize = format_pval(p_M_I)$size,
vjust    = format_pval(p_M_I)$adjust) +
geom_signif(xmin = 2.1, xmax = 5,
y_position = y_lower,
annotation = format_pval(p_I_F)$label,
color = 'black', tip_length = c(0, 0),
textsize = format_pval(p_I_F)$size,
vjust    = format_pval(p_I_F)$adjust) +
geom_signif(xmin = 1, xmax = 5,
y_position = y_upper,
annotation = format_pval(p_M_F)$label,
color = 'black', tip_length = c(0, 0),
textsize = format_pval(p_M_F)$size,
vjust    = format_pval(p_M_F)$adjust)
}



# Example usage:
plot_go_aucell('GO:0030182', sub_1, subcluster = '1_0')
model_go_aucell('GO:0030182', sub_1, subcluster = '1_0') %>% anova(test = 'Chisq') 
# oh its significant, but perhaps not in the direction I expected

neuron_differentiation_1_NP=  plot_go_aucell_manuscript('GO:0030182', sub_1, subcluster = '1_0')

#ggsave(plot = neuron_differentiation_1_NP,
#      file = "neuron_differentiation_1_NP.svg",
#    device = "svg",
# #    units = "in",
#  width = 2,
# height = 2,
#path = '~/Desktop/multiome_poa/Manuscript/Plots/Manuscript v1.2.1/RGC_modules')


plot_go_aucell('GO:0060853', sub_1, subcluster = '1_1')
# cell fate commitment, expect change in 1_1 where it is elevated in intermediates
# not gonna be significant 
model_go_aucell('GO:0060853', sub_1, subcluster = '1_1') %>% anova(test = 'Chisq') 


plot_go_aucell('GO:0032502', sub_1, subcluster = '1_1')
# dev process
model_go_aucell('GO:0032502', sub_1, subcluster = '1_1') %>% anova(test = 'Chisq') 


plot_go_aucell('GO:0022008', sub_1, subcluster = '1_1')# neurogenesis

model_go_aucell('GO:0022008', sub_1, subcluster = '1_1') %>% anova(test = 'Chisq') 


plot_go_aucell('GO:0019827', sub_1, subcluster = '1_1')# stem cell maintenence, 
# could be a result, show that this ISNT changing

plot_go_aucell('GO:0048843', sub_1, subcluster = '1_0') # neg reg of axon guidance


#plot_go_aucell('GO:0050681', sub_1, subcluster = '1_1')







