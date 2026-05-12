  go_module = function(term){
term2gene = readRDS("Function Scripts/Dependencies/Term2gene_clown_go2.rds")
term2name = readRDS('/Users/ggraham/Desktop/multiome_poa/Function Scripts/Dependencies/Term2name.rds')

go_terms = term2gene%>%
left_join(term2name, by = 'go_id')

term_genes = go_terms$aocellaris_name[ go_terms$go_id == term]
term_genes_in_obj = term_genes[term_genes %in% rownames(obj_6_only)]
print(paste0(length(term_genes_in_obj),' genes found'))

gene_pca_matrix = obj_6_only@assays$RNA$data[unique(term_genes_in_obj),]%>%t()%>%as.matrix()
gene_pca_matrix_no0 = gene_pca_matrix[,which(colSums(gene_pca_matrix)>0)]
pca = princomp(gene_pca_matrix_no0, cor = T)

  var_explained = pca$sdev^2
  max_var_pc = which.max(var_explained)
  print(paste0("PC", max_var_pc, " explains ", 
               round(var_explained[max_var_pc]/sum(var_explained)*100, 2), 
               "% of variance"))

if(mean(pca$loadings[,max_var_pc])>0){
  scores = pca$scores[,max_var_pc]
}else{
  scores=pca$scores[,max_var_pc]*-1
}
  return(scores)
  }


bd = go_module('GO:0007399') # nervous system development, 1.58 5 though...
obj_6_only$brain_dev = bd
ggplot(obj_6_only@meta.data%>%
         group_by(individual, Status)%>%
         summarize(mean_mod = mean(brain_dev)), aes(x = Status, y = mean_mod))+
  geom_boxplot()+
  geom_point()


bd = go_module('GO:0007158') # neuronal cell adhesion, ~10%
obj_6_only$brain_dev = bd
ggplot(obj_6_only@meta.data%>%
         group_by(individual, Status)%>%
         summarize(mean_mod = mean(brain_dev)), aes(x = Status, y = mean_mod))+
  geom_boxplot()+
  geom_point()

bd = go_module('GO:0048863') # stem cell differentiation 3.19
obj_6_only$brain_dev = bd
ggplot(obj_6_only@meta.data%>%
         group_by(individual, Status)%>%
         summarize(mean_mod = mean(brain_dev)), aes(x = Status, y = mean_mod))+
  geom_boxplot()+
  geom_point()
# interesting

#> in all these cases I think aucell is probably better
#> 
#> 

bd = go_module('GO:0050808') # synapse dev, 4% 
obj_6_only$brain_dev = bd
ggplot(obj_6_only@meta.data%>%
         group_by(individual, Status)%>%
         summarize(mean_mod = mean(brain_dev)), aes(x = Status, y = mean_mod))+
  geom_boxplot()+
  geom_point()



### aucell ####
library(AUCell)

  expr_matrix <- obj_6_only@assays$RNA$data
go_module_aucell <- function(term, obj) {
    set.seed(1)

  term2gene <- readRDS("Function Scripts/Dependencies/Term2gene_clown_go2.rds")
  term2name <- readRDS('/Users/ggraham/Desktop/multiome_poa/Function Scripts/Dependencies/Term2name.rds')
  
  go_terms <- term2gene %>% left_join(term2name, by = 'go_id')
  
  term_genes        <- go_terms$aocellaris_name[go_terms$go_id == term]
  term_name         <- unique(go_terms$go_name[go_terms$go_id == term])
  term_genes_in_obj <- unique(term_genes[term_genes %in% rownames(obj)])
  
  print(paste0(length(term_genes_in_obj), ' genes found for: ', term_name))
  
  if(length(term_genes_in_obj)<1){return(NULL)}
  # Build rankings from raw counts
  cell_rankings <- AUCell_buildRankings(expr_matrix, plotStats = FALSE, verbose = FALSE)
  
  # Score
gene_set <- setNames(list(term_genes_in_obj), term)
  
  auc_scores <- AUCell_calcAUC(gene_set, cell_rankings, verbose = FALSE)
  scores     <- as.numeric(getAUC(auc_scores)[1, ])
  
  return(list(scores = scores, term_name = term_name))
}

plot_go_aucell <- function(term, obj) {
    set.seed(1)

  result <- go_module_aucell(term, obj)
  
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
      x     = "Status",
      y     = "Mean AUC Score",
      title = result$term_name,
      subtitle = term
    ) +
    theme(
      axis.text.x  = element_text(angle = -45, vjust = 1, hjust = 0),
      legend.position = "none"
    )
}

model_go_aucell = function(term, obj){
  set.seed(1)
  
    result <- go_module_aucell(term, obj)
  
  obj$aucell_score <- result$scores
  
  mod_df <- obj@meta.data %>%
    group_by(individual, Status) %>%
    summarize(mean_mod = mean(aucell_score), .groups = "drop") %>%
    mutate(Status = factor(Status, levels = c('NRM', 'M', 'D', 'E', 'NF', 'F')))%>%
    subset(Status != 'NRM')
  
  regression = lm(mean_mod~Status, data = mod_df)
  return(regression)
}

model_go_aucell_mixed = function(term, obj){
  set.seed(1)
    result <- go_module_aucell(term, obj)
  
  obj$aucell_score <- result$scores
  
  mod_df <- obj@meta.data %>%
    subset(Status != 'NRM')
  
  regression = lmer(aucell_score~Status+(1|individual), data = mod_df)
  return(regression)
}

plot_go_aucell_manuscript <- function(term, obj) {
  set.seed(1)
  
  colors = c('#1965B0', '#4EB265', '#F7F056', '#DC050C', '#DC050C')
  
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
  
  # ── Score ──
  result           <- go_module_aucell(term, obj)
  obj$aucell_score <- result$scores
  
  # ── Aggregate per individual ──
  plot_df <- obj@meta.data %>%
    group_by(individual, Status) %>%
    summarize(mean_mod = mean(aucell_score), .groups = "drop") %>%
    filter(!Status %in% c('NRM')) %>%
    mutate(
      Phase = dplyr::recode(Status, 'D' = 'I', 'E' = 'LI'),
      Phase = factor(Phase, levels = c('M', 'I', 'LI','NF', 'F'))
    )
  
  # ── Stats ──
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
  
  # ── Y axis limits ──
  y_max      <- max(plot_df$mean_mod, na.rm = TRUE)
  y_lower    <- y_max * 1.10
  y_upper    <- y_max * 1.25
  y_axis_max <- y_max * 1.35
  
  plot_df_no_nf <- plot_df %>% filter(Phase != 'NF')

  # ── Plot ──
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
      subtitle = term
    ) +
    theme(legend.position = 'none') +
    # Lower tier: M vs I
    geom_signif(xmin = 1, xmax = 1.9,
                y_position = y_lower,
                annotation = format_pval(p_M_I)$label,
                color = 'black',
                tip_length = c(0, 0),
                textsize   = format_pval(p_M_I)$size,
                vjust      = format_pval(p_M_I)$adjust) +
    # Lower tier: I vs F
    geom_signif(xmin = 2.1, xmax = 5,
                y_position = y_lower,
                annotation = format_pval(p_I_F)$label,
                color = 'black',
                tip_length = c(0, 0),
                textsize   = format_pval(p_I_F)$size,
                vjust      = format_pval(p_I_F)$adjust) +
    # Upper tier: M vs F
    geom_signif(xmin = 1, xmax = 5,
                y_position = y_upper,
                annotation = format_pval(p_M_F)$label,
                color = 'black',
                tip_length = c(0, 0),
                textsize   = format_pval(p_M_F)$size,
                vjust      = format_pval(p_M_F)$adjust)
}
# Run the same  terms

# by Term
### Nervous system development
plot_go_aucell('GO:0007399', obj_6_only)  # nervous system development
model_go_aucell('GO:0007399', obj_6_only) %>%anova(test= 'Chisq')
model_go_aucell_mixed('GO:0007399', obj_6_only) %>%car::Anova(type= 3)



plot_go_aucell('GO:0007158', obj_6_only)  # neuronal cell adhesion
model_go_aucell('GO:0007158', obj_6_only) %>%anova(test= 'Chisq')
model_go_aucell_mixed('GO:0007158', obj_6_only) %>%car::Anova(type= 3)


plot_go_aucell('GO:0048863', obj_6_only)  # stem cell differentiation
model_go_aucell('GO:0048863', obj_6_only) %>%anova(test= 'Chisq') #0.07
model_go_aucell_mixed('GO:0048863', obj_6_only) %>%car::Anova(type= 3)

plot_go_aucell('GO:0050808', obj_6_only)  # synapse development
model_go_aucell('GO:0050808', obj_6_only) %>%anova(test= 'Chisq')
model_go_aucell_mixed('GO:0050808', obj_6_only) %>%car::Anova(type= 3)

plot_go_aucell('GO:0007416', obj_6_only)# synaptogenesis
model_go_aucell('GO:0007416', obj_6_only) %>%anova(test= 'Chisq')
model_go_aucell_mixed('GO:0007416', obj_6_only) %>%car::Anova(type= 3)


plot_go_aucell('GO:0022008', obj_6_only)# neurogenesis
model_go_aucell('GO:0022008', obj_6_only) %>%anova(test= 'Chisq')
model_go_aucell_mixed('GO:0022008', obj_6_only) %>%car::Anova(type= 3)


plot_go_aucell('GO:0007416', obj_6_only)# synapse assemblt
model_go_aucell('GO:0007416', obj_6_only) %>%anova(test= 'Chisq')
model_go_aucell_mixed('GO:0007416', obj_6_only) %>%car::Anova(type= 3)

plot_go_aucell('GO:0060074', obj_6_only)# synapse maturation
model_go_aucell('GO:0060074', obj_6_only) %>%anova(test= 'Chisq')
model_go_aucell_mixed('GO:0060074', obj_6_only) %>%car::Anova(type= 3)
#what!

plot_go_aucell('GO:0050767', obj_6_only)# regulation of neurogenesis
model_go_aucell_mixed('GO:0050767', obj_6_only) %>%car::Anova(type= 3)
# how are none of these significant literally just by the sheer numbers I have done 
# 18 tests ~ 1 should come back significant from type II error


plot_go_aucell('GO:0050769', obj_6_only)# positive regulation of neurogenesis
# interesting but not worth testing

plot_go_aucell('GO:0050772', obj_6_only)# positive regulation of axonogenesis
# if this isnt significant I will eat my shoe
model_go_aucell('GO:0050772', obj_6_only) %>%anova(test= 'Chisq')
#HOW HOW HOW WHAT HOW IS THIS NOT SIGNIFICANT
model_go_aucell_mixed('GO:0050772', obj_6_only) %>%car::Anova(type= 3)
plot_go_aucell('GO:0001764', obj_6_only) # neuron migration

###### Brain dev ####
plot_go_aucell('GO:0007420', obj_6_only)  # brain development this is SO significant
model_go_aucell('GO:0007420', obj_6_only) %>%anova(test= 'Chisq') #FINALLY

brain_dev = plot_go_aucell_manuscript('GO:0007420', obj_6_only) +
  labs(subtitle = 'GO:0007420 Brain Development')

    # ggsave(plot = brain_dev,
     #  file = paste0('brain_dev_go_module','.svg'),
      # device = "svg",
       #units = "in",
       #width = 2,
       #height = 2.25,
       #path = '/Users/ggraham/Desktop/multiome_poa/Manuscript/Plots/Manuscript v1.2.1/6_modules')



plot_go_aucell('GO:0007268', obj_6_only) # chemical synaptic transmission

plot_go_aucell('GO:0045664', obj_6_only) # reg of neuron diff
model_go_aucell('GO:0045664', obj_6_only) %>%anova(test= 'Chisq')


plot_go_aucell('GO:0007155', obj_6_only) # cell adhesion
model_go_aucell('GO:0007155', obj_6_only) %>%anova(test= 'Chisq')

plot_go_aucell('GO:0007162', obj_6_only) # neg cell adhesion

plot_go_aucell('GO:0045785', obj_6_only) # pos cell adhesion

plot_go_aucell('GO:0030155', obj_6_only) # reg cell adhesion

plot_go_aucell('GO:0031175', obj_6_only)  # neuron proj dev
model_go_aucell('GO:0031175', obj_6_only) %>%anova(test= 'Chisq')
# wow

plot_go_aucell('GO:0042551', obj_6_only) # neuron maturation
model_go_aucell('GO:0042551', obj_6_only) %>%anova(test= 'Chisq')






