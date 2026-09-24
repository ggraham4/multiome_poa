library(Seurat)
library(ggplot2)
library(tidyverse)
library(ggsignif)

obj  = readRDS("~/Desktop/optimal_clustering_rna_only.rds")
colors = c('#1965B0', '#4EB265', '#F7F056' , '#DC050C')



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

sub_6 = subset(obj, final_clusters == 6)

sub_6$gene_pos = colSums(
  sub_6@assays$RNA$data[degs_plasticity, ]
) > 0

plot_plas_prop <- sub_6@meta.data %>%
  group_by(individual, Status) %>%
  filter(Status != "NRM") %>%
  summarize(plas_score = mean(gene_pos), .groups = "drop") %>%
  mutate(
    Phase = case_when(
      Status == "D" ~ "I",
      Status == "E" ~ "LI",
      TRUE ~ Status
    ),
    Phase = factor(Phase, levels = c("M", "I", "LI", "NF", "F"))
  )

subset_data <- filter(plot_plas_prop, Status != "NF")  


plasticity_plot = ggplot(plot_plas_prop, aes(x = Phase, y = plas_score)) +
  geom_boxplot(
    data = subset_data,
    aes(x = Phase, y = plas_score, fill = Phase)
  ) +
      geom_point(size = 1,
               position = position_jitter(width = 0.1, seed = 42)) +

  scale_x_discrete(
    limits = c("M", "I", "LI", "NF", "F"),
    drop = FALSE
  ) +
  theme_classic() +
  scale_fill_manual(values = colors) +
  theme(legend.position = "none") +
  scale_y_continuous(labels = scales::percent, limits = c(min(plot_plas_prop$plas_score), 
                                                          1)) +
  labs(y = "% of 6_POA_Mixed")
plasticity_plot
  
ggsave(plot = plasticity_plot,
         file = 'plasticity_plot.svg',
         device = "svg",
         units = "in",
         width = 2,
         height = 2,
         path = paste0("Manuscript/Plots"))

# logistic
plas_logistic = sub_6@meta.data%>%
  group_by(individual, Status)%>%
  subset(Status %in% c("M",'D','F'))%>%
  summarize(n_cells = n())%>%
  right_join(sub_6@meta.data%>%
               group_by(individual)%>%
               summarize(n_plas = sum(gene_pos)), 'individual')%>%
  mutate(successes = n_plas,
         failures = n_cells - n_plas)

model = glm(cbind(plas_logistic$successes, plas_logistic$failures)~Status, 
              data = plas_logistic,
              family = 'binomial')  

anova(model, test = 'Chisq')
# not significant but I still think report it cause it makes the most sense to me

plas_summary <- sub_6@meta.data %>%
  filter(Status %in% c("M", "D", "F")) %>%
  group_by(individual, Status) %>%
  summarize(
    n_cells = n(),
    n_plas = sum(gene_pos),
    .groups = "drop"
  ) %>%
  mutate(
    prop = n_plas / n_cells
  ) %>%
  group_by(Status) %>%
  summarize(
    mean_prop = mean(prop*100, na.rm = TRUE),
    se_prop = sd(prop*100, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

plas_summary


gene_prop_summary_func <- function(gene){
  tem <- sub_6@meta.data
  tem$gene_pos <- sub_6@assays$RNA$data[gene, ] > 0
  
  tem %>%
    filter(Status %in% c("M", "D", "F")) %>%
    group_by(individual, Status) %>%
    summarize(
      n_cells = n(),
      prop = mean(gene_pos) * 100,
      .groups = "drop"
    ) %>%
    group_by(Status) %>%
    summarize(
      mean_prop = mean(prop, na.rm = TRUE),
      se_prop = sd(prop, na.rm = TRUE) / sqrt(n()),
      .groups = "drop"
    )
}

gene_prop_summary_func("drd3")

gene_prop_summary_func("tacr3a")

gene_prop_summary_func("npy7r")

gene_prop_summary_func("nmbr")

gene_prop_summary_func("LOC111571064")

gene_prop_summary_func("cckb")


#cyto = CytoTRACE(sub_6@assays$RNA$data %>% as.matrix())
sub_6$cyto = cyto$CytoTRACE

gene_pos_cyto_sum = function(gene){
  tem = sub_6@meta.data
  tem$gene = sub_6@assays$RNA$data[gene,]>0
  
  tem= subset(tem, gene >0)
  
  tem%>%
        filter(Status %in% c("M", "D", "F")) %>%
    group_by(Status)%>%
    summarize(mean_cyt =mean(cyto),
              se_cyt = sd(cyto)/sqrt(n()))
  
}

gene_pos_cyto_sum('cckb')

gene_pos_cyto_sum('tacr3a')

gene_pos_cyto_sum('drd3')


gex_summary = function(gene){
    tem = sub_6@meta.data
  tem$gene = sub_6@assays$RNA$data[gene,]
  
  tem%>%
        filter(Status %in% c("M", "D", "F")) %>%
    group_by(Status)%>%
    summarize(mean_gene =mean(gene),
              se_gene = sd(gene)/sqrt(n()))

}
gex_summary('drd3')
gex_summary('tacr3a')
gex_summary('npy7r')
gex_summary('nmbr')

gex_summary('LOC111571064')# gnrh

gex_summary('cckb')


### go module expression
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

tem = sub_6@meta.data
tem$scores = scores

return(  tem%>%
        filter(Status %in% c("M", "D", "F")) %>%
    group_by(Status)%>%
    summarize(mean_score =mean(scores),
              se_score = sd(scores)/sqrt(n()))%>%
        mutate(mean_score = format(mean_score, scientific = TRUE),
               se_score = format(se_score, scientific = TRUE))

)
  }
  
  go_module_aucell('GO:0007420', sub_6)
  go_module_aucell('GO:0097484', sub_6)



