sub1 = subset(obj, final_clusters==1)
sub1$arob = sub1@assays$RNA$data['LOC111577263',]

sub1@meta.data%>%
  group_by(Status)%>%
  summarize(mean_arob = mean(arob),
            se_arob = sd(arob)/sqrt(n()))

sub1 = FindSubCluster(sub1,
                     1, 'harmony.wsnn', resolution = 0.2, subcluster.name = 'sub_res0.2')

cells_ind = sub1@meta.data%>%
  group_by(individual)%>%
  summarize(n_cells = n())

cells_sub_ind = sub1@meta.data%>%
  group_by(individual, Status, sub_res0.2)%>%
  summarize(n_cells_in = n())%>%
    subset(Status%in%c('M','D','F'))

cells_total = cells_ind%>%
  right_join(cells_sub_ind, by = 'individual')
cells_total$prop = (cells_total$n_cells_in/cells_total$n_cells)*100

props =cells_total%>%
  group_by(Status, sub_res0.2)%>%
  summarize(mean_prop = mean(prop),
            se_prop = sd(prop)/sqrt(n()))

DimPlot(sub1, group.by = 'sub_res0.2')

sub_1 = sub1

term2gene <- readRDS("Function Scripts/Dependencies/Term2gene_clown_go2.rds")
term_genes <- term2gene$aocellaris_name[term2gene$go_id == "GO:0030182"]
term_genes_in_obj <- unique(term_genes[term_genes %in% rownames(sub_1)])

expr_matrix <- sub_1@assays$RNA$data
cell_rankings <- AUCell_buildRankings(
  expr_matrix,
  plotStats = FALSE,
  verbose = FALSE
)

gene_set <- setNames(list(term_genes_in_obj), "GO:0030182")
auc_scores <- AUCell_calcAUC(
  gene_set,
  cell_rankings,
  verbose = FALSE
)
sub_1$GO_0030182_score <- as.numeric(getAUC(auc_scores)[1, ])

GO_0030182_summary <- sub_1@meta.data %>%
  group_by(individual, Status) %>%
  summarise(
    mean_score = mean(GO_0030182_score, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(Status %in% c("M", "D", "F")) %>%
  group_by(Status) %>%
  summarise(
    mean = mean(mean_score, na.rm = TRUE),
    se = sd(mean_score, na.rm = TRUE) / sqrt(sum(!is.na(mean_score))),
    .groups = "drop"  
    )%>% mutate(
    mean = format(mean, scientific = TRUE),
    se = format(se, scientific = TRUE)
  )
GO_0030182_summary
