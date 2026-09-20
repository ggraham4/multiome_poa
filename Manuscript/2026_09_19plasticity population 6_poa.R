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



