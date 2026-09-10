# Chisq test for degs and dars
library(tidyverse)

degs =read.csv('/Users/ggraham/Desktop/multiome_poa/DEG Outputs/FINAL degs classified w singular.csv')
dars =read.csv("Collaboration/all_clusters_DARs_peak_level_classified_with_support.csv")

degs_clust = degs%>%
  group_by(cluster)%>%
  summarize(n = n())

dars_clust = dars%>%
  group_by(cluster_id)%>%
  summarize(n = n())

# chisq test
dar_chisq = chisq.test(dars_clust$n)
deg_chisq = chisq.test(degs_clust$n)


degs_clust <- degs_clust %>%
  mutate(
    expected = deg_chisq$expected,
    enrichment = n / expected,
    residual = (n - expected) / sqrt(expected),
    p_value = 2 * pnorm(-abs(residual)),
    p_adj = p.adjust(p_value, method = "BH")
  )%>%
  mutate(signif = p_adj < 0.05 & enrichment > 1)

dars_clust <- dars_clust %>%
  mutate(
    expected = dar_chisq$expected,
    enrichment = n / expected,
    residual = (n - expected) / sqrt(expected),
    p_value = 2 * pnorm(-abs(residual)),
    p_adj = p.adjust(p_value, method = "BH")
  )%>%
  mutate(signif = p_adj < 0.05 & enrichment > 1)
