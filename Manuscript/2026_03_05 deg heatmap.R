#install.packages("volcano3D")
library(pheatmap)
library(volcano3D)
degs =read.csv('/Users/ggraham/Desktop/multiome_poa/DEG Outputs/FINAL degs classified w singular.csv')

degs_6 = subset(degs, cluster ==6)

#degs_6_genes = degs_6$gene

heat_matrix = cbind(degs_6$f_m_estimate, degs_6$d_m_estimate, degs_6$d_f_estimate)
rownames(heat_matrix) = degs_6$gene
colnames(heat_matrix) = c('f_m_estimate', 'd_m_estimate', 'd_f_estimate')

# Generate letter matrix
get_letters = function(f_m_p, d_m_p, d_f_p, alpha = 0.05) {
  assigned = c(F = 1, M = 1, D = 1)
  contrasts = list(
    c("F", "M", f_m_p),
    c("D", "M", d_m_p),
    c("D", "F", d_f_p)
  )
  for (contrast in contrasts) {
    g1 = contrast[1]; g2 = contrast[2]; pval = as.numeric(contrast[3])
    if (pval < alpha && assigned[g1] == assigned[g2]) {
      assigned[g1] = max(assigned) + 1
    }
  }
  letter_map = setNames(letters[1:length(unique(assigned))],
                        as.character(sort(unique(assigned))))
  c(letter_map[as.character(assigned["F"])],
    letter_map[as.character(assigned["M"])],
    letter_map[as.character(assigned["D"])])
}

letter_matrix = t(mapply(get_letters,
                          degs_6$f_m_p.value,
                          degs_6$d_m_p.value,
                          degs_6$d_f_p.value))
rownames(letter_matrix) = degs_6$gene
colnames(letter_matrix) = c('f_m_estimate', 'd_m_estimate', 'd_f_estimate')

heat = pheatmap(heat_matrix, 
                cluster_cols = FALSE, 
                cluster_rows = FALSE,
                display_numbers = letter_matrix)
heat


