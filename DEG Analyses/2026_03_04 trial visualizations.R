degs = read.csv('~/Desktop/multiome_poa/DEG Outputs/FINAL degs classified w singular.csv')

degs$sum_abs_estimates = abs(degs$f_m_estimate)+abs(degs$d_m_estimate)+abs(degs$d_f_estimate)
degs$gene_name=ifelse(is.na(degs$gene_name), degs$gene, degs$gene_name)

ggplot(subset(degs, cluster == 6), aes(y = fct_reorder(gene_name,sum_abs_estimates) , x  = sum_abs_estimates, 
                                       color = av_q.value))+
  geom_point()+
  scale_color_gradientn(colors = c('red','yellow' ,'blue'))+
  labs(x = 'Absolute Value Sum of Estimates', title = '6_POA_2 DEGs', color = 'fdr', y ='')+
  theme_minimal()

library(ggrepel)
ggplot(subset(degs, cluster == 6), aes(y = f_m_estimate , x  = d_f_estimate, 
                                       color = d_m_estimate))+
  geom_point()+
    scale_color_gradientn(colors = c('darkred','orange' ,'blue'))+
  geom_text_repel(aes(label= gene))+
  labs(x = 'D-F Estimate', y = 'F-M Estimate', color = 'D-M Estimate')+
  theme_classic()

degs_6 = subset(degs, cluster ==6)
pca_mat = cbind(degs_6$f_m_estimate, degs_6$d_f_estimate, degs_6$d_f_estimate)

rownames(pca_mat) = degs_6$gene
p = princomp(pca_mat)
library(factoextra)
fviz_pca_ind(p, geom = 'text', habillage = degs_6$short_label)
