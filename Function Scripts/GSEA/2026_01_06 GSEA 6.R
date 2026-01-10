# Analysis of DEGs enriched in dominants

deg_data = read.csv("DEG Outputs/05_11_2025 Neg Bin w Doms New_clustering/cluster_6.csv")

deg_data_sorted_d_m = deg_data[['d_m_estimate']]
names(deg_data_sorted_d_m) = deg_data[['gene']]

deg_data_sorted_d_m <- sort(deg_data_sorted_d_m, decreasing = TRUE)

clown_gsea = readRDS('Functions/clown_gsea.rds')
t = clown_gsea(deg_data_sorted_d_m)
t%>%dotplot()

deg_data_sorted_d_f = deg_data[['d_f_estimate']]
names(deg_data_sorted_d_f) = deg_data[['gene']]

deg_data_sorted_d_f <- sort(deg_data_sorted_d_f, decreasing = TRUE)

t2 = clown_gsea(deg_data_sorted_d_f)
t2%>%dotplot()

deg_data_sorted_f_m = deg_data[['f_m_estimate']]
names(deg_data_sorted_f_m) = deg_data[['gene']]

deg_data_sorted_f_m <- sort(deg_data_sorted_f_m, decreasing = TRUE)

t3 = clown_gsea(deg_data_sorted_f_m)
t3%>%dotplot()
