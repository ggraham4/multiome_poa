#make data for cyp19a1b TF search

counts <- rgc@assays$RNA$data%>%as.matrix()
counts_no_0 <- counts[which(rowSums(counts)!=0),which(colSums(counts)!=0)]
counts_no_0_t = t(counts_no_0)

write.csv(counts_no_0_t, '/Users/ggraham/Desktop/RGC ML/counts.csv')
write.csv(colnames(counts_no_0_t), '/Users/ggraham/Desktop/RGC ML/counts_columns.csv')
write.csv(rownames(counts_no_0_t), '/Users/ggraham/Desktop/RGC ML/counts_rows.csv')
