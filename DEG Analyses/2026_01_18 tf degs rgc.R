deg_tfs = c('skida1',
            'gli1',
            'LOC111563338',# coup
            'LOC111568198', # pbnox2
            'LOC111587547', # dbx
            'sox6',
            'foxp2',
            'nkx2.2a',
            'sox1a'
            )

deg_tfs_name = c('skida1',
            'gli1',
            'COUPTFII',# coup
            'pbnox2', # 
            'dbx1bl', # 
            'sox6',
            'foxp2',
            'nkx2.2a',
            'sox1a'
            )


mecp  = readRDS('Functions/mean_expression_cluster_plot.rds')          
p = mecp(obj, 'znf710b',1)
for(i in 1:length(deg_tfs)){
  print(i)
  d = mecp(obj, deg_tfs[i], 1 )+labs(title = deg_tfs_name[i])
  p = p+d
  
}
p

mecp(obj, 'sox1a',1)


TFs = c('skida1',
            'gli1',
            'LOC111563338',# coup
            'LOC111568198', # pbnox2
            'LOC111587547', # dbx
            'sox6',
            'foxp2',
            'nkx2.2a'
            )

mat <- obj@assays$RNA$data[TFs, obj@meta.data$final_clusters==1] > 0
bin <- t(mat%>%as.matrix()) * 1

inter <- t(bin) %*% bin
uni <- outer(colSums(bin), colSums(bin), "+") - inter

jaccard_mat <- inter / uni

comp = prcomp(jaccard_mat)
library(factoextra)
fviz_pca_ind(comp)


conditional_probability = function(b, a){
  matrix = obj@assays$RNA$data[c(a, b), obj@meta.data$final_clusters==1]%>%as.matrix()%>%t()
  matrix_binary = matrix>0
  
  probability_a = sum(matrix_binary[,a])/nrow(matrix_binary)
  probability_b = sum(matrix_binary[,b])/nrow(matrix_binary)

  probability_b_given_a = sum(matrix_binary[,b][matrix_binary[,a]==T])/sum(matrix_binary[,a]==T)

  probability_a_given_b = (probability_b_given_a*probability_a)/probability_b
  return(probability_a_given_b)
}

conditional_probability('nkx2.2a','foxp2')
conditional_probability('sox6','gli1')

new_mat = matrix(NA, length(TFs), length(TFs))
rownames(new_mat)= TFs
colnames(new_mat)= TFs
for(a in TFs){
  for(b in TFs){
    out = conditional_probability(a,b)
    new_mat[a,b]=out
  }
}

prcomp(new_mat)%>%fviz_pca_ind()
# whelp what do we do here

