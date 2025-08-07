FeaturePlot(obj, c('lrp1bb',
                   'dock4b',
                   'abca2',
                   'srsf11',
                   'si:dkey-215k6.1'),
            reduction = 'harmony_wnn.umap')

FeaturePlot(obj, 'si:dkey-215k6.1',
            reduction = 'harmony_wnn.umap')

mecd = readRDS("Functions/mean_expression_cluster_data.rds")

mean_dat = data.frame()
for(i in c('lrp1bb',
           'dock4b',
           'abca2',
           'srsf11',
           'si:dkey-215k6.1')){
  exp = AverageExpression(obj,'RNA', i, group.by = 'res0.8_50nn_40PC_45LSI')[["RNA"]]@x
  mean_dat = rbind(mean_dat, exp)
  
}

colnames(mean_dat) = 0:26

rownames(mean_dat) = c('lrp1bb',
                       'dock4b',
                       'abca2',
                       'srsf11',
                       'si:dkey-215k6.1')

colSums(mean_dat)
rowSums(mean_dat)
