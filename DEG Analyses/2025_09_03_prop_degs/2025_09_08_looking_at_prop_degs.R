define_degs = readRDS('Functions/define_degs2')
clown_go = readRDS('Functions/clown_go2')

prop_degs = data.frame()
for(i in 0:26){
  print(i)
  if(i == 24){next}
  data = read.csv(paste0('/Users/ggraham/Desktop/multiome_poa/DEG Outputs/2025_09_04_proportion_DEGs/cluster_',i,'.csv'))
  data$cluster =i
  prop_degs = rbind(data, prop_degs)
}

prop_degs_defined = define_degs(prop_degs)
prop_degs_defined_grouped = prop_degs_defined%>%
  group_by(cluster, short_label)%>%
  summarize(n_degs = n())

ggplot(prop_degs_defined_grouped, aes(x = as.factor(as.numeric(cluster)), y = n_degs, fill = short_label))+
  geom_bar(stat = 'identity',position = 'stack')
#yet again, cluster 2 has the most degs, and I'm inclined to believe it. Cluster 26 is concerning though, 

clown_go(prop_degs_defined$gene[prop_degs_defined$cluster==2]) %>% dotplot()
# why is this showing up in oligodendrocytes??

clown_go(prop_degs_defined$gene[prop_degs_defined$cluster==6]) %>% dotplot()
clown_go(prop_degs_defined$gene[prop_degs_defined$cluster==1]) %>% dotplot()
clown_go(prop_degs_defined$gene[prop_degs_defined$cluster==11]) %>% dotplot()
clown_go(prop_degs_defined$gene[prop_degs_defined$cluster==0]) %>% dotplot()

markers = FindMarkers(object, ident.1 = 2)

clown_go(prop_degs_defined$gene[prop_degs_defined$cluster==2 & prop_degs_defined$short_label =='Late Downregulated']) %>% dotplot()
clown_go(prop_degs_defined$gene[prop_degs_defined$cluster==2 & prop_degs_defined$short_label =='Transiently Upregulated']) %>% dotplot()


