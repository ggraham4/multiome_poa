# does the GLU/GABA ratio change during sex change


obj$glut = ifelse(obj@assays$RNA$data['LOC111584103',] >0 |obj@assays$RNA$data['slc17a6b',] >0, T,  F)
obj$gaba = ifelse(obj@assays$RNA$data['gad1a',] >0 |obj@assays$RNA$data['LOC111588076',] >0, T, F)
obj$mixed = ifelse(obj$glut==T & obj$gaba == T, T,F)
obj$neuron = ifelse(obj@assays$RNA$data['elavl3',] >0 , T, NA )

neurons = subset(obj, neuron == T & mixed == F)

ratio = neurons@meta.data%>%
  group_by(individual, Status)%>%
  summarize(glut_gaba = sum(glut)/sum(gaba))

ggplot(ratio, aes(x = Status, y = glut_gaba))+
  geom_boxplot()+
  geom_point()
