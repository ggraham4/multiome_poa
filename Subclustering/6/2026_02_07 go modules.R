# generic module function where I input a GO term and it figures it out

go_module = function(term){
term2gene = readRDS("Function Scripts/Dependencies/Term2gene_clown_go2.rds")
term2name = readRDS('/Users/ggraham/Desktop/multiome_poa/Function Scripts/Dependencies/Term2name.rds')

go_terms = term2gene%>%
left_join(term2name, by = 'go_id')

term_genes = go_terms$aocellaris_name[ go_terms$go_id == term]
term_genes_in_obj = term_genes[term_genes %in% rownames(sub_6)]
print(paste0(length(term_genes_in_obj),' genes found'))

gene_pca_matrix = sub_6@assays$RNA$data[unique(term_genes_in_obj),]%>%t()%>%as.matrix()
gene_pca_matrix_no0 = gene_pca_matrix[,which(colSums(gene_pca_matrix)>0)]
pca = princomp(gene_pca_matrix_no0, cor = T)

  var_explained = pca$sdev^2
  max_var_pc = which.max(var_explained)
  print(paste0("PC", max_var_pc, " explains ", 
               round(var_explained[max_var_pc]/sum(var_explained)*100, 2), 
               "% of variance"))

if(mean(pca$loadings[,max_var_pc])>0){
  scores = pca$scores[,max_var_pc]
}else{
  scores=pca$scores[,max_var_pc]*-1
}

dupe = sub_6

dupe$module = scores

pca_ind =dupe@meta.data%>%
  group_by(Status, individual, sub.cluster)%>%
  summarize(mean_module = mean(module),
            se_module = sd(module)/sqrt(n()))

assign(paste0('module_', term), pca_ind, envir = .GlobalEnv)
p = ggplot(pca_ind, aes(x = Status, y= mean_module))+
  geom_boxplot()+
  geom_point()+
  facet_wrap(~sub.cluster, scales ='free')
return(p)
}

go_module('GO:0031012')
go_module('GO:0072534')

go_module('GO:0007416')# synapse assembly
#model

synapse_assembly = aov(mean_module~Status, data = subset(`module_GO:0007416`,Status!= 'NRM' & sub.cluster == '6_2' ))
summary(synapse_assembly) # truly cant believe the pp value is that high?
synapse_assembly_1 = aov(mean_module~Status, data = subset(`module_GO:0007416`,Status!= 'NRM' & sub.cluster == '6_1' ))
summary(synapse_assembly_1) # like how is that p value lower?


go_module('GO:0034056') #ERE binding
# 6_0 wooh

ere_0 = aov(mean_module~Status, data = subset(`module_GO:0034056`,Status!= 'NRM' & sub.cluster == '6_0' ))
summary(ere_0) 

go_module('GO:0001508') #action potential, some potentially interesting patterns here

go_module('GO:0022008')# neurogenesis, I dont think anything here

go_module('GO:0016081') # synaptic vesicle docking

go_module('GO:0007218') #neuropeptide signaling

go_module('GO:0030521') # androgen receptor
#woah a lot of shit looks to be going on here but also what'
# 6_3?

ar_3 = aov(mean_module~Status, data = subset(`module_GO:0030521`,Status!= 'NRM' & sub.cluster == '6_3' ))
summary(ar_3) #nope

go_module('GO:0007420') # brain development 6_1 maybe

bd1 = aov(mean_module~Status, data = subset(`module_GO:0007420`,Status!= 'NRM' & sub.cluster == '6_3' ))
summary(bd1) #how is that 0.6 what 

go_module('GO:0010975')# neuron projection development

go_module('GO:0042551') #neuron maturation

nm3 = aov(mean_module~Status, data = subset(`module_GO:0042551`,Status!= 'NRM' & sub.cluster == '6_3' ))
summary(nm3)
# maybe male gnRH neurons are immature? cyto disagrees BUT proportios say new 6_3 neurons are born so...
# i dont know the z-range is so low this is probably nothing




