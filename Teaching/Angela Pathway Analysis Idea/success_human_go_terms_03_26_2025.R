### steroid pathway via human go terms
library(biomaRt)
library(Seurat)
library(tidyverse)

obj <- readRDS("C:/Users/Gabe/Desktop/RNA Object.rds")
obj_human <- readRDS("C:/Users/Gabe/Desktop/RNA_object_human_names.rds")

human_to_ocellaris = read.csv("Reference/hsapiens_to_aocellaris.csv")

#make a biomart
human_mart <- useEnsembl(biomart = 'genes',
                             dataset = 'hsapiens_gene_ensembl')

att = attributes(human_mart)$attributes

human_go <- getBM(mart = human_mart,
                      attributes = c('entrezgene_accession',
                                     'go_id',
                                     'name_1006'
                      ))

joined = human_to_ocellaris%>%
  right_join(human_go, by =join_by('hsapiens_name' == 'entrezgene_accession'))%>%
  subset(!is.na(hsapiens_name)&
           hsapiens_name != ''&
           !is.na(aocellaris_name)&
           aocellaris_name != '')

steroid_genes <- joined%>%
  subset(name_1006 %in% c(
    'steroid metabolic process',
    'estrogen biosynthetic process',
    'androgen biosynthetic process',
    '5alpha-androstane-3beta,17beta-diol dehydrogenase activity',
    'estradiol 17-beta-dehydrogenase [NAD(P)] activity',
    'steroid catabolic process',
    'testosterone 17-beta-dehydrogenase (NADP+) activity',
    'aromatase activity',
    'steroid hydroxylase activity',
    'androgen catabolic process',
    'androgen metabolic process',
    'testosterone biosynthetic process'
  ))
steroid_metabolism_genes = unique(steroid_genes$aocellaris_name)

obj$endocrine <- NA
obj$endocrine <- apply(obj@assays$RNA$data[steroid_metabolism_genes,], 2, function(x) ifelse(any(x > 0), 'endocrine', 'no_endocrine'))

Idents(obj) <- 'harmony.wnn_res0.4_clusters'

DimPlot(obj, split.by = 'endocrine')
# now we're getting somewhere

Idents(obj) <- 'endocrine'

endocrine_data_cluster_level <- table(obj@meta.data$endocrine, obj@meta.data$harmony.wnn_res0.4_clusters)%>%
  as.data.frame()%>%
  pivot_wider(names_from = 'Var1', values_from = 'Freq')%>%
  mutate(prop_pos = endocrine/(endocrine+no_endocrine))

### 27 is the most +

###cursory look at sex differences in + before I do further analysis
library(lme4)
library(car)

data= data.frame()
for(i in 0:31){
  model_data = data.frame(Status = obj$Status[obj$harmony.wnn_res0.4_clusters ==i],
                          individual =obj$individual[obj$harmony.wnn_res0.4_clusters ==i],
                          endocrine = obj$endocrine[obj$harmony.wnn_res0.4_clusters ==i],
                          cluster =i)%>%
    subset(Status %in% c('M',"D","F"))%>%
    group_by(individual, Status)%>%
    summarize(endocrine_cells = sum(endocrine == 'endocrine'),
              non_endocrine_cells = sum(endocrine == 'no_endocrine'))
  
  glmer_model = glmer(cbind(model_data$endocrine_cells, model_data$non_endocrine_cells)~Status + (1|individual),
                      family = binomial(), data =model_data)
  
  av= car::Anova(glmer_model, type = 'III')%>%as.data.frame()
  
  out_data = data.frame(cluster = i,
                        p_val = av$`Pr(>Chisq)`[2],
                        singular = isSingular(glmer_model))
  data = rbind(data, out_data)  
  
}
data$signif = ifelse(data$p_val<0.05, '*', NA)
#all except 27 are significant lmao

#lets plot all clusters first
data_for_plot_1 = data.frame(Status = obj$Status,
                             individual =obj$individual,
                             endocrine = obj$endocrine)%>%
  subset(Status!= 'NRM')%>%
  group_by(individual, Status)%>%
  summarize(endocrine_cells = sum(endocrine == 'endocrine'),
            non_endocrine_cells = sum(endocrine == 'no_endocrine'))%>%
  mutate(prop_endocrine = endocrine_cells/(endocrine_cells+non_endocrine_cells))

data_for_plot_1$Status = factor(data_for_plot_1$Status, levels= c('M','D','E','NF','F'))

ggplot(data_for_plot_1, aes(x = Status, y = prop_endocrine))+
  geom_boxplot(aes(fill = Status))+
  geom_jitter(shape =1)


#4
data_for_plot_4 = data.frame(Status = obj$Status[obj$harmony.wnn_res0.4_clusters==4],
                             individual =obj$individual[obj$harmony.wnn_res0.4_clusters==4],
                             endocrine = obj$endocrine[obj$harmony.wnn_res0.4_clusters==4])%>%
  subset(Status!= 'NRM')%>%
  group_by(individual, Status)%>%
  summarize(endocrine_cells = sum(endocrine == 'endocrine'),
            non_endocrine_cells = sum(endocrine == 'no_endocrine'))%>%
  mutate(prop_endocrine = endocrine_cells/(endocrine_cells+non_endocrine_cells))

data_for_plot_4$Status = factor(data_for_plot_4$Status, levels= c('M','D','E','NF','F'))

ggplot(data_for_plot_4, aes(x = Status, y = prop_endocrine))+
  geom_boxplot(aes(fill = Status))+
  geom_jitter(shape =1)

#3
data_for_plot_3 = data.frame(Status = obj$Status[obj$harmony.wnn_res0.4_clusters==3],
                             individual =obj$individual[obj$harmony.wnn_res0.4_clusters==3],
                             endocrine = obj$endocrine[obj$harmony.wnn_res0.4_clusters==3])%>%
  subset(Status!= 'NRM')%>%
  group_by(individual, Status)%>%
  summarize(endocrine_cells = sum(endocrine == 'endocrine'),
            non_endocrine_cells = sum(endocrine == 'no_endocrine'))%>%
  mutate(prop_endocrine = endocrine_cells/(endocrine_cells+non_endocrine_cells))

data_for_plot_3$Status = factor(data_for_plot_3$Status, levels= c('M','D','E','NF','F'))

ggplot(data_for_plot_3, aes(x = Status, y = prop_endocrine))+
  geom_boxplot(aes(fill = Status))+
  geom_jitter(shape =1)





