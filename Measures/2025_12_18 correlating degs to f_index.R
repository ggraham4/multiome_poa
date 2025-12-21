{
  library(ggsignif)
  library(patchwork)
  library(ggplot2)
  library(Seurat)
  library(dplyr)
  library(tidyverse)
  library(CytoTRACE)
  library(BiocManager)
  library(lme4)
  library(enrichR)
  library(WGCNA)
  `%notin%` = Negate(`%in%`)
  library(car)
  library(emmeans)
  library(glmGamPoi)
  library(scran)
  library(parallel)
    library(parallel)
  library(reshape2)
library(ggplot2)
  library(Seurat)
  library(tidyr)
  library(lme4)
  library(dplyr)
  library(MASS)
  library(Signac)
  library('glmGamPoi')
  library(scran)
  library(emmeans)
  library(openxlsx)
  library(ggplot2)
  library(stringr)
  library(forcats)
  library(clusterProfiler)
  library(biomaRt)
  library(hdWGCNA)
  library(factoextra)
  
  mass_formula=readRDS( 'Functions/Theory/f_index_mass.rds')
kt_formula=readRDS( 'Functions/Theory/f_index_11kt.rds')
kt_e2_formula=readRDS( 'Functions/Theory/f_index_kt_e2.rds')
e2_formula=readRDS( 'Functions/Theory/f_index_e2.rds')
test_formula=readRDS( 'Functions/Theory/f_index_testicular.rds')
beh_formula=readRDS( 'Functions/Theory/f_index_behavior.rds')
time_formula=readRDS( 'Functions/Theory/f_index_time.rds')
mass_formula =readRDS('Functions/Theory/f_index_mass.rds')

scale_2 <- function(x) {
  x <- as.numeric(x)
  rng <- range(x, na.rm = TRUE)
  (x - rng[1]) / (rng[2] - rng[1])
}
  
}
obj = readRDS("~/Desktop/optimal_clustering_rna_only.rds")


measures = read.csv('Measures/all_data.csv')
measures$Time_Day_2= as.numeric(measures$Time_Day_2)
measures$Behaviors_Day_2= as.numeric(measures$Behaviors_Day_2)
measures$Log_11KT = as.numeric(measures$Log_11KT)

status_to_phase = list('D'='I',
                       'E' = 'LI',
                       'EP' = 'LIP',
                       'S' = 'IP',
                       'M'='M',
                       'F'='F',
                       'NF'='NF',
                       'NM'='NM')
measures$Phase = unlist(status_to_phase[measures$Status])
measures$Phase = factor(measures$Phase, levels =c('F','M','I','IP','LI','LIP','NF','NM'))
measures$Dominance = ifelse(measures$Phase %in%c('F','I','LI','NF'), 'Dominant', 'Subordinate')
measures$Behaviors_Day_2 = as.numeric(measures$Behaviors_Day_2)

# f index
vars4= c(#'Percent_Testicular',
  'Log_11KT',
  'mass_final_cm'
)

##kt + e2
pc_data = measures[complete.cases(measures[,vars4]),]%>%
  subset(!Phase %in% c('S',
                       'EP',
                       'IDK',
                       'NM'))
X4 <- scale(pc_data[, vars4])
pca4 <- prcomp(X4, center = TRUE, scale. = TRUE)
fviz_pca_biplot(pca4, habillage = pc_data$Phase)
# similar proportion of variation

pc_data$f_index <- pca4$x[,1] *-1

pc_data$individual = pc_data$Fish
obj@meta.data$individual

obj@meta.data$individual%in%pc_data$individual%>%table()

# named vector: names = individual, values = f_index
f_index_vec <- pc_data$f_index
names(f_index_vec) <- pc_data$individual

# map individual → cell
obj@meta.data$f_index <- f_index_vec[obj@meta.data$individual]

pc_data$Phase = factor(pc_data$Phase, levels =c('IP','LIP','NM','M','I','LI','NF','F'))

ggplot(pc_data, aes(x = Phase, y = f_index))+
  geom_boxplot()+
    geom_point()

### Any patterns? ###
obj@meta.data$f_index = as.numeric(as.character(obj@meta.data$f_index))
hist(obj@meta.data$f_index)

FeaturePlot(obj, features = 'f_index', reduction = 'harmony_wnn.umap')
FeaturePlot(obj,
            features = 'f_index',
            reduction = 'harmony_wnn.umap', 
            cells=rownames(obj@meta.data[obj@meta.data$final_clusters==6,]))

FeaturePlot(obj,
            features = 'f_index',
            reduction = 'harmony_wnn.umap', 
            cells=rownames(obj@meta.data[obj@meta.data$final_clusters==1,]))

FeaturePlot(obj,
            features = 'f_index',
            reduction = 'harmony_wnn.umap', 
            cells=rownames(obj@meta.data[obj@meta.data$final_clusters==1,]))


FeaturePlot(obj,
            features = 'f_index',
            reduction = 'harmony_wnn.umap', 
            cells=rownames(obj@meta.data[obj@meta.data$final_clusters==0,]))

### aromatase expression ##
obj@meta.data$aromatase = obj@assays$RNA$data['LOC111577263',]

arom_group = obj@meta.data%>%
  group_by(individual, Status, final_clusters)%>%
  subset(final_clusters == 1)%>%
  summarize(mean_arom = mean(aromatase),
            f_index =mean(f_index))
ggplot(arom_group, aes(x = f_index, y = mean_arom))+
  geom_point(aes(color = Status))+
  geom_smooth(method = 'lm')

obj@meta.data$esr2b = obj@assays$RNA$data['esr2b',]
esr2b_group = obj@meta.data%>%
  group_by(individual, Status, final_clusters)%>%
  subset(final_clusters == 6)%>%
  summarize(mean_esr2b= mean(esr2b),
            f_index =mean(f_index))
ggplot(esr2b_group, aes(x = f_index, y = mean_esr2b))+
  geom_point(aes(color = Status))

obj@meta.data$tacr3a = obj@assays$RNA$data['tacr3a',]
tacr3a_group = obj@meta.data%>%
  group_by(individual, Status, final_clusters)%>%
  subset(final_clusters == 6)%>%
  summarize(mean_tacr3a= mean(tacr3a),
            f_index =mean(f_index))

ggplot(tacr3a_group, aes(x = f_index, y = mean_tacr3a))+
  geom_point(aes(color = Status))

obj@meta.data$drd3 = obj@assays$RNA$data['drd3',]
drd3_group = obj@meta.data%>%
  group_by(individual, Status, final_clusters)%>%
  subset(final_clusters == 6)%>%
  summarize(mean_drd3= mean(drd3),
            f_index =mean(f_index))

ggplot(drd3_group, aes(x = f_index, y = mean_drd3))+
  geom_point(aes(color = Status))

#read in DEG data
degs = read.csv("DEG Outputs/FINAL_degs_classified_08_11_2025.csv")

cor.out = data.frame()
for(i in 1:nrow(degs)){
  print(i)
  temp = degs[i,]
  gene = temp$gene
  cluster = temp$cluster
  
  cells = rownames(obj@meta.data[obj@meta.data$final_clusters==cluster,])
  
  obj@meta.data$temp_gene= NA
  obj@meta.data$temp_gene[obj@meta.data$final_clusters==cluster] = obj@assays$RNA$data[paste0(gene),cells]
  
  mean_exp = obj@meta.data%>%
  group_by(individual, Status, final_clusters)%>%
  subset(final_clusters == cluster & !is.na(f_index))%>%
  summarize(mean_gene= mean(temp_gene),
            f_index =mean(f_index))

  correlation = cor(mean_exp$mean_gene, mean_exp$f_index)
  
  newd = data.frame(gene = gene,
                    cluster = cluster,
                    correlation = correlation)
  
  cor.out= rbind(cor.out, newd)
  
}

obj@meta.data$rragd = obj@assays$RNA$data['rragd',]
rragd_group = obj@meta.data%>%
  group_by(individual, Status, final_clusters)%>%
  subset(final_clusters == 0)%>%
  summarize(mean_rragd= mean(rragd),
            f_index =mean(f_index))

ggplot(rragd_group, aes(x = f_index, y = mean_rragd))+
  geom_point(aes(color = Status))

obj@meta.data$igf2bp1 = obj@assays$RNA$data['igf2bp1',]
igf2bp1_group = obj@meta.data%>%
  group_by(individual, Status, final_clusters)%>%
  subset(final_clusters == 0)%>%
  summarize(mean_igf2bp1= mean(igf2bp1),
            f_index =mean(f_index))

ggplot(igf2bp1_group, aes(x = f_index, y = mean_igf2bp1))+
  geom_point(aes(color = Status))
#perhaps not surprising that an insulin like growth factor binding protein
#is negatively correlated with growth and testosterone

obj@meta.data$esr2b = obj@assays$RNA$data['esr2b',]
esr2b_group = obj@meta.data%>%
  group_by(individual, Status, final_clusters)%>%
  subset(final_clusters == 9)%>%
  summarize(mean_esr2b= mean(esr2b),
            f_index =mean(f_index))

ggplot(esr2b_group, aes(x = f_index, y = mean_esr2b))+
  geom_point(aes(color = Status))


obj@meta.data$coup = obj@assays$RNA$data['LOC111563338',]
coup_group = obj@meta.data%>%
  group_by(individual, Status, final_clusters)%>%
  subset(final_clusters == 1)%>%
  summarize(mean_coup= mean(coup),
            f_index =mean(f_index))

ggplot(coup_group, aes(x = f_index, y = mean_coup))+
  geom_point(aes(color = Status))

obj@meta.data$grm8d = obj@assays$RNA$data['grm8b',]
grm8d_group = obj@meta.data%>%
  group_by(individual, Status, final_clusters)%>%
  subset(final_clusters == 19)%>%
  summarize(mean_grm8d= mean(grm8d),
            f_index =mean(f_index))

ggplot(grm8d_group, aes(x = f_index, y = mean_grm8d))+
  geom_point(aes(color = Status))


### proportions
n_cells=obj@meta.data%>%
  group_by(individual, Status)%>%
  summarize(n_cells = n(),
            f_index = mean(f_index))

n_cells_cluster=obj@meta.data%>%
  group_by(individual, final_clusters)%>%
  summarize(n_cells_cluster = n())

tog = n_cells_cluster%>%
  right_join(n_cells, by = 'individual')

tog$prop = tog$n_cells_cluster/tog$n_cells

ggplot(subset(tog, final_clusters ==1), aes(x = f_index, y = prop))+
  geom_point(aes(color = Status))+
  labs(title = 'Proportion of Radial Glia')

ggplot(subset(tog, final_clusters ==6), aes(x = f_index, y = prop))+
  geom_point(aes(color = Status))+
  labs(title = 'Proportion of 6')

ggplot(subset(tog, final_clusters ==0), aes(x = f_index, y = prop))+
  geom_point(aes(color = Status))+
  labs(title = 'Proportion of 0')

ggplot(subset(tog, final_clusters ==9), aes(x = f_index, y = prop))+
  geom_point(aes(color = Status))+
  labs(title = 'Proportion of 9')

ggplot(subset(tog, final_clusters ==22), aes(x = f_index, y = prop))+
  geom_point(aes(color = Status))+
  labs(title = 'Proportion of 22')

mod = gam(data = subset(tog, final_clusters ==22), f_index~prop)
summary(mod)
plot(mod)

ggplot(subset(tog, final_clusters ==24), aes(x = f_index, y = prop))+
  geom_point(aes(color = Status))+
  labs(title = 'Proportion of 24')

### comparing to curves
tog$time_prediction = sapply(tog$f_index, time_formula)
tog$beh_prediction = sapply(tog$f_index, beh_formula)
tog$test_prediction = sapply(tog$f_index, test_formula)
tog$e2_prediction = sapply(tog$f_index, e2_formula)
tog$kt_prediction = sapply(tog$f_index, kt_formula)
tog$mass_prediction = sapply(tog$f_index, mass_formula)

ggplot(subset(tog, final_clusters ==1), aes(x = f_index, y = scale_2(prop)))+
  geom_point(aes(color = Status))+
  labs(title = 'Proportion of 1')+
  geom_line(aes(color = 'inverted % testicular',x = f_index, y = 1+(-1*scale_2(test_prediction))))+
    geom_line(aes(color = 'e2',x = f_index, y = (scale_2(e2_prediction))))+
      geom_line(aes(color = 'kt',x = f_index, y = (scale_2(kt_prediction))))

ggplot(subset(tog, final_clusters ==1), aes(x = f_index, y = (prop)/max(prop)))+
  geom_point(aes(color = Status))+
  labs(title = 'Proportion of 1')+
  geom_line(aes(color = 'inverted % testicular',x = f_index, y = 1+(-1*(test_prediction)/max(test_prediction, na.rm =T))))+
    geom_line(aes(color = 'e2',x = f_index, y = ((e2_prediction)/max(e2_prediction, na.rm = T))))+
      geom_line(aes(color = 'kt',x = f_index, y = kt_prediction/max(kt_prediction, na.rm =T)))

ggplot(subset(tog, final_clusters ==24), aes(x = f_index, y = scale_2(prop)))+
  geom_point(aes(color = Status))+
  labs(title = 'Proportion of 24')+
  geom_line(aes(color = 'inverted % testicular',x = f_index, y = 1+(-1*scale_2(test_prediction))))+
    geom_line(aes(color = 'e2',x = f_index, y = (scale_2(e2_prediction))))+
      geom_line(aes(color = 'kt',x = f_index, y = (scale_2(kt_prediction))))

ggplot(subset(tog, final_clusters ==24), aes(x = f_index, y = scale_2(prop)))+
 # geom_point(aes(color = Status))+
    geom_line(aes(color = 'e2',x = f_index, y = ((e2_prediction))))+
      geom_line(aes(color = 'mass',x = f_index, y = ((mass_prediction))))+
        geom_line(aes(color = 'kt',x = f_index, y = ((kt_prediction))))+
          geom_line(aes(color = 'testicular',x = f_index, y = ((test_prediction))))+
            geom_line(aes(color = 'behavior',x = f_index, y = ((beh_prediction))/600))+
              geom_line(aes(color = 'time',x = f_index, y = ((time_prediction))/600))













