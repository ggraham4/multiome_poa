obj  = readRDS("~/Desktop/optimal_clustering_rna_only.rds")
sub_6 = subset(obj,res0.8_50nn_40PC_45LSI==6)

sub_6$gnrh1 = sub_6@assays$RNA$data['LOC111571064',]

gnrh_dat =sub_6@meta.data%>%
  group_by(individual,Status)%>%
  summarize(mean_gnrh = mean(gnrh1))%>%
  subset(Status!='NRM')
mod =lm(mean_gnrh~Status, data = gnrh_dat)

anova(mod, test ='Chisq')
pairs(emmeans(mod, 'Status'), adjust ='none')
# pairs are significant BUT type 3 is not so we ignore