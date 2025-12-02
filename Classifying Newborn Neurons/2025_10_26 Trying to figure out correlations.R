marks_9_22 = FindMarkers(obj, ident.1 = 9, ident.2 = 22)

clown_go(rownames(marks_9_22[marks_9_22$p_val_adj<0.05,]))%>%dotplot()

obj_9_22 = subset(obj, final_clusters %in% c(9, 22, 6))

matrix_9_22 = obj_9_22@assays$RNA$data%>% as.matrix()

cyto = CytoTRACE(matrix_9_22)

obj_9_22$cyto =cyto$CytoTRACE

ggplot(obj_9_22@meta.data, aes(x = final_clusters, y = cyto))+
  geom_boxplot()

marks_9 = FindMarkers(obj, ident.1 = 9, only.pos = T)
clown_go(rownames(marks_9[marks_9$p_val_adj<0.05,]))%>%dotplot()

DotPlot(obj, 'mex3a')

obj_2 = FindSubCluster(obj, 1, graph.name = 'harmony.wsnn')
obj_3 = FindSubCluster(obj, 6, graph.name = 'harmony.wsnn')
obj_2$sub.cluster = as.character(obj_2$sub.cluster)
obj_2$sub.cluster = ifelse(obj_2$final_clusters ==6,obj_3$sub.cluster, obj_2$sub.cluster)
Idents(obj_2) <- 'sub.cluster'
DimPlot(obj_2)

dap2 = obj_2@meta.data%>%
  group_by(Status, individual)%>%
  summarize(n = n())%>%
  select(individual, n)

dap=obj_2@meta.data%>%
  group_by(Status, individual, sub.cluster)%>%
  summarize(n_cells = n())%>%
  right_join(dap2, by = 'individual')%>%
  mutate(prop = n_cells/n)

dap$Status.x = factor(dap$Status.x, levels = c('NRM','M','D','E','NF','F'))

ggplot(subset(dap, sub.cluster %in%c('1_2', '9', '6_3'
                                     #,'22'
                                     ) ), aes(x = factor(sub.cluster, levels = c('1_2',
                                                                                                     '9',
                                                                                                     #'22',
                                                                                                     '6_3')), y = n_cells, group = individual))+
  geom_point()+
  geom_line(aes(color = individual))

ggplot(subset(dap, sub.cluster %in%c('1_2', '9', '6_3') ), aes(x = factor(sub.cluster, levels = c('1_2','9','6_3')),
                                                             y = n_cells, group = Status.x, color = Status.x))+
  geom_point()+
  geom_smooth(aes(color = Status.x), method ='lm', se = F)

ggplot(subset(dap, sub.cluster %in%c('1_2', '9', '6_3') ), aes(x = factor(sub.cluster, levels = c('1_2','9','6_3')),
                                                             y = n_cells, group = individual))+
  geom_point()+
  geom_line(aes(color = individual))

ggplot(subset(dap, sub.cluster %in%c('1_2','22', '9', '6_3') ), aes(x = factor(sub.cluster, levels = c('1_2','22','9','6_3')),
                                                             y = prop, group = individual))+
  geom_point()+
  geom_line(aes(color = individual))

ggplot(subset(dap, sub.cluster %in%c('1_2','22', '9', '6_2') ), aes(x = factor(sub.cluster, levels = c('1_2','22','9','6_3')),
                                                             y = prop, fill = Status.x))+
  geom_boxplot()

ggplot(subset(dap, sub.cluster %in%c('1_1','1_3', '1_2')&Status.x %in% c('M','D','F')  ),
       aes(x = factor(sub.cluster, levels = c('1_1','1_3', '1_2')),
           y = prop))+
  geom_boxplot(aes(fill = Status.x))+
    stat_summary(geom = 'line', fun = 'mean', aes(group = Status.x), color = 'black', size = 2)+
  stat_summary(geom = 'line', fun = 'mean', aes(group = Status.x, color = Status.x))

ggplot(subset(dap, sub.cluster %in%c('1_1','1_3', '1_2', '6_3')), aes(x = factor(sub.cluster, levels = c('1_1','1_3', '1_2', '6_3')),
                                                             y = prop, group = individual))+
    geom_line(aes(color = individual))

ggplot(subset(dap, sub.cluster %in%c('1_1','1_3', '1_2', '6_3') ), aes(x = factor(sub.cluster, levels = c('1_1','1_3', '1_2', '6_3')),
                                                             y = prop, group = Status.x))+
    geom_smooth(aes(color = Status.x), se = F)

subset_dap =subset(dap, sub.cluster %in%c('1_2','22', '9', '6_3') ) 

prop_wide = dap %>%
  select(individual, sub.cluster, prop) %>%
  pivot_wider(names_from = sub.cluster, values_from = prop,  values_fill = 0)

cor_matrix = (cor((prop_wide[, -1]), use = "pairwise.complete.obs"))

cor_matrix
diag(cor_matrix) =0
heatmap(cor_matrix)

ggplot(prop_wide, aes(x = `9`, y = `1_2`, color = individual))+
  geom_point()

m = lm(data = prop_wide, `9` ~`1_2`)
summary(m)

m = lm(data = prop_wide, `9` ~`26`)
summary(m)

m = lm(data = prop_wide, `6_3` ~`9`)
summary(m)
m = lm(data = prop_wide, `6_3` ~`22`)
summary(m)

m = lm(data = prop_wide, `6_2` ~`9`)
summary(m)
m = lm(data = prop_wide, `6_2` ~`22`)
summary(m)

m = lm(data = prop_wide, `6_3` ~`1_2`)
summary(m)
m = lm(data = prop_wide, `6_2` ~`1_2`)
summary(m)

m = lm(data = prop_wide, `6_3` ~`1_3`)
summary(m)
m = lm(data = prop_wide, `6_2` ~`1_3`)
summary(m)

ggplot(subset(dap, sub.cluster %in%c('1_1','1_3','1_2', '6_2', '6_3')& Status.x %in%c('M','D',"F")), aes(x = factor(sub.cluster, levels = c('1_1','1_3','1_2', '6_2', '6_3')),
                                                             y = prop))+
  geom_boxplot(aes(fill = Status.x))+
  geom_smooth(method = 'loess', se =F, aes(color = Status.x, group = Status.x))
# 1_2 and 6_2 are the most highly correlated


## partial correlation
library(ppcor)
pcor.test(prop_wide$`6_3`, prop_wide$`1_2`, prop_wide$`9`)  # control for 9
pcor.test(prop_wide$`6_3`, prop_wide$`1_2`, prop_wide$`22`) # control for 22


summary(lm(`9` ~ `1_2`, data = prop_wide))
summary(lm(`22` ~ `1_2`, data = prop_wide))

# Step 2: Does 9 or 22 predict 6_3 after accounting for 1_2?
summary(lm(`6_3` ~ `1_2` + `9`, data = prop_wide))
summary(lm(`6_3` ~ `1_2` + `22`, data = prop_wide))

m9  <- lm(`6_3` ~ `1_2` + `9`, data = prop_wide)
m22 <- lm(`6_3` ~ `1_2` + `22`, data = prop_wide)
m9_22 <- lm(`6_3` ~ `1_2` + `9` + `22`, data = prop_wide)

AIC(m9, m22, m9_22)
summary(m9)$adj.r.squared
summary(m22)$adj.r.squared
summary(m9_22)$adj.r.squared

library(bnlearn)
network_data = prop_wide[, c("1_1",'1_3', "1_2", "9", "22", "6_3", '6_2', '6_1','6_0')]
bn = hc(network_data)
plot(bn)

### im realizing its kind of silly to try to valdate my developmental lineage with proportions 
# when that is not how I came up with it, I used gene expression...

