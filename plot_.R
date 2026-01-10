library(Seurat)
library(tidyverse)
library(lme4)
library(ggplot2)

latent = read.csv('Measures/2025_12_26 all_data.csv')
obj <- readRDS("~/Desktop/nemo.orig_harmony.integration_all_testd_clusters.rds") 

rownames = rownames(obj@meta.data)
obj@meta.data=obj@meta.data%>%
  left_join(latent, by = join_by('individual'=='Fish'))
rownames(obj@meta.data) = rownames

obj$aromatase = obj@assays$RNA$data['LOC111577263',]
dat = obj@meta.data%>%
  subset(res0.8_50nn_40PC_45LSI==1)%>%
  group_by(individual, Status.x)%>%
  summarize(mean_arom = mean(aromatase),
            mean_sex = mean(SexState))

dat$Status.x = factor(dat$Status.x, levels = c('NRM',
                                                                   'M',
                                                                   'D',
                                                                   'E',
                                                                   'NF',
                                                                   'F'))



ggplot(subset(dat, Status.x != 'NRM'), aes(x = mean_sex, y = scale(mean_arom), color = Status.x))+
  geom_point()

plot_ = function(gene, cluster){
  temp = obj@meta.data
  temp$data = obj@assays$RNA$data[gene,]
  
  dat = temp%>%
  subset(res0.8_50nn_40PC_45LSI==cluster)%>%
  group_by(individual, Status.x)%>%
  summarize(mean_gene = mean(data),
            mean_sex = mean(SexState))
  
  p = ggplot(subset(dat, Status.x != 'NRM'), aes(x = mean_sex, y = scale(mean_gene), color = Status.x))+
  geom_point()+
    labs(title = gene, subtitle = cluster)

return(p)
}

plot_('tacr3a', 6)

# prime triangel candidate
obj$tacr3a = obj@assays$RNA$data['tacr3a',]
tacr3a_ = obj@meta.data%>%
  subset(res0.8_50nn_40PC_45LSI==6)%>%
  group_by(individual, Status.x)%>%
  summarize(mean_tacr3a = mean(tacr3a),
            mean_sex = mean(SexState))


mod1 = lm(mean_tacr3a~mean_sex, data = subset(tacr3a_, Status.x %in%c('M','D')))
mod2 = lm(mean_tacr3a~mean_sex, data = subset(tacr3a_, Status.x %in%c('D','F')))

m1 = mod1$coefficients[2]
b1 = mod1$coefficients[1]

m2 = mod2$coefficients[2]
b2 = mod2$coefficients[1]

#intersection 
offset = mean(tacr3a_$mean_sex[tacr3a_$Status.x=='D'])
#calculate intersection
#(m1x)+ b1 = (m2x)+b2
#(m1-m2)*x =-b1+b2
xmid = (-b1+b2)/(m1-m2)


ggplot(subset(tacr3a_, Status.x != 'NRM'), aes(x = mean_sex, y =as.numeric(mean_tacr3a), color = Status.x))+
  geom_point()+
  geom_line(data= data.frame(x = seq(-1,2, by =0.1)),aes(x =x , y = (m1*x)+b1), inherit.aes = F)+
  geom_line(data= data.frame(x = seq(-1,2, by =0.1)),aes(x =x , y = (m2*x)+b2), inherit.aes = F)+
  geom_vline(xintercept = xmid)

#### ok lets consider this conceptually ####

#>When F is present
#>m1 + Gravity + Ff > alpha
#>when F is absent
#>m1+ Gravity < alpha
#>
#> where
#> m1 is the slope of the line
#> Constant can be thought of as gravity
#> Ff is the force the female is applying against the slope
#> Ff == 
#> qlphq = force necessary to overcome  I
#> 
#> delta_g = male_mean (gex) - female mean (gex)
#> 

#if m1 + gravity < alpha 
# we can pin gravity to 1
gravity = 1
m1 = m1 # 0.87
alpha = 1.09
m1 + gravity < alpha
Ff = 0.1
m1 + gravity + Ff > alpha

# but really we can consider any number of values to be valid in this system
#equation 1
#G + m1 < alpha
# m1 < alpha - G
# 0.87 < x - y
# y = x-0.87
# G < alpha -m1


#equation 2
#m1+gravity + FF > alpha
# m1 > alpha - G - FF
# 0.87 > x - y - z
#G > alpha - m1 - Ff
#y = -m1 + x - z

# together
#a - m1 - Ff < G < a - m1

# G < a - m1
# G > 0  & a > 0 & Ff > 0 

# for shits and giggles lets pin G = 1
#G = 1
#G>0
#1 < a-m1
G =1
a_low =m1+G

#G > alpha - m1 - Ff
#G > a_low- m1 - Ff
# Ff > a_low - m1 - g
Ff_low =a_low - m1 - G

ggplot(data.frame(x = seq(-100, 100, by = 0.1)), aes(x = x, y = x))+
  geom_point(size = 0, alpha = 0)+
  geom_abline( linetype = 'dashed', aes(slope = m1, intercept  = G,color = 'values for a, F absent above'))+
  geom_abline( linetype = 'dashed', aes(slope = -m1, intercept= -G-a_low,color = 'values for Ff, F present above') )


