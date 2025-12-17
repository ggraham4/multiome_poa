library(tidyverse)
library(mgcv)
library(ggplot2)
library(factoextra)


coalesced = read.csv('Measures/coalesced_data.csv')

ggplot(coalesced, aes(x = mass_final, y = 10^Log_11KT))+
  geom_point(aes(color = Phase))+
  geom_smooth()

ggplot(coalesced, aes(x = mass_final, y = 10^Log_11KT))+
  geom_point(aes(color = Phase))

ggplot(coalesced, aes(x = mass_final, y = 10^Log_11KT))+
  geom_point(aes(color = Phase))+
  geom_smooth(aes(color=Phase), method= 'gam')

ggplot(coalesced, aes(x = mass_final, y = 10^Log_11KT))+
  geom_point(aes(color = Phase))+
  geom_smooth(aes(color=Phase), method= 'glm')

ggplot(coalesced, aes(x = mass_final, y = 10^Log_E2))+
  geom_point(aes(color = Phase))+
  geom_smooth()

ggplot(coalesced, aes(x = mass_final, y = Log_E2))+
  geom_point(aes(color = Phase))+
    geom_smooth(aes(color=Phase), method= 'gam')

ggplot(coalesced, aes(x = mass_final, y = Log_11KT))+
  geom_point(aes(color = Phase))+
  geom_smooth(aes(color=Phase), method= 'gam')

ggplot(subset(coalesced, Phase != 'S'), aes(x = mass_final, y = 10^Log_11KT))+
  geom_point(aes(color = Phase))+
  geom_smooth()

ggplot(subset(coalesced, Phase != 'S'), aes(x = mass_final, y = Log_11KT))+
  geom_point(aes(color = Phase))+
  geom_smooth()

ggplot(coalesced, aes(x = 10^Log_11KT, y = 10^Log_E2))+
  geom_point(aes(shape = Phase))+
  geom_smooth()

ggplot(coalesced, aes(x = Log_11KT, y = Log_E2))+
  geom_point(aes(shape = Phase))+
  geom_smooth(aes(color=Phase), method= 'gam')

ggplot(coalesced, aes(x = Log_11KT, y = Log_E2))+
  geom_point(aes(color = Phase))+
  geom_smooth()

ggplot(coalesced, aes(x = mass_final, y = Percent_Testicular))+
  geom_point(aes(color = Phase))+
  geom_smooth()

ggplot(coalesced, aes(x = mass_final, y = Percent_Testicular))+
  geom_point(aes(color = Phase))+
  geom_smooth(aes(color=Phase), method= 'gam')

ggplot(coalesced, aes(x = Log_E2, y = Percent_Testicular))+
  geom_point(aes(shape = Phase))
# no relationship at all

ggplot(coalesced, aes(x = Log_11KT, y = Percent_Testicular))+
  geom_point(aes(color = Phase))

### theory time #### 

# pca 1
vars1 = c('Percent_Testicular',
         'Log_11KT',
         #'Log_E2',
         'mass_final'
         )

pc_data1 = coalesced[complete.cases(coalesced[,vars1]),]%>%
  subset(!Phase %in% c('S',
                       'EP',
                       'IDK',
                       'NM'))
X1 <- scale(pc_data1[, vars1])
pca1 <- prcomp(X1, center = TRUE, scale. = TRUE)
fviz_pca_biplot(pca1, habillage = pc_data1$Phase)

pc_data1$pc1 <- pca$x[,1]
pc_data1$pc2 <- pca$x[,2] *-1

pc_data1$Phase= factor(pc_data$Phase, levels = c('S','EP','NM','M','D','E','NF','F'))
ggplot(pc_data1, aes(x = Phase, y = pc1))+
  geom_boxplot()+
  geom_jitter()

ggplot(pc_data1, aes(x = Phase, y = pc2))+
  geom_boxplot()+
  geom_jitter()


z1_1 = pca1$x[,1]
z2_1 = pca1$x[,2]

#> here , pc1 explains sex change extrordinarily well, but not pc2 
#pc1 loads mass and 11kt

###pca2
vars2 = c(#'Percent_Testicular',
          'Log_11KT',
          'Log_E2',
          'mass_final'
)

pc_data2 = coalesced[complete.cases(coalesced[,vars2]),]%>%
  subset(!Phase %in% c('S',
                       'EP',
                       'IDK',
                       'NM'))
X2 <- scale(pc_data2[, vars2])
pca2 <- prcomp(X2, center = TRUE, scale. = TRUE)
fviz_pca_biplot(pca2, habillage = pc_data2$Phase)

pc_data2$pc1 <- pca2$x[,1]
pc_data2$pc2 <- pca2$x[,2] *-1

pc_data2$Phase= factor(pc_data2$Phase, levels = c('S','EP','NM','M','D','E','NF','F'))
ggplot(pc_data2, aes(x = Phase, y = pc1))+
  geom_boxplot()+
  geom_jitter()

ggplot(pc_data2, aes(x = Phase, y = pc2))+
  geom_boxplot()+
  geom_jitter()


z1_2 = pca1$x[,1]
z2_2 = pca1$x[,2]

# again, 1 explains sex change and not 2
#pc1 loads mass and 11kt

# parsimonious pc3
vars3= c(#'Percent_Testicular',
  'Log_11KT',
  #'Log_E2',
  'mass_final'
)

pc_data3 = coalesced[complete.cases(coalesced[,vars3]),]%>%
  subset(!Phase %in% c('S',
                       'EP',
                       'IDK',
                       'NM'))
X3 <- scale(pc_data3[, vars3])
pca3 <- prcomp(X3, center = TRUE, scale. = TRUE)
fviz_pca_biplot(pca3, habillage = pc_data3$Phase)

pc_data3$pc1 <- pca3$x[,1]
pc_data3$pc2 <- pca3$x[,2] *-1

pc_data3$Phase= factor(pc_data3$Phase, levels = c('S','EP','NM','M','D','E','NF','F'))
ggplot(pc_data3, aes(x = Phase, y = pc1))+
  geom_boxplot()+
  geom_jitter()

ggplot(pc_data3, aes(x = Phase, y = pc2))+
  geom_boxplot()+
  geom_jitter()
# yeah just pc1 again


#some more theory
pc_data3$Z = pc_data3$pc1

mod1 = lm(Z~Phase, pc_data3)
anova(mod1, test = 'Chisq')

mod2 = gam(Z~Phase, data =pc_data3)
summary(mod2)
anova(mod2)

mod3 =gam(Z~Log_E2,data= pc_data3)
plot(mod3, all.terms=T)

ggplot(pc_data3, aes(x =Z, y = Log_E2))+
  geom_point(aes(color=Phase))

ggplot(pc_data3, aes(x =Z, y = Percent_Testicular))+
  geom_point(aes(color=Phase))

#modeling
testicular_nls_data = subset(pc_data3, !is.na(Percent_Testicular))

nls_testis <- nls(
  Percent_Testicular ~ L / (1 + exp(k * (Z - Z0))),
  data = testicular_nls_data,
  start = list(
    L = max(testicular_nls_data$Percent_Testicular, na.rm = TRUE),
    k = 2,
    Z0 = median(testicular_nls_data$Z, na.rm = TRUE)
  ),
  control = nls.control(maxiter = 200)
)
summary(nls_testis)

nls_tesis_curve = function(Z){
  output =  2.6679/(1 + exp(0.8805 * (Z - (-2.8581))))
  return(output)
}

pc_data3$testis_curve = sapply( pc_data3$Z,nls_tesis_curve)

ggplot(pc_data3, aes(x =Z, y = Percent_Testicular))+
  geom_point(aes(color=Phase))+
  geom_line(aes(x = Z, y = testis_curve, color = 'NLS Line'), size =2, color ='blue')




