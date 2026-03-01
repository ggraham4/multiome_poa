library(tidyverse)
library(mgcv)
library(ggplot2)
library(factoextra)
library(patchwork)
library(ggeffects)

coalesced = read.csv('Measures/coalesced_data.csv')
multiome_measures = read.csv('Measures/all_data.csv')



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

pc_data1$pc1 <- pca1$x[,1]
pc_data1$pc2 <- pca1$x[,2] *-1

pc_data1$Phase= factor(pc_data1$Phase, levels = c('S','EP','NM','M','D','E','NF','F'))
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

# parsimonious pc4
vars4= c(#'Percent_Testicular',
  'Log_11KT',
  'Log_E2'
)

##kt + e2
pc_data4 = coalesced[complete.cases(coalesced[,vars4]),]%>%
  subset(!Phase %in% c('S',
                       'EP',
                       'IDK',
                       'NM'))
X4 <- scale(pc_data4[, vars4])
pca4 <- prcomp(X4, center = TRUE, scale. = TRUE)
fviz_pca_biplot(pca4, habillage = pc_data4$Phase)

pc_data4$pc1 <- pca4$x[,1]
pc_data4$pc2 <- pca4$x[,2] *-1

pc_data4$Phase= factor(pc_data4$Phase, levels = c('S','EP','NM','M','D','E','NF','F'))
ggplot(pc_data4, aes(x = Phase, y = pc1))+
  geom_boxplot()+
  geom_jitter()

ggplot(pc_data4, aes(x = Phase, y = pc2))+
  geom_boxplot()+
  geom_jitter()



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


a=ggplot(pc_data3, aes(x =Z, y = Log_11KT))+
  geom_point(aes(color=Phase))

b=ggplot(pc_data3, aes(x =Z, y = Log_E2))+
  geom_point(aes(color=Phase))

pc_data3$kt_e2 = pc_data3$Log_11KT/pc_data3$Log_E2
c=ggplot(pc_data3, aes(x =Z, y = kt_e2))+
  geom_point(aes(color=Phase))+
  labs(y = "LogKT / Log E2")

d=ggplot(pc_data3, aes(x =Z, y = mass_final))+
  geom_point(aes(color=Phase))+
  labs(y = "Mass (g)")

e=ggplot(pc_data3, aes(x =Z, y = Percent_Testicular))+
  geom_point(aes(color=Phase))

female_z = mean(pc_data3$Z[pc_data3$Phase=='F'])
dom_z = mean(pc_data3$Z[pc_data3$Phase=='D'])

f=ggplot(pc_data3, aes(x =Z, y = days))+
  labs(y = "Days Since Pairing")+
  geom_point(aes(color=Phase))#
# geom_vline(xintercept  = female_z)

dat_behavior = pc_data3%>%
  right_join(multiome_measures, by = 'Log_11KT')

g=ggplot(dat_behavior, aes(x =Z, y = Behaviors_Day_2))+
  labs(y = "Parental Behaviors")+
  geom_point(aes(color=Phase))#


a+b+c+d+e+g

# fit a line to measures #
# 11kt
kt_line = nls(Z~m*Log_11KT+b, data = pc_data3,
                start = list(
    m =-1,
    b = 4
  ),
  control = nls.control(maxiter = 200)
)

summary(kt_line)

kt_data = pc_data3
kt_data$nls = -1.946*kt_data$Log_11KT+4.561
pc_data3$kt_line = kt_data$nls

ggplot(kt_data, aes(x =Log_11KT, y =Z ))+
  geom_point(aes(color=Phase))+
  geom_line(aes(x = Log_11KT, y = nls))

# compare to linear regression
kt_model = lm(Z~Log_11KT, data = pc_data3)
summary(kt_model)
#exact same 

#mass
mass_model = lm(Z~mass_final, data = pc_data3)
summary(mass_model)

pc_data3$mass_line = (mass_model$coefficients[2]*pc_data3$mass_final)+mass_model$coefficients[1]

ggplot(pc_data3, aes(x =mass_final, y =Z ))+
  geom_point(aes(color=Phase))+
  geom_line(aes(x = mass_final, y = mass_line))

#e2 
e2_model = lm(Z~Log_E2, data = pc_data3)
summary(e2_model)

pc_data3$e2_line = (e2_model$coefficients[2]*pc_data3$Log_E2)+e2_model$coefficients[1]
#z = 1.4x*E2 -3.7

ggplot(subset(pc_data3, !is.na(Log_E2)), aes(x =Log_E2, y =Z ))+
  geom_point(aes(color=Phase))+
  geom_line(data =subset(pc_data3, !is.na(Log_E2)), aes(x = Log_E2, y = e2_line))

#kt_e2
kt_e2_model = lm(Z~kt_e2, data = pc_data3)
summary(kt_e2_model)

pc_data3$kt_e2_line = (kt_e2_model$coefficients[2]*pc_data3$kt_e2)+kt_e2_model$coefficients[1]
#z = -2.4x*E2 2.4

ggplot(subset(pc_data3, !is.na(Log_E2)), aes(x =kt_e2, y =Z ))+
  geom_point(aes(color=Phase))+
  geom_line(data =subset(pc_data3, !is.na(Log_E2)), aes(x = kt_e2, y = kt_e2_line))

# trying to fit to testicular
ggplot(pc_data3, aes(x =Percent_Testicular, y =Z ))+
  geom_point(aes(color=Phase))

test_line <- nls(
  Z ~ a - b * log(Percent_Testicular),
  data = subset(pc_data3, Percent_Testicular > 0),
  start = list(a = 0, b = 1)
)

summary(test_line)
test_line_coef = summary(test_line)$coefficients

pc_data3$test_line = test_line_coef[1,1] -  test_line_coef[2,1] * log(pc_data3$Percent_Testicular)

ggplot(pc_data3, aes(x =Percent_Testicular, y =Z ))+
  geom_point(aes(color=Phase))+
    geom_line( aes(x = Percent_Testicular, y = test_line))

## behavior
ggplot(dat_behavior, aes(x =Behaviors_Day_2, y =Z ))+
  geom_point(aes(color=Phase))

beh_line <- nls(
  Z ~ a - b * log(Behaviors_Day_2),
  data = subset(dat_behavior),
  start = list(a = 0, b = 1)
)
beh_line_coef = summary(beh_line)$coef

dat_behavior$beh_line = beh_line_coef[1,1]-beh_line_coef[2,1]*log(dat_behavior$Behaviors_Day_2)

ggplot(dat_behavior, aes(x =Behaviors_Day_2, y =Z ))+
  geom_point(aes(color=Phase))+
    geom_line( aes(x = Behaviors_Day_2, y = beh_line))


c_start <- with(dat_behavior,
                mean(Z[Behaviors_Day_2 > 500], na.rm = TRUE))
a_start <- with(dat_behavior,
                mean(Z[Behaviors_Day_2 < 100], na.rm = TRUE) - c_start)
b_start <- 1 / 250   # ≈ 0.004

beh_line_2 <- nls(
  Z ~ a * exp(-b * Behaviors_Day_2) + c,
  data = subset(dat_behavior, Behaviors_Day_2 >= 0),
  start = list(
    a = a_start,
    b = b_start,
    c = c_start
  ),
  control = nls.control(maxiter = 200)
)
beh_line_2_coef = summary(beh_line_2)$coefficients
summary(beh_line_2)

dat_behavior$beh_line_2 = beh_line_2_coef[1,1] *exp(-beh_line_2_coef[2,1]*dat_behavior$Behaviors_Day_2)+  beh_line_2_coef[3,1]

ggplot(dat_behavior, aes(x =Behaviors_Day_2, y =Z ))+
  geom_point(aes(color=Phase))+
    geom_line( aes(x = Behaviors_Day_2, y = beh_line_2))


anova(beh_line, beh_line_2)
# anova says no difference and I think equation 1 is better


a_fit = ggplot(kt_data, aes(x =Log_11KT, y =Z ))+
  geom_point(aes(color=Phase))+
  geom_line(aes(x = Log_11KT, y = nls))

b_fit = ggplot(pc_data3, aes(x =mass_final, y =Z ))+
  geom_point(aes(color=Phase))+
  geom_line(aes(x = mass_final, y = mass_line))

c_fit = ggplot(subset(pc_data3, !is.na(Log_E2)), aes(x =Log_E2, y =Z ))+
  geom_point(aes(color=Phase))+
  geom_line(data =subset(pc_data3, !is.na(Log_E2)), aes(x = Log_E2, y = e2_line))

d_fit = ggplot(subset(pc_data3, !is.na(Log_E2)), aes(x =kt_e2, y =Z ))+
  geom_point(aes(color=Phase))+
  geom_line(data =subset(pc_data3, !is.na(Log_E2)), aes(x = kt_e2, y = kt_e2_line))

e_fit =ggplot(pc_data3, aes(x =Percent_Testicular, y =Z ))+
  geom_point(aes(color=Phase))+
    geom_line( aes(x = Percent_Testicular, y = test_line))

f_fit =ggplot(dat_behavior, aes(x =Behaviors_Day_2, y =Z ))+
  geom_point(aes(color=Phase))+
    geom_line( aes(x = Behaviors_Day_2, y = beh_line))

a_fit+b_fit+c_fit+d_fit+e_fit+f_fit

#time
time_model = lm(Z~Time_Day_2, data = dat_behavior)
summary(time_model)

dat_behavior$time_line = (time_model$coefficients[2]*dat_behavior$Time_Day_2)+time_model$coefficients[1]

ggplot(dat_behavior, aes(x =Time_Day_2, y =Z ))+
  geom_point(aes(color=Phase))+
  geom_line(aes(x = Time_Day_2, y = time_line))

### other misc recordings
ggplot(dat_behavior, aes(x =log(Estimated_Volume_2.5x+0.1), y =Z ))+
  geom_point(aes(color=Phase))

ggplot(dat_behavior, aes(x =Behaviors_Per_Time, y =Z ))+
  geom_point(aes(color=Phase))

ggplot(dat_behavior, aes(x =Percent_Ovarian, y =Z ))+
  geom_point(aes(color=Phase))

ggplot(dat_behavior, aes(x =Percent_Ovarian, y =Percent_Testicular.x ))+
  geom_point(aes(color=Phase))


### together
#> so right now all of my lines are predicting z, but now I think I want z to predict them

#z =-1.946*kt_data$Log_11KT+4.561 
#kt = (z-4.561)/-1.946
kt_formula = function(Z){
  out = (Z-4.561)/-1.946
  return(out)
  }
mass_formula = function(Z){
  #z=(mass_model$coefficients[2]*pc_data3$mass_final)+mass_model$coefficients[1]
  out = (Z-mass_model$coefficients[1])/mass_model$coefficients[2]
  return(out)
}
e2_formula = function(Z){
    out = (Z-e2_model$coefficients[1])/e2_model$coefficients[2]
  return(out)
}
kt_e2_formula = function(Z){
    out = (Z-kt_e2_model$coefficients[1])/kt_e2_model$coefficients[2]
  return(out)
}
test_formula = function(Z){
 a <- test_line_coef[1,1]
  b <- test_line_coef[2,1]
  out <- exp((a - Z) / b)
  return(out)
}
beh_formula = function(Z){
   a <- beh_line_coef[1,1]
  b <- beh_line_coef[2,1]
  out <- exp((a - Z) / b)
  return(out)
}
time_formula = function(Z){
  #z = (time_model$coefficients[2]*dat_behavior$Time_Day_2)+time_model$coefficients[1]
  #(z-time_coef1).time_coef2
  out = ((Z-time_model$coefficients[1] )/time_model$coefficients[2])
  return(out)
}


line_data = data.frame(Z = seq(range(pc_data3$Z)[1], range(pc_data3$Z)[2], by =0.1))

line_data$kt =sapply(line_data$Z, kt_formula)
line_data$e2 =sapply(line_data$Z, e2_formula)
line_data$kt_e2 =sapply(line_data$Z, kt_e2_formula)
line_data$testicular =sapply(line_data$Z, test_formula)
line_data$behavior =sapply(line_data$Z, beh_formula)
line_data$time =sapply(line_data$Z, time_formula)
line_data$mass =sapply(line_data$Z, mass_formula)

scale_2 <- function(x) {
  x <- as.numeric(x)
  rng <- range(x, na.rm = TRUE)
  (x - rng[1]) / (rng[2] - rng[1])
}

ggplot(line_data, aes(x = Z, y = scale_2(10^(kt))))+
  geom_line(aes(color = '11KT (pg/ml)'))+
  geom_line(aes(color = 'E2 (pg/ml)', y = scale_2(10^(e2))))+
    geom_line(aes(color = 'Mass (g)', y = scale_2((mass))))+
 # geom_line(aes(color = 'KT/ E2', y = scale(kt_e2)))+
  geom_line(data = subset(line_data, testicular <= 1),aes(color = 'Testicular %', y = scale_2(testicular)))+
  geom_line(data = subset(line_data, behavior <= 600),aes(color = 'Behaviors (10m)', y = scale_2(behavior)))+
  geom_line(data = subset(line_data, time <= 600),aes(color = 'Time in Nest (10m)', y = scale_2(time)))+
  labs(y= 'Scaled Value', x = 'F Index', color = 'Measurement')+
  theme_minimal()

ggplot(line_data, aes(x = Z, y = scale_2((kt))))+
  geom_line(aes(color = 'Log10 11KT (pg/ml)'))+
  geom_line( aes(linetype = 'Log10 E2 (pg/ml)', color = 'Log10 E2 (pg/ml)', y = scale_2((e2))))+
      geom_line(color = 'black',aes(linetype = 'Mass (g)',color = 'Mass (g)', y = scale_2((mass))))+
 # geom_line(aes(color = 'KT/ E2', y = scale(kt_e2)))+
  geom_line(data = subset(line_data, testicular <= 1),aes(color = 'Testicular %', y = scale_2(testicular)))+
  geom_line(data = subset(line_data, behavior <= 600),aes(color = 'Behaviors (10m)', y = scale_2(behavior)))+
  geom_line(data = subset(line_data, time <= 600),aes(color = 'Time in Nest (10m)', y = scale_2(time)))+
  labs(y= 'Scaled Value', x = 'F Index', color = 'Measurement')+
  theme_minimal()

pc_data3$Phase= factor(pc_data3$Phase, levels = c('S','EP','NM','M','D','E','NF','F'))
boxes = ggplot(pc_data3, aes(x = Phase, y = pc1))+
  geom_boxplot()+
  geom_jitter()+coord_flip()+
  theme_minimal()


exponentiated = ggplot(line_data, aes(x = Z, y = scale_2(exp(kt))))+
  geom_line(aes(color = '11KT (pg/ml)'))+
  geom_line(aes(color = 'E2 (pg/ml)', y = scale_2(exp(e2))))+
    geom_line(aes(color = 'Mass (g)', y = scale_2((mass))))+
 # geom_line(aes(color = 'KT/ E2', y = scale(kt_e2)))+
  geom_line(data = subset(line_data, testicular <= 1),aes(color = 'Testicular %', y = scale_2(testicular)))+
  geom_line(data = subset(line_data, behavior <= 600),aes(color = 'Behaviors (10m)', y = scale_2(behavior)))+
  geom_line(data = subset(line_data, time <= 600),aes(color = 'Time in Nest (10m)', y = scale_2(time)))+
  labs(y= 'Scaled Value', x = 'F Index', color = 'Measurement')+
  theme_minimal()

(exponentiated+boxes)+
  plot_layout(ncol = 1)

logged = ggplot(line_data, aes(x = Z, y = scale_2((kt))))+
  geom_line(aes(color = 'Log10 11KT (pg/ml)'))+
  geom_line( aes(linetype = 'Log10 E2 (pg/ml)', color = 'Log10 E2 (pg/ml)', y = scale_2((e2))))+
      geom_line(color = 'black',aes(linetype = 'Mass (g)',color = 'Mass (g)', y = scale_2((mass))))+
 # geom_line(aes(color = 'KT/ E2', y = scale(kt_e2)))+
  geom_line(data = subset(line_data, testicular <= 1),aes(color = 'Testicular %', y = scale_2(testicular)))+
  geom_line(data = subset(line_data, behavior <= 600),aes(color = 'Behaviors (10m)', y = scale_2(behavior)))+
  geom_line(data = subset(line_data, time <= 600),aes(color = 'Time in Nest (10m)', y = scale_2(time)))+
  labs(y= 'Scaled Value', x = 'F Index', color = 'Measurement')+
  theme_minimal()
(logged+boxes)+
  plot_layout(ncol = 1)

### save_functions
time_formula = function(f_index){
  # inputs f index outputs predicted time spent in nest out of 600s
  out = ((f_index-2.268618)/-0.004926462)
  return(out)
}
saveRDS(time_formula, 'Functions/Theory/f_index_time.rds')


beh_formula = function(f_index){
   # inputs f index outputs predicted sum of parental care behaviors performed in 600s
  out <- exp((5.141225 - f_index) / 0.9294855)
  return(out)
}
saveRDS(beh_formula, 'Functions/Theory/f_index_behavior.rds')


test_formula = function(f_index){
  #inputs f index and outputs predicted % Testicular tissue in the gonad between 0 and 1
  out <- exp((-1.094307 - f_index) / 0.3785067)
  return(out)
}
saveRDS(test_formula, 'Functions/Theory/f_index_testicular.rds')

e2_formula = function(f_index){
  #inputs f index and outputs predicted Log10 transformed E2 in pg/ml
    out = (f_index- -3.687469 )/1.408077 
  return(out)
}
saveRDS(e2_formula, 'Functions/Theory/f_index_e2.rds')


kt_e2_formula = function(f_index){
    #inputs f index and outputs predicted Log10 transformed 11KT pg/ml/ log10 E2 pg/ml

    out = (f_index-2.386484 )/-2.434733 
  return(out)
}
saveRDS(kt_e2_formula, 'Functions/Theory/f_index_kt_e2.rds')

kt_formula = function(f_index){
    #inputs f index and outputs predicted Log10 transformed 11KT in pg/ml

  out = (f_index-4.561)/-1.946
  return(out)
}
saveRDS(kt_formula, 'Functions/Theory/f_index_11kt.rds')

mass_formula = function(f_index){
      #inputs f index and outputs mass in grams

  out = (f_index+2.285363) / 0.4422991
  return(out)
}
saveRDS(mass_formula, 'Functions/Theory/f_index_mass.rds')

mass_formula=readRDS( 'Functions/Theory/f_index_mass.rds')
kt_formula=readRDS( 'Functions/Theory/f_index_11kt.rds')
kt_e2_formula=readRDS( 'Functions/Theory/f_index_kt_e2.rds')
e2_formula=readRDS( 'Functions/Theory/f_index_e2.rds')
test_formula=readRDS( 'Functions/Theory/f_index_testicular.rds')
beh_formula=readRDS( 'Functions/Theory/f_index_behavior.rds')
time_formula=readRDS( 'Functions/Theory/f_index_time.rds')



### I talked to justin and he said neurons in dodd are pPOA (parvo)
# and the lsmeans is predict with the mean weight
library(readxl)
dodd_2019 = read_excel('Measures/dodd_2019_data.xlsx')
dodd_2019$time_better = ifelse(str_detect(dodd_2019$group2, '3.wks'), 21, NA )
dodd_2019$time_better = ifelse(str_detect(dodd_2019$group2, '6.mon'), 365/2, dodd_2019$time_better )
dodd_2019$time_better = ifelse(str_detect(dodd_2019$group2, '1.yr'), 365, dodd_2019$time_better )
dodd_2019$time_better = ifelse(str_detect(dodd_2019$group2, '3.yr'), 365*3, dodd_2019$time_better )
dodd_2019$time_better = ifelse(dodd_2019$group2 =='M', 0, dodd_2019$time_better )
dodd_2019$time_better = ifelse(dodd_2019$group2 =='F', 365*4, dodd_2019$time_better )

dodd_2019$phase = paste0(dodd_2019$Domint, dodd_2019$sex)
dodd_2019$phase2 = ifelse(dodd_2019$phase == 'Damb', 'D', NA)
dodd_2019$phase2 = ifelse(dodd_2019$phase == 'Samb', 'S', dodd_2019$phase2)
dodd_2019$phase2 = ifelse(dodd_2019$phase == 'SM', 'M', dodd_2019$phase2)
dodd_2019$phase2 = ifelse(dodd_2019$phase == 'DF', 'F', dodd_2019$phase2)
dodd_2019$Phase = dodd_2019$phase2

dodd_2019_good = subset(dodd_2019, phase2!='S')
dodd_2019_good$Log_11KT = log10(dodd_2019_good$KT11)
dodd_2019_good$Log_E2 = log10(dodd_2019_good$E2)
dodd_2019_good$mass_final = dodd_2019_good$`weight 2`



# lets write a function to produce f_index rq

f_indexer =function(data.frame){
  vars= c(
  'Log_11KT',
  'mass_final'
)

  newd = data.frame[complete.cases(data.frame[,vars]),]
newd = newd%>%
  subset(!Phase%in% c('S'))

X <- scale(newd[, vars])
pca <- prcomp(X, center = TRUE, scale. = TRUE)


ifelse(pca$rotation['mass_final','PC1']<0, newd$f_index <- pca$x[,1]*-1, newd$f_index <- pca$x[,1])

return(newd)
}
dodd_2019_good=f_indexer(dodd_2019_good)

ggplot(dodd_2019_good, aes(x = f_index, y = pPOAmedium))+
  geom_point(aes(color = Phase))

# regress out mass
mass_mod = lm(pPOAmedium~Phase+mass_final, data = dodd_2019_good)
anova(mass_mod, test ='Chisq')
dodd_2019_good$pred_pPOA = predict(mass_mod, newdata = data.frame(Phase = dodd_2019_good$Phase, 
                                                      mass_final =mean(dodd_2019_good$mass_final)))

ggplot(dodd_2019_good, aes(x = f_index, y = pred_pPOA))+
  geom_point(aes(color = Phase))

###

ggplot(dodd_2019_good, aes(x = f_index, y = pPOAarea))+
  geom_point(aes(color = Phase))

mass_mod2 = lm(pPOAarea~Phase+mass_final, data = dodd_2019_good)
anova(mass_mod2, test ='Chisq')
dodd_2019_good$pred_pPOA_area = predict(mass_mod2, newdata = data.frame(Phase = dodd_2019_good$Phase, 
                                                      mass_final =mean(dodd_2019_good$mass_final)))

ggplot(dodd_2019_good, aes(x = f_index, y = pred_pPOA_area))+
  geom_point(aes(color = Phase))

mass_mo = lm(pPOAarea~Phase, data = dodd_2019_good)
anova(mass_mo, test = 'Chisq')

##

parker_2023 = read.csv('Measures/parker_2023_data.csv')

coltan_dict = c('male'='M',
                'female' ='F',
                'dominant' ='D',
                'dom+eggs'='NF',
                'sub+' = 'S',
                'subordinate' = 'S')
parker_2023$Phase = coltan_dict[parker_2023$status2]
parker_2023$status2[is.na(parker_2023$Phase)]

parker_2023$mass_final = parker_2023$mass_sack_g
parker_2023$Log_11KT = log10(parker_2023$androgen_pg_mL)
parker_2023$Log_E2 = log10(parker_2023$estradiol_pg_mL)
parker_2023$days = parker_2023$week_sack_c*7
parker_2023$source = 'Parker'

ggplot(parker_2023, aes(x = Phase, y =poa_ant_brdu_avg))+
  geom_boxplot()+
  geom_point()

parker_2023 = f_indexer(parker_2023)

ggplot(parker_2023, aes(x = f_index, y =poa_ant_brdu_avg))+
  geom_point()

ggplot(parker_2023, aes(x = Phase, y =poa_ant_brdu_total))+
  geom_boxplot()+
  geom_point()

ggplot(parker_2023, aes(x = Phase, y =poa_ant_vol_mm3))+
  geom_boxplot()+
  geom_point()

ggplot(subset(parker_2023), aes(x = week_sack_c, y =poa_ant_brdu_total, color = Phase, group = interaction(Phase, week_sack_c)))+
  geom_boxplot()

ggplot(subset(parker_2023), aes(x = week_sack_c, y =f_index, color = Phase))+
  geom_point()+
  geom_smooth()

ggplot(parker_2023, aes(x = f_index, y =poa_mid_brdu_avg))+
  geom_point()

ggplot(parker_2023, aes(x = f_index, y =poa_post_brdu_avg))+
  geom_point()




#regress out volume
mass_mod_3 =lm(data = parker_2023, poa_ant_brdu_avg~Phase+poa_ant_vol_mm3)
anova(mass_mod_3, test = "Chisq")

parker_2023$pred_poa_ant_brdu_avg= predict(mass_mod_3,
                                   newdata = data.frame(Phase = parker_2023$Phase, 
                                                        poa_ant_vol_mm3 =mean(parker_2023$poa_ant_vol_mm3, na.rm =T)))

ggplot(parker_2023, aes(x = f_index, y =pred_poa_ant_brdu_avg))+
  geom_point(aes(color = Phase))

ggplot(parker_2023, aes(x = Phase, y =pred_poa_ant_brdu_avg))+
  geom_boxplot()

ggplot(parker_2023, aes(x = Phase, y =f_index))+
  geom_boxplot()

# can I git s vs D in my model??

ggplot(subset(parker_2023, !is.na(Phase)), aes(x = f_index, y =poa_ant_brdu_avg/poa_ant_vol_mm3))+
  geom_point(aes(color = Phase))

# i mean tbh I dont see anything here 

ggplot(subset(parker_2023, !is.na(Phase)), aes(x = f_index, y =poa_ant_vol_mm3/mass_final))+
  geom_point(aes(color = Phase))

#### lets look at logans data
ggplot(dodd_2019_good, aes(x = f_index, y = pPOAmedium/mass_final))+
  geom_point(aes(color =Phase))+
  geom_smooth(method = 'gam')

# this could be something
ggplot(dodd_2019_good, aes(x = f_index, y = pPOAmedium))+
  geom_point(aes(color =Phase, shape = 'Dodd', size =4))+
  geom_point(data = parker_2023, aes(x =f_index, y = poa_ant_brdu_total, color =Phase, shape = 'Parker', size =4))
# I think its tood different


# ill try my best

ggplot(subset(parker_2023, Phase %in% c("M","D","F")), aes(x = f_index, y =poa_full_brdu_total))+
  geom_point(aes(color =Phase))

attempt_gam_1 = gam(poa_full_brdu_total~s(f_index)+mass_final, data = parker_2023)
summary(attempt_gam_1)
plot(attempt_gam_1)

attempt_gam_1

attempt_glm_1 = glm(poa_full_brdu_total~(f_index)+mass_final, data = parker_2023)
summary(attempt_glm_1)

glm_equation =function(f_index#, mass_final
                       ){
 out= 643+(44.69*f_index)+(-19.04*mean(parker_2023$mass_final, na.rm =T))
return(out)
 }

parker_2023$fitted_glm = mapply(glm_equation, parker_2023$f_index#, parker_2023$mass_final
                                )

ggplot(subset(parker_2023, Phase %in% c("M","D","F")), aes(x = f_index, y =poa_full_brdu_total))+
  geom_point(aes(color =Phase))+
  geom_line(aes(x = f_index, y= fitted_glm))

# idk this feels like a waste 

nls_value_tester = function(a=1, b = 1, c = 1){

 p= ggplot(subset(parker_2023, Phase %in% c("M","D","F")), aes(x = f_index, y =poa_full_brdu_total))+
  geom_point(aes(color =Phase))+
  geom_line(aes(x = f_index, y= a*(f_index^2)+(b*f_index)+c))

 return(p)
  }

nls_value_tester()
nls_value_tester(a=-40, 50, 400+min(parker_2023$poa_ant_brdu_total))


nls_attempt = nls(poa_full_brdu_total~ a*(f_index^2)+(b*f_index)+c,
                  data = parker_2023,
                  start = list(
                    a = -40,
                    b = 50,
                    c=400+min(parker_2023$poa_ant_brdu_total)
                  ))

summary(nls_attempt)

nls_value_tester(a=-17.23, 33.62, 584.65)

anova(nls_attempt, attempt_glm_1)

# ok this is my best attempt 

poa_brdu_formula = function(f_index){
 
  out = ( -17.23*(f_index^2)+(33.62*f_index)+584.65)
return(out)
}

saveRDS(poa_brdu_formula, 'Functions/Theory/f_index_neurogenesis.rds')

line_data$neurogenesis = sapply(line_data$Z, poa_brdu_formula)

logged2 = ggplot(line_data, aes(x = Z, y = scale_2((kt))))+
  geom_line(aes(color = 'Log10 11KT (pg/ml)'))+
  geom_line( aes(linetype = 'Log10 E2 (pg/ml)', color = 'Log10 E2 (pg/ml)', y = scale_2((e2))))+
      geom_line(color = 'black',aes(linetype = 'Mass (g)',color = 'Mass (g)', y = scale_2((mass))))+
 # geom_line(aes(color = 'KT/ E2', y = scale(kt_e2)))+
  geom_line(data = subset(line_data, testicular <= 1),aes(color = 'Testicular %', y = scale_2(testicular)))+
  geom_line(data = subset(line_data, behavior <= 600),aes(color = 'Behaviors (10m)', y = scale_2(behavior)))+
  geom_line(data = subset(line_data, time <= 600),aes(color = 'Time in Nest (10m)', y = scale_2(time)))+
    geom_line(data = subset(line_data),aes(color = 'Neurogenesis', y = scale_2(neurogenesis)))+
  labs(y= 'Scaled Value', x = 'F Index', color = 'Measurement')+
  theme_minimal()
(logged2+boxes)+
  plot_layout(ncol = 1)




#### proportion of neurons from snRNA-seq data ####
obj = readRDS("~/Desktop/optimal_clustering_rna_only.rds")
obj@meta.data$Status = factor(obj@meta.data$Status, levels = c('NRM',
                                                               'M',
                                                               'D',
                                                               'E',
                                                               'NF',
                                                               'F'))

status_to_phase = list('D'='I',
                       'E' = 'LI',
                       'EP' = 'LIP',
                       'S' = 'IP',
                       'M'='M',
                       'F'='F',
                       'NF'='NF',
                       'NM'='NM')
multiome_measures$Phase = unlist(status_to_phase[multiome_measures$Status])

neur = obj@meta.data%>%
  group_by(individual, Status)%>%
  summarize(n_cells = n())

joint=multiome_measures%>%right_join(neur, by = join_by('Fish'=='individual'))%>%
  subset(!is.na(n_cells))
  
joint$Phase = factor(joint$Phase, levels = c('M','I','LI','NF','F'))

ggplot(joint, aes(x = Phase, y = n_cells))+
  geom_boxplot()+
  geom_point()+
  labs(y = 'Neurons')
  

obj@meta.data$neuron = ifelse(!obj@meta.data$final_clusters %in% c(1,2,11,13,20,15,26 ), 'neuron', 'not')

prop_neuron =obj@meta.data%>%
  group_by(individual, Status)%>%
  subset(individual != 'GH')%>%
  summarize(n_neuron = sum(neuron =='neuron'),
            prop_neuron = sum(neuron=='neuron')/n())

joint2 = joint%>%
  left_join(prop_neuron,by = join_by('Fish'=='individual'))

joint2$mass_final = joint2$mass_final_cm
joint2 = f_indexer(joint2)

ggplot(joint2, aes(x = f_index, y = prop_neuron))+
  geom_point(aes(color= Phase))+
  labs(y = 'Neurons')

# regress out mass
mod_mass_neuron= lm(prop_neuron~f_index+mass_final, data = joint2)
mod_mass_pred = predict(mod_mass_neuron, newdata = data.frame(f_index = joint2$f_index, 
                                                              mass_final = mean(joint2$mass_final, na.rm=T)))
summary(mod_mass_neuron)

joint2$pred_prop_neuron = mod_mass_pred

ggplot(joint2, aes(x = f_index, y = mod_mass_pred))+
  geom_point(aes(color= Phase, shape = 'fitted', size = 4))+
  geom_point(aes(color = Phase, y = prop_neuron, shape = 'Real', size =4))+
  labs(y = 'Neurons')

# not convinced
multiome_measures$mass_final = multiome_measures$mass_final_cm
multiome_measures= f_indexer(multiome_measures)
mod_mass_neuron2= glm(cbind(n_neuron, n_cells-n_neuron)~I(f_index^2), data = joint2, family = 'binomial')
summary(mod_mass_neuron2)

mod_mass_pred2 = predict(mod_mass_neuron2, newdata = data.frame(f_index = joint2$f_index, 
                                                              mass_final = mean(joint2$mass_final, na.rm=T))
                         )

test_dat = data.frame(f_index = seq(-1.5,3, by = 0.1)
                     # , mass_final = mean(multiome_measures$mass_final)
                      )

out_dat <- ggpredict(
  mod_mass_neuron2,
  terms = test_dat,
  scale = 'response'
)


ggplot(joint2, aes(x = f_index, y = pred_prop_neuron2))+
  geom_point(aes(color = Phase, y = prop_neuron, shape = 'Real'), size =4)+
  geom_line(data = out_dat, aes(x = x, y = predicted, color ='fitted'), color= 'black')+
  geom_ribbon(data = out_dat, aes(x = x, ymin = conf.low, ymax= conf.high), inherit.aes = F, alpha = 0.2)+
  labs(y = 'Neurons')

summary(mod_mass_neuron2)$dispersion
# doesnt make any sense
# chatgpt says to use quasibinomial with prop neuron but that is stupid I think
# there is nothing here

ggplot(joint2, aes(x = Phase, y = prop_neuron))+
  geom_boxplot(aes(color = Phase))+
  geom_point()
  ylim(0,1)












