library(ggplot2)
library(lavaan)
library(mgcv)
library(patchwork)

scale_2 <- function(x) {
  x <- as.numeric(x)
  rng <- range(x, na.rm = TRUE)
  (x - rng[1]) / (rng[2] - rng[1])
}

coalesced = read.csv('Measures/coalesced_data.csv')
measures = read.csv('Measures/all_data.csv')
measures$mass_final = measures$mass_final_cm
status_to_phase = list('D'='I',
                       'E' = 'LI',
                       'EP' = 'LIP',
                       'S' = 'IP',
                       'M'='M',
                       'F'='F',
                       'NF'='NF',
                       'NM'='NM')
measures$Phase = unlist(status_to_phase[measures$Status])
measures$mass_final = measures$mass_final_cm
measures$Percent_Testicular[measures$Phase %in% c('F','NF')] =0 
# true by definition

measurements = c('Log_11KT', 'mass_final', 'Percent_Testicular') 
#these are much more objective than behaviors so any theory based on them will
#> much more consistent than any model based on behavior, we can see how behav-
#> ior relates to the latent f after the fact, but should not use it in the mo-
#> del I think
#> 
#> We can also test this model on other data to see how well it correalates
#> with priors, if dodd or parker data is correctly defined as male, female, I
#> etc, then we are good

# vars must be scaled
measures_of_interest = subset(measures, Phase %in% c('M','I','LI','NF','F'))
#latent variable model

sex_model <- '
  # latent variable
  SexState =~ 
    1*Log_11KT +
    Percent_Testicular +
    mass_final
    
    SexState ~~ SexState

'
fit <- lavaan::sem(
  sex_model,
  data = measures_of_interest,
  missing = "fiml"
)

summary(fit, standardized = TRUE)
inspect(fit, "coef")

measures_of_interest$Phase = factor(measures_of_interest$Phase, levels = c('M','I','LI','NF','F'))

measures_of_interest$SexScore = lavaan::lavPredict(fit)[,1]*-1 # i like female being higher

# test if loadings are all large and same sized
summary(fit, standardized = TRUE)

#standardized loadings
std <- standardizedSolution(fit)
std[std$op == "=~", ]

# residual variances
std[std$op == "~~" & std$lhs == std$rhs, ]

#measures of fit
fitMeasures(fit, c("cfi", "tli", "rmsea", "srmr"))

# mod indicies
modindices(fit)

#r2
inspect(fit, "r2")

#measurement invariance across phase
coalesced_good = subset(coalesced, Phase %in% c('M','D','E', 'F'))
table(coalesced_good$Phase)
coalesced_good$Percent_Testicular[coalesced_good$Phase == 'F']

#configural invariance, test with larger dataset
fit_configural <- lavaan::sem(
  sex_model,
  data = coalesced_good,
  group = "Phase",
  missing = "fiml"
)

#metric invariance
fit_metric <- lavaan::sem(
  sex_model,
  data = measures_of_interest,
  group = "Phase",
  group.equal = c("loadings"),
  missing = "fiml"
)

#scalar invariance
fit_scalar <- lavaan::sem(
  sex_model,
  data = measures_of_interest,
  group = "Phase",
  group.equal = c("loadings", "intercepts"),
  missing = "fiml"
)

# I think this cannot work because percent testicular does not vary in females

ggplot(measures_of_interest, aes(x =Phase, y = SexScore))+
  geom_boxplot()+
    geom_point(aes(color = Phase))

ggplot(subset(measures_of_interest, !Phase %in% c('NM','LIP','IP')), aes(y =Percent_Testicular, x = SexScore))+
  geom_point(aes(color = Phase))

ggplot(subset(measures_of_interest, !Phase %in% c('NM','LIP','IP')), aes(y =mass_final, x = SexScore))+
  geom_point(aes(color = Phase))

ggplot(subset(measures_of_interest, !Phase %in% c('NM','LIP','IP')), aes(y =Log_11KT, x = SexScore))+
  geom_point(aes(color = Phase))


ggplot(subset(measures_of_interest, !Phase %in% c('NM','LIP','IP')), aes(y =(Time_Day_2), x = SexScore))+
  geom_point(aes(color = Phase))

ggplot(subset(measures_of_interest, !Phase %in% c('NM','LIP','IP')), aes(y =(Behaviors_Day_2), x = SexScore))+
  geom_point(aes(color = Phase))

### Describe trajectories ####
#testicular
eq_function_test = function(formula_string){
  
  dat = data.frame(SexScore = seq(min(measures_of_interest$SexScore), max(measures_of_interest$SexScore), by =0.1))
  
  dat$line = with(dat, eval(parse(text = formula_string)))
  
  p = ggplot(subset(measures_of_interest, !Phase %in% c('NM','LIP','IP')), 
             aes(y = Percent_Testicular, x = SexScore )) +
    geom_point(aes(color = Phase)) +
    geom_line(data = dat, aes(x = SexScore, y = line, color = formula_string))
  
  return(p)  
}

eq_function_test('(0.5*(SexScore)^(2))-(.4*SexScore)+0.07')

sex_testic = nls(Percent_Testicular ~ (a*(SexScore)^2)-(b*SexScore)+c,
                 start = list(
                   a = 0.5, 
                   b = .4,
                   c = 0.07
                 ),
                 data = measures_of_interest)
summary(sex_testic)

eq_function_test('(0.4012*(SexScore)^(2))-(.52565*SexScore)+0.16047')

#with a logarithmic decay
eq_function_test('(-.4*log10(SexScore+0.6)+0.1)')
#y = a logb(x - h) + k
eq_function_test('1/(1.5+exp(5)^(SexScore+0.3))')


sex_testic_2 = nls(Percent_Testicular ~ 1/(a+ exp(c*(SexScore + d))),
                 start = 
                   list(
                     a = 1.5,
                     c = 5,
                     d = 0.3
                     ),
                 data = measures_of_interest)
summary(sex_testic_2)

eq_function_test('1/(1.1663+ exp(4.6027*(SexScore + 0.3996)))')

anova(sex_testic, sex_testic_2) # second is barely better


### 11KT 
eq_function_kt = function(formula_string){
  
  dat = data.frame(SexScore = seq(min(measures_of_interest$SexScore), max(measures_of_interest$SexScore), by =0.1))
  
  dat$line = with(dat, eval(parse(text = formula_string)))
  
  p = ggplot(subset(measures_of_interest, !Phase %in% c('NM','LIP','IP')), 
             aes(y = Log_11KT, x = SexScore )) +
    geom_point(aes(color = Phase)) +
    geom_line(data = dat, aes(x = SexScore, y = line, color = formula_string))
  
  return(p)  
}
eq_function_kt('(SexScore*-1)+2')

#get from linear regression
kt_func = lm(Log_11KT ~ SexScore, data = measures_of_interest)
summary(kt_func)

eq_function_kt('(SexScore*-1.11394)+2.32032')
#im not super happy with how Is are systematically underpredicted

eq_function_kt('(-1*(-1*SexScore)^(2))-(SexScore)+2.7')

kt_func2 = nls(Log_11KT ~ (a*(-1*SexScore)^2)-(b*SexScore)+c,
               data = measures_of_interest,
               start = list(
                 a = -1, 
                 b = 1,
                 c = 2.7
               ))
summary(kt_func2)

eq_function_kt('(-0.7413*(-1*SexScore)^(2))-(1.1285*SexScore)+2.4970')
t = summary(kt_func) 
sum(sqrt(t$residuals^2))
t2 = summary(kt_func2)
sum(sqrt(t2$residuals^2))
# second is better
### mass
eq_function_mass = function(formula_string){
  
  dat = data.frame(SexScore = seq(min(measures_of_interest$SexScore), max(measures_of_interest$SexScore), by =0.1))
  
  dat$line = with(dat, eval(parse(text = formula_string)))
  
  p = ggplot(subset(measures_of_interest, !Phase %in% c('NM','LIP','IP')), 
             aes(y = mass_final, x = SexScore )) +
    geom_point(aes(color = Phase)) +
    geom_line(data = dat, aes(x = SexScore, y = line, color = formula_string))
  
  return(p)  
}
eq_function_mass('(SexScore*5)+6')

mass_func = lm(mass_final~SexScore, data = measures_of_interest)
summary(mass_func)
eq_function_mass('(SexScore*5.4986)+6.2182')


### behavior, this one is gonna be tough
eq_function_beh = function(formula_string){
  
  dat = data.frame(SexScore = seq(min(measures_of_interest$SexScore), max(measures_of_interest$SexScore), by =0.1))
  
  dat$line = with(dat, eval(parse(text = formula_string)))
  
  p = ggplot(subset(measures_of_interest, !Phase %in% c('NM','LIP','IP')), 
             aes(y = Behaviors_Day_2, x = SexScore )) +
    geom_point(aes(color = Phase)) +
    geom_line(data = dat, aes(x = SexScore, y = line, color = formula_string))
  
  return(p)  
}
eq_function_beh('(SexScore*-400)+300')

# it looks linear to me though
beh_mod = lm(Behaviors_Day_2 ~ SexScore, data = measures_of_interest)
summary(beh_mod)

eq_function_beh('(SexScore*-252.32)+276.09')

# time
eq_function_time = function(formula_string){
  
  dat = data.frame(SexScore = seq(min(measures_of_interest$SexScore), max(measures_of_interest$SexScore), by =0.1))
  
  dat$line = with(dat, eval(parse(text = formula_string)))
  
  p = ggplot(subset(measures_of_interest, !Phase %in% c('NM','LIP','IP')), 
             aes(y = Time_Day_2, x = SexScore )) +
    geom_point(aes(color = Phase)) +
    geom_line(data = dat, aes(x = SexScore, y = line, color = formula_string))
  
  return(p)  
}
eq_function_time('(-100*(-1*SexScore)^(2))-(200*SexScore)+500')
# looks more parabolic to me?

tim_model = nls(data = measures_of_interest, 
                formula = Time_Day_2 ~ (a*(-1*SexScore)^2)-(b*SexScore)+c,
                start = list(
                  a =-100,
                  b = 200,
                  c = 500
                ))
summary(tim_model)

eq_function_time('(-127.86*(-1*SexScore)^(2))-(194.19*SexScore)+435.27')

tim_model2 = lm(Time_Day_2~SexScore, data = measures_of_interest)

anova(tim_model, tim_model2)

#e2
coalesced$SexScore = predict(fit, coalesced)*-1
# looks linear 
e2_mod = lm(Log_E2~I(SexScore*-1), data = coalesced)
summary(e2_mod)

e2_func = function(SexScore){
  out = (SexScore*1.23989)+ 2.64937
  return(out)
}

### functions
test_func = function(SexScore){
  out = 1/(1.1663+ exp(4.6027*(SexScore + 0.3996)))
return(out)
}

kt_func = function(SexScore){
  out = (-0.7413*(-1*SexScore)^(2))-(1.1285*SexScore)+2.4970
  return(out)
}

mass_func= function(SexScore){
  out = ((SexScore*5.4986)+6.2182)
  return(out)
}

beh_func = function(SexScore){
  out = ((SexScore*-252.32)+276.09)
return(out)
}


time_func = function(SexScore){
  out = (-127.86*(-1*SexScore)^(2))-(194.19*SexScore)+435.27
  return(out)
}

e2_func = function(SexScore){
  out = (SexScore*1.23989)+ 2.64937
  return(out)
}

test_dat = data.frame(SexScore = seq(-1, 1, by = 0.1))
test_dat$Percent_Testicular = sapply(test_dat$SexScore, 
                                     test_func)
test_dat$Log_11KT = sapply(test_dat$SexScore, 
                                     kt_func)
test_dat$mass_final = sapply(test_dat$SexScore, 
                                     mass_func)
test_dat$Behaviors_Day_2 = sapply(test_dat$SexScore, 
                                     beh_func)
test_dat$Time_Day_2 = sapply(test_dat$SexScore, 
                                     time_func)

ggplot(test_dat, aes(x = SexScore, y = scale(Log_11KT)))+
  geom_line(aes(color = 'Log_11KT'))+
  geom_line(aes(color = 'Percent Testicular', y = scale(Percent_Testicular)))+
  geom_line(aes(color = 'Mass', y = scale(mass_final)))+
  geom_line(aes(color = 'Behavior', y = scale(Behaviors_Day_2)))+
  geom_line(aes(color = 'Time', y = scale(Time_Day_2)))

ggplot(test_dat, aes(x = SexScore, y = scale_2(Log_11KT)))+
  geom_line(aes(color = 'Log_11KT'))+
  geom_line(aes(color = 'Percent Testicular', y = scale_2(Percent_Testicular)))+
  geom_line(aes(color = 'Mass', y = scale_2(mass_final)))+
  geom_line(aes(color = 'Behavior', y = scale_2(Behaviors_Day_2)))+
  geom_line(aes(color = 'Time', y = scale_2(Time_Day_2)))

# read in other data
coalesced$SexScore = predict(fit, coalesced)*-1

ggplot(coalesced, aes(x= SexScore, y =Percent_Testicular))+
  geom_point(aes(color = Phase))

ggplot(subset(coalesced, source!='DeAngelis 2016' & Phase %in% c('M','D','E','NF','F')), aes(x= SexScore, y =Log_11KT))+
  geom_point(aes(color = Phase))+
  geom_line(data = test_dat, aes(x = SexScore, y= Log_11KT))

ggplot(subset(coalesced, source=='Graham' & Phase %in% c('M','D','E','NF','F')), aes(x= SexScore, y =Log_11KT))+
  geom_point(aes(color = Phase))+
  geom_line(data = test_dat, aes(x = SexScore, y= Log_11KT))

# clearly my data is a little special and I need improvement, this is better modeled by a sigmoid
eq_function_kt2 = function(formula_string){
  
  dat = data.frame(SexScore = seq(min(coalesced$SexScore, na.rm = T), max(coalesced$SexScore,na.rm = T), by =0.1))
  
  dat$line = with(dat, eval(parse(text = formula_string)))
  
  p = ggplot(subset(coalesced,source != 'DeAngelis 2016'& Phase %in% c('M','D','E','NF','F')), 
             aes(y = Log_11KT, x = SexScore )) +
    geom_point(aes(color = Phase)) +
    geom_line(data = dat, aes(x = SexScore, y = line, color = formula_string))
  
  return(p)  
}
eq_function_kt2('1/( exp(1*(SexScore-0.6)))')

kt_3 =nls(Log_11KT ~ 1/(a+ exp(c*(SexScore + d))),
                 start = 
                   list(
                     a = 0,
                     c = 1,
                     d = -0.6
                     ),
                 data = subset(coalesced, source!= 'DeAngelis 2016'&
                                  Phase %in% c('M','D','E','NF','F')))
summary(kt_3)

eq_function_kt2('1/( 0.1392+exp(0.7603*(SexScore-1.5265 )))')

kt_func3 = function(SexScore){
  out = 1/( 0.1392+exp(0.7603*(SexScore-1.5265 )))
  return(out)
}

ggplot(subset(coalesced, Phase != 'S'), aes(x= SexScore, y =days))+
  geom_point(aes(color = Phase))
# uh oh problem

ggplot(subset(coalesced, Phase != 'S'), aes(x= SexScore, y =Log_E2))+
  geom_point(aes(color = Phase))
# not a super great correlation 

test_dat$Log_E2 = sapply(test_dat$SexScore, e2_func)
test_dat$Percent_Testicular2 = sapply(test_dat$SexScore, test_func)


eq_function_kt('1/( 0.1392+exp(0.7603*(SexScore-1.5265 )))') # fit with other data
eq_function_test('1/(1.1663+ exp(4.6027*(SexScore + 0.3996)))')
eq_function_mass('(SexScore*5.4986)+6.2182')
eq_function_beh('(SexScore*-252.32)+276.09') # i mean this is rly pushin it lol
eq_function_time('(-127.86*(-1*SexScore)^(2))-(194.19*SexScore)+435.27')

test_dat$Log_11KT = sapply(test_dat$SexScore, 
                                     kt_func3)


plot1=ggplot(test_dat, aes(x = SexScore, y = scale_2(Log_11KT)))+
  geom_line(aes(color = 'Log_11KT'))+
  geom_line(aes(color = 'Percent Testicular', y = scale_2(Percent_Testicular)))+
  geom_line(aes(color = 'Mass', y = scale_2(mass_final)))+
  geom_line(aes(color = 'Behavior', y = scale_2(Behaviors_Day_2)))+
  geom_line(aes(color = 'Time', y = scale_2(Time_Day_2)))+
    geom_line(aes(color = 'Log_E2', y = scale_2(Log_E2)+0.05))

plot2 =ggplot(measures_of_interest, aes(x = SexScore, y = Phase))+
  geom_boxplot()+
  geom_point()

(plot1+plot2)+
  plot_layout(ncol = 1)

#save progress
saveRDS(test_func, 'Functions/Theory/latent_test.rds')
saveRDS(e2_func, 'Functions/Theory/latent_e2.rds')
saveRDS(mass_func,'Functions/Theory/latent_mass.rds' )
saveRDS(beh_func,'Functions/Theory/latent_beh.rds' )
saveRDS(time_func,'Functions/Theory/latent_time.rds' )
saveRDS(kt_func3, 'Functions/Theory/latent_kt.rds' )

saveRDS(fit, 'Functions/Theory/latent_fit.rds')
saveRDS(sex_model, 'Functions/Theory/latent_sex_model.rds')


# read back in 
test_func=readRDS('Functions/Theory/latent_test.rds')
e2_func=readRDS( 'Functions/Theory/latent_e2.rds')
mass_func=readRDS('Functions/Theory/latent_mass.rds' )
beh_func=readRDS('Functions/Theory/latent_beh.rds' )
time_func=readRDS('Functions/Theory/latent_time.rds' )
kt_func =readRDS('Functions/Theory/latent_kt.rds' )

write.csv(measures_of_interest, 'Measures/measures_SexScore.csv')

####
coalesced_good = subset(coalesced, Phase %in% c('M',
                                                'D',
                                                'E',
                                                'NF',
                                                'F'))

coalesced_pred = predict(fit, newdata=coalesced_good)
coalesced_good$SexScore =coalesced_pred*-1

ggplot(coalesced_good, aes(x = SexScore, y = Percent_Testicular))+
  geom_point()+
  geom_line(data = test_dat, aes(color = 'Percent Testicular', y = (Percent_Testicular)))

ggplot(coalesced_good, aes(x = SexScore, y = Log_11KT))+
  geom_point()+
  geom_line(data = test_dat, aes(color = 'KT', y = (Log_11KT)))

ggplot(coalesced_good, aes(x = SexScore, y = mass_final))+
  geom_point(aes(colour = Phase))+
  geom_line(data = test_dat, aes(color = 'Mass', y = (mass_final)))


### looking at proportions of celltypes ####
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
measures_of_interest$Phase = unlist(status_to_phase[measures_of_interest$Status])

neur = obj@meta.data%>%
  group_by(individual, Status)%>%
  summarize(n_cells = n())

clust = obj@meta.data%>%
  group_by(individual, Status, final_clusters)%>%
  summarize(i_cells = n())

joint = neur%>%
  left_join(clust, by = 'individual')%>%
  mutate(prop_cells = i_cells/n_cells)%>%
  left_join(measures_of_interest, by = join_by('individual'=='Fish'))

ggplot(subset(joint,final_clusters == 1 ), aes(x = SexScore, y = prop_cells, color = Phase))+
  geom_point()

ggplot(subset(joint,final_clusters ==6 ), aes(x = SexScore, y = prop_cells, color = Phase))+
  geom_point()

ggplot(subset(joint,final_clusters ==2 ), aes(x = SexScore, y = prop_cells, color = Phase))+
  geom_point()

ggplot(subset(joint,final_clusters ==10), aes(x = SexScore, y = prop_cells, color = Phase))+
  geom_point()

for(i in 0:26){
  print(i)
  p = ggplot(subset(joint,final_clusters ==i), aes(x = SexScore, y = prop_cells, color = Phase))+
  geom_point()
  print(p)
  Sys.sleep(5)
}
#3,6,11,24,26

ggplot(subset(joint,final_clusters ==3), aes(x = SexScore, y = (prop_cells), color = Phase))+
  geom_point()
ggplot(subset(joint,final_clusters ==6), aes(x = SexScore, y = prop_cells, color = Phase))+
  geom_point()
ggplot(subset(joint,final_clusters ==11), aes(x = SexScore, y = prop_cells, color = Phase))+
  geom_point()
ggplot(subset(joint,final_clusters ==24), aes(x = SexScore, y = prop_cells, color = Phase))+
  geom_point()
ggplot(subset(joint,final_clusters ==26), aes(x = SexScore, y = (prop_cells), color = Phase))+
  geom_point()

####
coalesced$SexScore= predict(fit, coalesced)*-1

ggplot(coalesced, aes(SexScore,days ))+
  geom_point()

days_mod = lm(days~SexScore, data = coalesced)
summary(days_mod)

# i think a correlation matrix would be cool here but Im too lazy to do that RN
# i guess lets just test

out = data.frame()
for(i in 0:26){
mod = lm(prop_cells~SexScore, data = subset(joint, final_clusters ==i))
m = summary(mod)
p.val = m$coefficients[2,4]

newd = data.frame(cluster =i,
                  p.val = p.val)
out = rbind(out, newd)
}
out$q.val = p.adjust(out$p.val, 'fdr', 27)
# nuffink






