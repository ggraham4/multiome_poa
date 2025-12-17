library(forcats)
library(tidyverse)
library(ggplot2)
library(lme4)
library(readxl)
library(patchwork)
library(mgcv)

#all the clownfish data

dodd_2019 = read_excel('Measures/dodd_2019_data.xlsx')

dodd_2019$time_better = ifelse(str_detect(dodd_2019$group2, '3.wks'), 21, NA )
dodd_2019$time_better = ifelse(str_detect(dodd_2019$group2, '6.mon'), 365/2, dodd_2019$time_better )
dodd_2019$time_better = ifelse(str_detect(dodd_2019$group2, '1.yr'), 365, dodd_2019$time_better )
dodd_2019$time_better = ifelse(str_detect(dodd_2019$group2, '3.yr'), 365*3, dodd_2019$time_better )
dodd_2019$time_better = ifelse(dodd_2019$group2 =='M', 0, dodd_2019$time_better )

dodd_2019$phase = paste0(dodd_2019$Domint, dodd_2019$sex)
dodd_2019$phase2 = ifelse(dodd_2019$phase == 'Damb', 'D', NA)
dodd_2019$phase2 = ifelse(dodd_2019$phase == 'Samb', 'S', dodd_2019$phase2)
dodd_2019$phase2 = ifelse(dodd_2019$phase == 'SM', 'M', dodd_2019$phase2)
dodd_2019$phase2 = ifelse(dodd_2019$phase == 'DF', 'F', dodd_2019$phase2)

dodd_2019$days = dodd_2019$time_better

dodd_2019$Phase = dodd_2019$phase2
dodd_2019$phase2[is.na(dodd_2019$Phase)]

dodd_2019$Log_11KT= log10(dodd_2019$KT11)
dodd_2019$Log_E2= log10(dodd_2019$E2)
dodd_2019$mass_final= (dodd_2019$`weight 2`)
dodd_2019$source ='Dodd'
dodd_2019$Percent_Testicular = dodd_2019$pertest

parker_2023 = read.csv('Measures/parker_2023_data.csv')

coltan_dict = c('male'='M',
                'female' ='F',
                'dominant' ='D',
                'dom+eggs'='NF',
                'sub+' = 'S',
                'subordinate' = 'S')
parker_2023$Phase = coltan_dict[parker_2023$status2]
parker_2023$status2[is.na(parker_2023$Phase)]
parker_2023$Percent_Testicular = NA

parker_2023$mass_final = parker_2023$mass_sack_g
parker_2023$Log_11KT = log10(parker_2023$androgen_pg_mL)
parker_2023$Log_E2 = log10(parker_2023$estradiol_pg_mL)
parker_2023$days = parker_2023$week_sack_c*7
parker_2023$source = 'Parker'


multiome_measures = read.csv('Measures/all_data.csv')
multiome_measures$Time_Day_2= as.numeric(multiome_measures$Time_Day_2)
multiome_measures$Behaviors_Day_2= as.numeric(multiome_measures$Behaviors_Day_2)
multiome_measures$Log_11KT = as.numeric(multiome_measures$Log_11KT)

multiome_measures$Phase = multiome_measures$Status

graham_2026 = multiome_measures
graham_2026$mass_final =graham_2026$mass_final
graham_2026$source ='Graham'

graham_2026$days = ifelse(graham_2026$Phase %in% c('E','I','NF','NM','IP','S'), 180,NA)
graham_2026$days = ifelse(graham_2026$Phase %in% c('M'), 0, graham_2026$days)


### add in gonzalez data
gonzalez_2021 = read_excel("Measures/gonzalez_2021_data.xlsx")
gonzalez_2021$Percent_Testicular = gonzalez_2021$pertest
gonzalez_2021$mass_final = gonzalez_2021$mass6
gonzalez_2021$Log_11KT = log10(gonzalez_2021$`11kt6`%>%as.numeric())
gonzalez_2021$Log_E2 = log10(gonzalez_2021$E26%>%as.numeric())
gonzalez_2021$days= 180
gonzalez_2021$Phase = ifelse(gonzalez_2021$status == 'dom', 'D','S')
gonzalez_2021 = subset(gonzalez_2021, treat == 'Control')
gonzalez_2021$source = 'Gonzalez'

### deangelis 2018
deangelis_2018 = read.csv("Measures/deangelis_2018_data.csv")

deangelis_2018$Log_E2 = log10(deangelis_2018$E22)
deangelis_2018$Log_11KT = log10(deangelis_2018$X11kt2)
deangelis_2018$mass_final = deangelis_2018$Weight
deangelis_2018$source = 'DeAngelis 2018'
deangelis_2018=subset(deangelis_2018, Sex %in% c('f', 'm'))
deangelis_2018$Phase = ifelse(deangelis_2018$Sex=='f', 'F','M' )
deangelis_2018$days = ifelse(deangelis_2018$Phase=='M', 0,NA )
deangelis_2018$Percent_Testicular = NA


# deangelis 2016
deangelis_2016 = read.csv("Measures/deangelis_2016_data.csv")
deangelis_2016$Log_11KT = log10(deangelis_2016$X11KT%>%as.numeric())
deangelis_2016$Log_E2 = log10(deangelis_2016$E2.Level%>%as.numeric())
deangelis_2016$mass_final = deangelis_2016$Weight
deangelis_2016$source = 'DeAngelis 2016'
deangelis_2016$Phase = ifelse(deangelis_2016$MF=='f', 'F','M')
deangelis_2016$days = ifelse(deangelis_2016$Phase=='M', 0,NA )
deangelis_2016$Percent_Testicular = NA

common_cols = c('source',
                'Phase',
                'Log_11KT',
                'Log_E2',
                'days',
                'mass_final',
                'Percent_Testicular')

coalesced = rbind(deangelis_2018[,common_cols],
                  dodd_2019[,common_cols],
                  parker_2023[,common_cols],
                  graham_2026[,common_cols],
                  gonzalez_2021[,common_cols],
                  deangelis_2016[,common_cols])

coalesced$mass_final = as.numeric(coalesced$mass_final)


ggplot(subset(coalesced, Phase != 'S'), aes(x = days, y = mass_final))+
  geom_point(position = position_dodge(1),aes( shape = Phase, color = source), size =3)+
  xlim(-1, 380)+
  geom_smooth()

coalesced$kt_e2 = coalesced$Log_11KT/coalesced$Log_E2

ggplot(subset(coalesced, Phase != 'S'& Phase!= 'EP'), aes(x = kt_e2, y = mass_final))+
  geom_point(position = position_dodge(1),aes( shape = Phase, color = source), size = 3)+
  geom_smooth()

p1 = ggplot(subset(coalesced, Phase != 'S'), aes(x = days, y = kt_e2))+
  geom_point(position = position_dodge(1),aes( shape = Phase, color = source), size = 3)+
  geom_smooth()+
      ylim(0, 2)+
  xlim(-1, 400)


p2=ggplot(subset(coalesced, Phase == 'F'), aes(x='F', y = kt_e2))+
  geom_point(aes( shape = Phase, color = source), shape=3, size = 4)+
      ylim(0, 2)

p1+p2

ggplot(subset(coalesced, Phase != 'S'& Phase!= 'EP'), aes(x = Log_E2, y = mass_final))+
  geom_point(position = position_dodge(1),aes( shape = Phase, color = source), size = 4)+
  geom_smooth()

ggplot(subset(coalesced, Phase != 'S'& Phase!= 'EP'), aes(x = Log_11KT, y = mass_final))+
  geom_point(position = position_dodge(1),aes( shape = Phase, color = source), size = 3)+
  geom_smooth()

ggplot(subset(coalesced, Phase != 'S' & Phase!= 'EP'), aes(x = Log_11KT, y = Percent_Testicular))+
  geom_point(position = position_dodge(1),aes( shape = Phase, color = source), size = 3)+
  geom_smooth()

ggplot(subset(coalesced, Phase != 'S'), aes(x = Log_E2, y = Percent_Testicular))+
  geom_point(position = position_dodge(1),aes( shape = Phase, color = source), size = 3)


#pca analysis
coalesced$Phase[is.na(coalesced$Phase)] = 'IDK'

pca_1 = coalesced[,c(2,3,6,7)]%>%na.omit()%>%
  subset(!Phase %in% c('IDK',
                       'S',
                       'NM',
                       'EP'))

pca_kt_mass_test = prcomp(pca_1[,c(2,3,4)])
library(factoextra)

fviz_pca_biplot(pca_kt_mass_test, habillage = pca_1$Phase)

pca_2 = coalesced[,c(2,3,6)]%>%na.omit()%>%
  subset(!Phase %in% c('IDK',
                       'S',
                       'NM',
                       'EP'))

pca_kt_mass = prcomp(pca_2[,c(2,3)])

fviz_pca_biplot(pca_kt_mass, habillage = pca_2$Phase)

ggplot(subset(coalesced, Phase != 'S'& Phase!= 'EP'), aes(x = Log_11KT, y = mass_final))+
  geom_point(position = position_dodge(1),aes( shape = Phase, color = source), size = 3)+
  geom_smooth()

# model time
gam_progress <- gam(
  Phase %in% c('F','D','E') ~ s(Log_11KT) +s(Log_E2)+ s(mass_final),
  family = binomial,
  data = subset(coalesced, Phase %in% c('M','D','F'))
)
plot(gam_progress)

coalesced$pfemale <- predict(
  gam_progress,
  newdata = coalesced,
  type = "response"
)

ggplot(coalesced, aes(x = Phase, y = pfemale))+
  geom_point()+
  geom_boxplot(alpha = 0)

gam_progress_noE2 <- gam(
  Phase %in% c('F','D') ~ 
    s(Log_11KT) + s(mass_final),
  family = binomial,
  data = subset(coalesced, Phase %in% c('M','D','F') & !is.na(Log_11KT))
)
coalesced$pfemale_noE2 <- predict(
  gam_progress_noE2,
  newdata = coalesced,
  type = "response"
)
gam_E2 <- gam(
  Log_E2 ~ s(pfemale_noE2) + s(mass_final),
  data = subset(coalesced, !is.na(Log_E2)),
  method = "REML"
)
coalesced$Log_E2_imputed <- coalesced$Log_E2

coalesced$Log_E2_imputed <- predict(
  gam_E2,
  newdata = coalesced
)

ggplot(subset(coalesced, Phase %in% c('F','M','D','E')), aes(x = Phase, y = Log_E2))+
  geom_boxplot(alpha = 0, aes(y = Log_E2, color ='true'),position = position_nudge(x = -0.21), width = 0.4)+
    geom_boxplot(alpha = 0, aes(y = Log_E2_imputed, color ='imputed'),position = position_nudge(x = 0.21), width = 0.4)

### train test split 
set.seed(123)

# only rows that could ever be used for E2 evaluation
coalesced$female_regime <- coalesced$Phase == 'F'

idx <- which(!is.na(coalesced$Log_11KT) & 
             coalesced$Phase %in% c('M','D','F','E'))

train_idx <- sample(idx, size = floor(0.7 * length(idx)))
test_idx  <- setdiff(idx, train_idx)

train <- coalesced[train_idx, ]
test  <- coalesced[test_idx, ]


train$Phase <- factor(train$Phase, levels = c('M','D','F'))
train$female_regime <- as.numeric(train$Phase == 'F')
test$female_regime  <- as.numeric(test$Phase == 'F')

gam_E2 <- gam(
  Log_E2 ~ 
    s(Log_11KT) +                 
    s(mass_final),
  data = subset(train, !is.na(Log_E2)),
  method = "REML"
)

train$Log_E2_imputed <- predict(gam_E2, newdata = train)
test$Log_E2_imputed  <- predict(gam_E2, newdata = test)

eval_test <- subset(test, !is.na(Log_E2))

with(eval_test, cor(Log_E2, Log_E2_imputed, use = 'pairwise.complete.obs'))
with(subset(eval_test, !is.na(Log_E2) & !is.na(Log_E2_imputed)), mean((Log_E2 - Log_E2_imputed)^2))
#MSE = 0.2 and Correlation = 0.9

ggplot(eval_test, aes(Log_E2, Log_E2_imputed, shape = Phase, color = Phase)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  coord_equal()

eval_test$resid <- eval_test$Log_E2 - eval_test$Log_E2_imputed

eval_test%>%group_by(Phase)%>%
  subset(!is.na(Log_E2) &
           !is.na(Log_E2_imputed)&
           !is.na(resid))%>%
  summarize(mean_resid = mean(resid))

ggplot(eval_test, aes(pfemale_noE2, resid, color = Phase)) +
  geom_hline(yintercept = 0) +
  geom_point() 

#predict E & D E2
graham_2026$Phase = factor(graham_2026$Phase, levels = c('S','EP','NM',
                                                   'M',
                                                   "D",
                                                   'E',
                                                   'NF',
                                                   "F"))

graham_2026$Log_E2_imputed  <- predict(gam_E2, newdata = graham_2026)

ggplot(graham_2026, aes(Phase, Log_E2_imputed)) +
  geom_boxplot()+
    geom_point() 


graham_2026$kt_e2=graham_2026$Log_11KT/graham_2026$Log_E2_imputed

ggplot(graham_2026, aes(Phase, kt_e2)) +
  geom_boxplot()+
    geom_point() 

graham_2026$Log_E2 = graham_2026$Log_E2_imputed

imputation_babey = rbind(coalesced[,c(common_cols,'kt_e2')], graham_2026[,c(common_cols,'kt_e2')])
imputation_babey$imputed = ifelse(imputation_babey$source=='Graham', 'imputed', 'true')
imputation_babey$imputed = factor(imputation_babey$imputed, levels =c('true',
                                                                      'imputed'))
imputation_babey$Phase = factor(imputation_babey$Phase, levels = c('S',
                                                             'EP',
                                                             'NM',
                                                             'M',
                                                             "D",
                                                             "E",
                                                             'NF',
                                                             'F'))
ggplot(subset(imputation_babey, !is.na(Phase)), aes(Phase, kt_e2, color = imputed)) +
  geom_boxplot()+
    geom_point(position = position_dodge(0.75)) +
  labs(y = 'Log10 11KT/Log10 E2', color = '')

ggplot(subset(imputation_babey, !is.na(Phase)), aes(Phase, Log_E2, color = imputed)) +
  geom_boxplot()+
    geom_point(position = position_dodge(0.75)) +
  labs(y = 'Log10 E2', color = '')

write.csv(coalesced, 'Measures/coalesced_data.csv')


#theoretical model of sex change
gam_progress_noE2 <- gam(
  Phase %in% c('F','D') ~ 
    s(Log_11KT) + s(mass_final),
  family = binomial,
  data = subset(train, Phase %in% c('M','D','F') & !is.na(Log_11KT)),
  method = "REML"
)

train$pfemale_noE2 <- predict(
  gam_progress_noE2,
  newdata = train,
  type = "response"
)

test$pfemale_noE2 <- predict(
  gam_progress_noE2,
  newdata = test,
  type = "response"
)

ggplot(test, aes(x = Phase, y = pfemale_noE2))+
  geom_boxplot()+
  geom_point()
