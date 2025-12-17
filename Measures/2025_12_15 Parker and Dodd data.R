library(forcats)
library(tidyverse)
library(ggplot2)
library(lme4)
library(readxl)
library(patchwork)
library(mgcv)

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

dodd_2019_good = subset(dodd_2019, phase2!='S')
dodd_2019_good$Log_11KT = log10(dodd_2019_good$KT11)
dodd_2019_good$Log_E2 = log10(dodd_2019_good$E2)
dodd_2019_good$weight2 = dodd_2019_good$`weight 2`


multiome_measures = read.csv('Measures/all_data.csv')
multiome_measures$Status = factor(multiome_measures$Status, levels = c('NRM',
                                                                       'S',
                                                                       'EP',
                                                                       'NM',
                                                                       'M',
                                                                       'D',
                                                                       'E',
                                                                       'NF',
                                                                       'F'))


### gam

# exclude females when fitting hormone relationships
dodd_prog <- subset(dodd_2019_good, phase2 %in% c("M","D",'F'))

# define weak labels
dodd_prog$progress = ifelse(dodd_prog$phase2=="M", 0,
           ifelse(dodd_prog$phase2=="F", 1, NA))

gam_data = dodd_prog
male_data_mult = subset(multiome_measures, Status %in% c('M', 'D', "E",'F'))
male_data_mult$progress = ifelse(male_data_mult$Status=="M", 0,
           ifelse(male_data_mult$Status=="F", 1, NA))

mult_prog = data.frame(phase2 = male_data_mult$Status,
                       Log_11KT= male_data_mult$Log_11KT,
                       pertest = male_data_mult$Percent_Testicular,
                       weight2 = male_data_mult$mass_final_cm,
                       progress = male_data_mult$progress,
                       source = 'multiome',
                       fish = male_data_mult$Fish
                       )

dodd_prog$source = 'Dodd'
dodd_prog$fish = dodd_prog$FishID

gam_train = rbind(mult_prog, dodd_prog[colnames(dodd_prog)%in%colnames(mult_prog)])

gam_train = gam_train%>%
  subset(!is.na(Log_11KT))

nrow(gam_train[gam_train$phase2%in%c('M','F'),])
sum(!is.na(gam_train$progress))
gam_train$pertest[gam_train$phase2=='F']=0

gam_prog <- gam(progress ~ 
                s(Log_11KT) + 
                s(pertest, k = 5) + 
                s(weight2),
                data = gam_train,
                family = binomial,
                  method = "REML")

multiome_DE <- subset(gam_train, phase2 %in% c("D","E") & source =='multiome' )

# predict probabilities
multiome_DE$predicted_progress <- predict(gam_prog, newdata = multiome_DE, type = "response")


ggplot(multiome_DE, aes(x = phase2, y = predicted_progress, color  = phase2))+
  geom_point(position = position_jitterdodge(2,0.04,seed = 2))+
  geom_text(aes(label = fish),position= position_jitterdodge(2,0.04, seed = 2), color = 'black')


# i'm goin beast mode with it and bringing parker 2023 into the mix ####
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

graham_2026 =  read.csv('Measures/all_data.csv')
graham_2026$mass_final =graham_2026$mass_final_cm

common_cols <- intersect(colnames(parker_2023), colnames(dodd_2019))
parker_dodd = rbind(parker_2023[,common_cols], dodd_2019[,common_cols])

parker_dodd$Phase = factor(parker_dodd$Phase, levels = c('S',
                                                         'M',
                                                         'D',
                                                         'NF',
                                                         'F'))
ggplot(parker_dodd, aes(x = Phase, y = Log_11KT, shape = source, color = source))+
  geom_point(position = position_dodge(1))+
  geom_boxplot(alpha=0)

ggplot(subset(parker_dodd, Phase =='D'), aes(x = days, y = Log_11KT, shape = source, color = source))+
  geom_point(position = position_dodge(1))+
  xlim(0, 350)

ggplot(subset(parker_dodd, Phase =='D'), aes(x = days, y = Log_E2, shape = source, color = source))+
  geom_point(position = position_dodge(1))+
  xlim(0, 350)

ggplot(subset(parker_dodd, Phase =='D'), aes(x = Log_11KT, y = Log_E2, shape = source, color = source, size = days))+
  geom_point(position = position_dodge(1))

parker_dodd$kt_e2 =parker_dodd$Log_11KT/parker_dodd$Log_E2

parker_dodd$days[parker_dodd$Phase=='M']=0
parker_dodd$days[parker_dodd$Phase=='F']=NA

plop_1 = ggplot(subset(parker_dodd, Phase != 'S'), aes(x = days, y = kt_e2))+
  geom_point(position = position_dodge(1),aes( shape = Phase, color = source))+
  xlim(0, 380)+
  geom_smooth()+
    ylim(0, 2)

plop_2 = ggplot(subset(parker_dodd, Phase == 'F'), aes(x='F', y = kt_e2))+
  geom_point(aes( shape = Phase, color = source), shape=3, size = 4)+
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.x = element_blank())+
  ylim(0, 2)

plop_1+plop_2+ plot_layout(guides = "collect") 

### I wonder if I can somehow coerce a curve ####

#add in multiome_data
multiome_measures = read.csv('Measures/all_data.csv')
multiome_measures$Time_Day_2= as.numeric(multiome_measures$Time_Day_2)
multiome_measures$Behaviors_Day_2= as.numeric(multiome_measures$Behaviors_Day_2)
multiome_measures$Log_11KT = as.numeric(multiome_measures$Log_11KT)

status_to_phase = list('D'='I',
                       'E' = 'LI',
                       'EP' = 'LIP',
                       'S' = 'IP',
                       'M'='M',
                       'F'='F',
                       'NF'='NF',
                       'NM'='NM')
multiome_measures = multiome_measures[!is.na(multiome_measures$Status),]

#multiome_measures$Phase = unlist(status_to_phase[multiome_measures$Status])
#multiome_measures$Phase = factor(multiome_measures$Phase, levels =c('F','M','I','IP','LI','LIP','NF','NM'))

multiome_measures$Phase = multiome_measures$Status

graham_2026 = multiome_measures
graham_2026$mass_final =graham_2026$mass_final
graham_2026$source ='Graham'

graham_2026$days = ifelse(graham_2026$Phase %in% c('E','I','NF','NM','IP','S'),
                          180,
                          NA)
graham_2026$days = ifelse(graham_2026$Phase %in% c('M'),
                          0,
                          graham_2026$days)


parker_dodd_graham = rbind(parker_dodd[,common_cols], graham_2026[,common_cols])%>%
  subset(!Phase %in% c('S','EP', 'NM'))

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
deangelis_2018$source = 'DeAngelis'
deangelis_2018=subset(deangelis_2018, Sex %in% c('f', 'm'))
deangelis_2018$Phase = ifelse(deangelis_2018$Sex=='f', 'F','M' )
deangelis_2018$days = ifelse(deangelis_2018$Phase=='M', 0,NA )

parker_dodd_graham_gonzalez_deangelis = rbind(parker_dodd_graham[,common_cols], 
                                    gonzalez_2021[,common_cols],
                                    deangelis_2018[,common_cols])

parker_dodd_graham_gonzalez_deangelis$mass_final = as.numeric(parker_dodd_graham_gonzalez_deangelis$mass_final)

coalesced = parker_dodd_graham_gonzalez_deangelis

ggplot(subset(coalesced, Phase != 'S'), aes(x = days, y = mass_final))+
  geom_point(position = position_dodge(1),aes( shape = Phase, color = source), size =3)+
  xlim(-1, 380)+
  geom_smooth()

coalesced$kt_e2 = coalesced$Log_11KT/coalesced$Log_E2

ggplot(subset(coalesced, Phase != 'S'), aes(x = kt_e2, y = mass_final))+
  geom_point(position = position_dodge(1),aes( shape = Phase, color = source), size = 3)+
  geom_smooth()

p1 = ggplot(subset(coalesced, Phase != 'S'), aes(x = days, y = kt_e2))+
  geom_point(position = position_dodge(1),aes( shape = Phase, color = source), size = 3)+
  geom_smooth()+
      ylim(0, 2)


p2=ggplot(subset(coalesced, Phase == 'F'), aes(x='F', y = kt_e2))+
  geom_point(aes( shape = Phase, color = source), shape=3, size = 4)+
      ylim(0, 2)

p1+p2

ggplot(subset(coalesced, Phase != 'S'), aes(x = Log_E2, y = mass_final))+
  geom_point(position = position_dodge(1),aes( shape = Phase, color = source), size = 4)+
  geom_smooth()

ggplot(subset(coalesced, Phase != 'S'), aes(x = Log_11KT, y = mass_final))+
  geom_point(position = position_dodge(1),aes( shape = Phase, color = source), size = 3)+
  geom_smooth()

