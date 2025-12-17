dodd_2019 = read_excel('Measures/dodd_2019_data.xlsx')
library(forcats)


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

days_to_ymwd <- function(d) {
  y <- d %/% 365
  d <- d %% 365
  
  m <- d %/% 30
  d <- d %% 30
  
  w <- d %/% 7
  d <- d %% 7
  
  paste0(y, "y ", m, "m ", w, "w ", d, "d")
}

dodd_2019$time_label <- days_to_ymwd(dodd_2019$time_better)

ggplot(subset(dodd_2019, phase2!='S'), aes(x = time_label, y = log10(KT11), shape = phase2, color = phase2))+
  geom_point(size =3)

ggplot(subset(dodd_2019, phase2!='S'), aes(x = time_label, y = log10(Cort), shape = phase2, color = phase2))+
  geom_point(size =3)

ggplot(subset(dodd_2019, phase2!='S'), aes(x = time_label, y = log10(E2), shape = phase2, color = phase2))+
  geom_point(size =3)

ggplot(subset(dodd_2019, phase2%in%c('D','M') & E2 < 1000), aes(x = (E2), y = (KT11), shape = phase2, color = phase2))+
  geom_point(size =3)
# not correlated?

ggplot(subset(dodd_2019, phase2!='S'), aes(x = time_label, y = gPOAmedium, shape = phase2, color = phase2))+
  geom_point(size =3)
# i think this is supposed to be aPOA

dodd_2019$length2= dodd_2019$`length 2`
mop_gpoa = lm(gPOAarea ~phase2+length2, data = dodd_2019)
corrected_gpoa = predict(mop_gpoa, data.frame(phase2 = dodd_2019$phase2,
                                              length2 = mean(dodd_2019$`length 2`)))

dodd_2019$corrected_gpoa_area = corrected_gpoa

dodd_2019_mean = dodd_2019%>%
  group_by(phase2)%>%
  summarize(mean_gpoa = mean(corrected_gpoa_area, na.rm = T))


dodd_2019_mean$phase2 = factor(dodd_2019_mean$phase2, levels = c('M', 'D','F', 'S', 'NA'))

ggplot(subset(dodd_2019_mean, phase2%in% c('M','D','F')), aes(x = phase2, y = mean_gpoa, fill =phase2))+
  geom_bar(stat = 'identity', position = 'dodge')
# i think gPOA is a poa but need to talk to justin

#######
multiome_measures = read.csv('Measures/all_data.csv')

# step 1, fit a linear regression with days as timepoints assuming F = 4 yrs which is questionable
# lets do 11kt first

dodd_2019$time_better

dodd_2019_good = subset(dodd_2019, phase2!='S')

fm1 = lm(time_better~(KT11) , data = dodd_2019_good)
summary(fm1)

#prediction time

predicted_days = predict(fm1, data.frame(KT11 = 10^multiome_measures$Log_11KT ))
multiome_measures$predicted_days_kt = predicted_days
multiome_measures$percent_female_kt = predicted_days/(365*4)




ggplot(multiome_measures, aes(x = percent_female_kt, y =10^Log_11KT, color = Status))+
  geom_point()
multiome_measures$Status = factor(multiome_measures$Status, levels = c('NRM',
                                                                       'S',
                                                                       'EP',
                                                                       'NM',
                                                                       'M',
                                                                       'D',
                                                                       'E',
                                                                       'NF',
                                                                       'F'))
ggplot(multiome_measures, aes(x = Status, y = percent_female_kt))+
  geom_boxplot()+
  geom_point()

#confused as to why males are not 0

dodd_2019_good$Log_11KT = log10(dodd_2019_good$KT11)
dodd_2019_good$Log_E2 = log10(dodd_2019_good$E2)

range(dodd_2019_good$Log_11KT[!is.na(dodd_2019_good$Log_11KT)])
range(multiome_measures$Log_11KT[!is.na(multiome_measures$Log_11KT)])
# would it be wrong here to apply a scale factor? They are close but i feel like
# you could coerce them to be 1-1 

hist(dodd_2019_good$Log_11KT[!is.na(dodd_2019_good$Log_11KT)])
hist(multiome_measures$Log_11KT[!is.na(multiome_measures$Log_11KT)])
# they are very similar

# lets just scale
normalize <- function(x) {
  return((x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE)))
}
multiome_measures$coerced_percent = normalize(multiome_measures$predicted_days_kt)

ggplot(multiome_measures, aes(x = Status, y = coerced_percent))+
  geom_boxplot()+
  geom_point()

# can we predict estradiol from this?

fm_e2 = lm(log10(E2)~time_better , data = dodd_2019_good)

predicted_e2 = predict(fm_e2, data.frame(time_better = multiome_measures$predicted_days_kt ))
multiome_measures$predicted_e2_fromkt = predicted_e2


ggplot(multiome_measures, aes(x = Status, y = predicted_e2_fromkt))+
  geom_boxplot()+
  geom_point()
# yeah its just wrong, it does fit the true distribution from dodd

# ok the problem is that males and females are very similar, and the 4yr timepoint is problematic
# so, lets fit a model without females and then try again


fm1_nof = lm(time_better~log10(KT11) , data = subset(dodd_2019_good, phase2!='F'))
summary(fm1_nof)

predicted_days_nof = predict(fm1_nof, data.frame(KT11 = multiome_measures$Log_11KT ))

multiome_measures$predicted_days_kt_nof = predicted_days_nof

ggplot(multiome_measures, aes(x = Status, y = predicted_days_kt_nof))+
  geom_boxplot()+
  geom_point()

# still not great but I think its statistically more accurate

# i think I should probably include another variable, like percent testicular

hist((dodd_2019$pertest))
hist(multiom_measures$Percent_Testicular)

library(patchwork)
a = ggplot(multiome_measures, aes(x = Status, y = Percent_Testicular))+
  geom_boxplot()+
  geom_point()+
    ylim(0,1)

b=ggplot(dodd_2019,aes(x = phase2, y = pertest) )+
    geom_boxplot()+
  geom_point()+
  ylim(0,1)
a+b

mean(subset(dodd_2019, phase2 == 'M')$pertest, na.rm = TRUE)
mean(multiome_measures$Percent_Testicular[multiome_measures$Status=='M'])
# consistent


mean(subset(dodd_2019, phase2 == 'D')$pertest, na.rm = TRUE)
mean(multiome_measures$Percent_Testicular[multiome_measures$Status=='D'])
# closeish

# lets try 
dodd_2019_good$weight2 = dodd_2019_good$`weight 2`
fm_3 = lm(time_better~Log_11KT*pertest*weight2, data = subset(dodd_2019_good, phase2!='F'))
predicted_days_fm3 = predict(fm_3, data.frame(Log_11KT = multiome_measures$Log_11KT,
                                              pertest = multiome_measures$Percent_Testicular,
                             weight2 = multiome_measures$mass_final_cm))

multiome_measures$predicted_fm3 = predicted_days_fm3

 ggplot(multiome_measures, aes(x = Status, y = predicted_days_fm3))+
  geom_boxplot()+
  geom_point()
 
mean(subset(multiome_measures, Status == 'E')$predicted_fm3, na.rm = TRUE)
mean(subset(multiome_measures, Status == 'D')$predicted_fm3, na.rm = TRUE)
# man I wish I could get it to be close to 180
# you could lie and pyt pertest = 0 for f and nf and it wouldnt be wrong

multiome_measures$Percent_Testicular_lie = ifelse(multiome_measures$Status=='F' |multiome_measures$Status=='NF' , 0, multiome_measures$Percent_Testicular)
predicted_days_fm3_2 = predict(fm_3, data.frame(Log_11KT = multiome_measures$Log_11KT,
                                              pertest = multiome_measures$Percent_Testicular_lie,
                             weight2 = multiome_measures$mass_final_cm))

multiome_measures$predicted_fm3_2 = predicted_days_fm3_2

 ggplot(multiome_measures, aes(x = Status, y = predicted_days_fm3_2))+
  geom_boxplot()+
  geom_point()
 
# I mean this is alright
 female_mean =mean(subset(multiome_measures, Status == 'F')$predicted_fm3_2, na.rm = TRUE)
 
 multiome_measures$percent_female_complex = (multiome_measures$predicted_fm3_2/female_mean)

 ggplot(multiome_measures, aes(x = Status, y = percent_female_complex))+
  geom_boxplot()+
  geom_point()
# i guess its just not gonna work this way unfort
 
 
 
 


 
