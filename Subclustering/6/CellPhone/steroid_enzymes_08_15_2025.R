mecp = readRDS("Functions/mean_expression_cluster_plot.rds")

mecd = readRDS("Functions/mean_expression_cluster_data.rds")

mecp(obj, 'hsd11b1la',6, 'res0.8_50nn_40PC_45LSI')
mecp(obj, 'hsd11b1la',1, 'res0.8_50nn_40PC_45LSI')
#irreleveant, but should also make 11kt

mecp(obj, 'LOC111587899',1, 'res0.8_50nn_40PC_45LSI')
mecp(obj, 'LOC111587899',6, 'res0.8_50nn_40PC_45LSI')
#chol - preg

mecp(obj, 'hsd3b1',1, 'res0.8_50nn_40PC_45LSI')
mecp(obj, 'hsd3b1',6, 'res0.8_50nn_40PC_45LSI')
#  preg - prog
#7-hydroxypregnenolone into 17hydroxyprogesterone.

mecp(obj, 'LOC111579695',1, 'res0.8_50nn_40PC_45LSI')
mecp(obj, 'LOC111579695',6, 'res0.8_50nn_40PC_45LSI') #first expressed enzyme, but no diff
# preg to 17-hpreg

#In a next step, 17-hydroxypregnenolone and 17-hydroxyprogesterone can be 
#converted by CYP17 into dehydroepiandrosterone (DHEA) and androstenedione, respectively.

mecp(obj, 'hsd17b12b',1, 'res0.8_50nn_40PC_45LSI')
mecp(obj, 'hsd17b12a',1, 'res0.8_50nn_40PC_45LSI')

mecp(obj, 'hsd17b12b',6, 'res0.8_50nn_40PC_45LSI')
mecp(obj, 'hsd17b12a',6, 'res0.8_50nn_40PC_45LSI')
# make testosterone

mecp(obj, 'hsd11b2',6, 'res0.8_50nn_40PC_45LSI')
mecp(obj, 'hsd11b2',1, 'res0.8_50nn_40PC_45LSI')
#makes 11kt
mecp(obj, 'hsd17b14',6, 'res0.8_50nn_40PC_45LSI')
#makes dhea

mecp(obj, 'LOC111577263',1, 'res0.8_50nn_40PC_45LSI')
mecp(obj, 'LOC111577263',6, 'res0.8_50nn_40PC_45LSI') #AHA!, the enzyme that makes testosterone goes up, but SO DOES cyp
#cyp19a1b, makes e1 and e2

aro_6 = mecd(obj, 'LOC111577263',6, 'res0.8_50nn_40PC_45LSI')
test_6 = mecd(obj, 'hsd17b12b',6, 'res0.8_50nn_40PC_45LSI')
cor.test(aro_6$mean, test_6$mean) # aw man
plot(aro_6$mean, test_6$mean)

aro_6$test= test_6$mean
ggplot(aro_6,aes(x= test, y = mean, color =Status))+
  geom_point()+
  geom_smooth(method ='lm')
#truly no correlation lmfao

aro_1 = mecd(obj, 'LOC111577263',1, 'res0.8_50nn_40PC_45LSI')
aro_1$test= test_6$mean
ggplot(aro_1,aes(x= test, y = mean, color =Status))+
  geom_point()+
  geom_smooth(method ='lm')

#none of this is to say anything about the enzymes activities also so this is as far as i can tell a waste of time


tacr3a_6 = mecd(obj, 'tacr3a',6, 'res0.8_50nn_40PC_45LSI')
tacr3a_6$test = test_6$mean

ggplot(tacr3a_6,aes(x= test, y = mean, color =Status))+
  geom_point()+
  geom_smooth(method ='lm')
#seems to be correlated in males nad not others? my theory isnt lookin too hot

  test2 = mecd(obj, 'hsd17b12a',6, 'res0.8_50nn_40PC_45LSI')
  tacr3a_6$test2 =test2$mean
ggplot(tacr3a_6,aes(x= test2, y = mean, color =Status))+
  geom_point()+
  geom_smooth(method ='lm')

pgr_6 = mecd(obj, 'pgr',6, 'res0.8_50nn_40PC_45LSI')
tacr3a_6$pgr = pgr_6$mean

ggplot(tacr3a_6,aes(x= pgr, y = mean, color =Status))+
  geom_point()+
  geom_smooth(method ='lm')

hsd17b14_6 = mecd(obj, 'hsd17b14',6, 'res0.8_50nn_40PC_45LSI')
 mecp(obj, 'ptpra',6, 'res0.8_50nn_40PC_45LSI')

tacr3a_6$hsd = hsd17b14_6$mean
ggplot(tacr3a_6,aes(x= hsd, y = mean, color =Status))+
  geom_point()+
  geom_smooth(method ='lm')

