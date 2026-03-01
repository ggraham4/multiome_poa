library(ggplot2)
library(patchwork)
library(tidyverse)
library(car)

multiome = read.csv('Measures/2025_12_26 all_data.csv')
multiome_good = subset(multiome, Phase %in% c('M','I','LI','NF','F'))
multiome_good$logit_test = car::logit(multiome_good$Percent_Testicular, adjust = 0.01)

mod1 = lm(SexState ~ logit_test, data = multiome_good)

multiome_good$mod1_resid = resid(mod1)
multiome_good$mod1_resid_squared = multiome_good$mod1_resid^2

ggplot(multiome_good, aes(x = SexState, y = mod1_resid_squared))+
  geom_point(aes(color = Phase))+
  geom_smooth(method ='gam')

ggplot(multiome_good, aes(x = SexState, y = mod1_resid))+
  geom_point(aes(color = Phase))

ggplot(multiome_good, aes(x = SexState, y = logit_test))+
  geom_point(aes(color = Phase))

# hm maybe this isnt the right approach,
# but I think I want to try anyway

### triangle 1 
i_residSquared = multiome_good$mod1_resid_squared[multiome_good$Phase=='I']
m_residSquared= multiome_good$mod1_resid_squared[multiome_good$Phase=='M']
i_sexState = multiome_good$SexState[multiome_good$Phase=='I']
m_sexState = multiome_good$SexState[multiome_good$Phase=='M']

m1 = (sd(i_residSquared)-sd(m_residSquared))/(mean(i_sexState)-mean(m_sexState))

# I just want to confirm I got the right slope
mod2 = lm(mod1_resid_squared~SexState, data = subset(multiome_good,Phase %in% c('M','I')))
mod2$coefficients[2] # aw man I got the wrong answer 
# i also forgot how to calculate yintercept

line_datm1 = data.frame(SexState = seq(min(multiome_good$SexState), max(multiome_good$SexState), by = 0.1))
line_datm1$y = mod2$coefficients[2]* line_datm1$SexState

ggplot(subset(multiome_good,Phase %in% c('M','I')), aes(x = SexState, y = mod1_resid_squared))+
  geom_point(aes(color = Phase))+
  geom_line(data =line_datm1, aes(x = SexState, y = y+mod2$coefficients[1]))

### anyway triangle time #####
i_residSquared = multiome_good$mod1_resid_squared[multiome_good$Phase=='I']
m_residSquared= multiome_good$mod1_resid_squared[multiome_good$Phase=='M']
f_residSquared =  multiome_good$mod1_resid_squared[multiome_good$Phase=='F']
i_sexState = multiome_good$SexState[multiome_good$Phase=='I']
m_sexState = multiome_good$SexState[multiome_good$Phase=='M']
f_sexState = multiome_good$SexState[multiome_good$Phase=='F']

mod2 = lm(mod1_resid_squared~SexState, data = subset(multiome_good,Phase %in% c('M','I')))
m1 = mod2$coefficients[2] 
b1 = mod2$coefficients[1] 

mod3 = lm(mod1_resid_squared~SexState, data = subset(multiome_good,Phase %in% c('I', 'F')))
m2 = mod3$coefficients[2]
b2= mod3$coefficients[1] 

x1 = mean(i_sexState)- mean(i_sexState)
x2 = mean(f_sexState) - mean(i_sexState)


sequence =seq(min(multiome_good$SexState), max(multiome_good$SexState), by = 0.1)
triangle_dat = data.frame()
for(i in sequence){
  if(i < mean(i_sexState)){
    y = (m1*i)+b1
  }
  if(i > mean(i_sexState)){
        y = (m2*i)+b2
  }
  if(i> mean(f_sexState)){break}
  newd = data.frame(SexState= i, y = y)
  triangle_dat = rbind(newd, triangle_dat)
}

ggplot(triangle_dat, aes(x =SexState, y = y))+
  geom_line()

ggplot(subset(multiome_good,Phase %in% c('I','F')), aes(x = SexState, y = mod1_resid_squared))+
  geom_point(aes(color = Phase))+
  geom_line(data =line_datm1, aes(x = SexState, y = (SexState*m2+b2)))
# this one mf is ruining my triangle


#variance squared?
multiome_good$logit_test_var2 = (multiome_good$logit_test- mean(multiome_good$logit_test, na.rm = T))^2

ggplot(multiome_good, aes(x = SexState, y = logit_test_var2))+
  geom_point()

# i think I am becoming schizophrenic
# what if we just look at the variance of the sex state

multiome_good$sexstate_var2 = (multiome_good$SexState- mean(multiome_good$SexState, na.rm = T))^2

ggplot(multiome_good, aes(x = SexState, y = sexstate_var2))+
  geom_point() # woah

# pretty sure IDK what Im doing and need to actually look up these equations

