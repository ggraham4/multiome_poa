library(ggplot2)
library(lme4)
library(tidyverse)
library(emmeans)
library(ggsignif)
library(multcomp)
library(factoextra)
`%notin%` = Negate(`%in%`)
set.seed(0)

measures = read.csv('Measures/all_data.csv')
measures$Time_Day_2= as.numeric(measures$Time_Day_2)
measures$Behaviors_Day_2= as.numeric(measures$Behaviors_Day_2)
measures$Log_11KT = as.numeric(measures$Log_11KT)


### including this so it runs i guess.  #########
library(lavaan)
library(car)
library(ggplot2)
library(patchwork)

set.seed(0)
coalesced = read.csv('Measures/coalesced_data.csv')
multiome = read.csv('Measures/all_data.csv')

coalesced_good = subset(coalesced, Phase %in% c('M','D','F'))
coalesced_good$Percent_Testicular[coalesced_good$Phase=='F']=0

coalesced_good =coalesced_good[complete.cases(coalesced_good$Log_11KT,
                                              coalesced_good$Percent_Testicular, 
                                              coalesced_good$mass_final),]


multiome$mass_final =multiome$mass_final_cm
multiome$logit_testis = car::logit(multiome$Percent_Testicular, adjust = 0.01)
multiome$z_11kt <- scale(multiome$Log_11KT)
multiome$z_testis <- scale(multiome$logit_testis)
multiome$z_mass <- scale(multiome$mass_final)
multiome_good = subset(multiome, Status %in% c('M', 'D','E','NF','F'))
multiome_good$z_behavior <- scale(multiome_good$Behaviors_Day_2)
multiome_good$z_time<- scale(multiome_good$Time_Day_2)


# The "Squeezed Logit": adjust=0.01 handles the 0s and 1s safely.
coalesced_good$logit_testis <- car::logit(coalesced_good$Percent_Testicular, adjust = 0.01)

coalesced_good$z_11kt <- scale(coalesced_good$Log_11KT)
coalesced_good$z_testis <- scale(coalesced_good$logit_testis)
coalesced_good$z_mass <- scale(coalesced_good$mass_final)
# 2. Define the Model
# This makes comparing individuals much easier.
sex_model <- '
  # Measurement Model
  SexState =~ 1*z_11kt + z_testis + z_mass 
  
  # Variance of the latent variable fixed to 1 for identification
  SexState ~~ SexState
'

# 3. Fit the Model
# MLR (Robust Maximum Likelihood) is essential here to handle the 
# non-normal distribution of the logit-transformed data.
fit <- lavaan::sem(
  model = sex_model,
  data = coalesced_good,
  missing = "fiml",
  estimator = "MLR"
)
fit1a <- lavaan::sem( # for comparison to fit2
  model = sex_model,
  data = multiome_good,
  missing = "fiml",
  estimator = "MLR"
)


# validate the model is not bs
summary(fit, fit.measures = TRUE,standardized = TRUE)
# the scores are not great but god enough

#extract scores
latent_scores <- lavPredict(fit)
coalesced_good$SexState<- latent_scores[,1]*-1


####plotting #### 

ggplot(coalesced_good, aes(x = Phase, y = SexState))+
  geom_boxplot()+
  geom_point(aes(shape = source))

ggplot(coalesced_good, aes(x = SexState, y = logit_testis))+
  geom_point(aes(shape = source))

ggplot(coalesced_good, aes(x = SexState, y = Percent_Testicular))+
  geom_point(aes(shape = source))


ggplot(coalesced_good, aes(x = SexState, y = Log_11KT, color = Phase)) +
  geom_point(aes(shape = source))

ggplot(coalesced_good, aes(x = SexState, y = mass_final, color = Phase)) +
  geom_point(aes(shape = source))


# multiome
multiome_good$SexState = lavPredict(fit, multiome_good)*-1

ggplot(multiome_good, aes(x = SexState, y = Time_Day_2, color = Status)) +
  geom_point(aes())

ggplot(multiome_good, aes(x = SexState, y = Behaviors_Day_2, color = Status)) +
  geom_point(aes())

#include behavior?

sex_model2 <- '
  # Measurement Model
  SexState =~ 1*z_11kt + z_testis + z_mass +z_behavior + z_time
  
  # Variance of the latent variable fixed to 1 for identification
  SexState ~~ SexState
'

fit2 <- lavaan::sem(
  model = sex_model2,
  data = multiome_good,
  missing = "fiml",
  estimator = "MLR"
)

summary(fit2, fit.measures = TRUE,standardized = TRUE)

multiome_good$SexState2 = predict(fit2, multiome_good)*-1

ggplot(multiome_good, aes(x = SexState2, y = Time_Day_2, color = Status)) +
  geom_point(aes())
ggplot(multiome_good, aes(x = SexState, y = Time_Day_2, color = Status)) +
  geom_point(aes())

ggplot(multiome_good, aes(x = SexState2, y = Behaviors_Day_2, color = Status)) +
  geom_point(aes())
ggplot(multiome_good, aes(x = SexState, y = Behaviors_Day_2, color = Status)) +
  geom_point(aes())

anova(fit1a, fit2) # 1a is better than 2

# get regression lines
vars = c("z_11kt",
         "z_testis",
         "z_mass")

func_dat = data.frame()
for(i in vars){
  formula =as.formula(paste0('SexState~',i))
  model = lm(formula, data = coalesced_good)# there does seem to be an effect of
  #source with some of ross's females, so it may be prudent to 
  # regress out source, though also this could just be real who knows
  
  m =model$coefficients [2]
  b = model$coefficients[1]
  
  dat = data.frame(var = i,
                   slope = m,
                   intercept = b)
  
  func_dat = rbind(dat, func_dat)
  
}

state_from_zmass = function(z_mass){
  out = ((0.6016658*z_mass) -3.988931e-16)
  return(out)
}

ggplot(coalesced_good, aes(x = SexState, y = z_mass))+
  geom_point(aes(color = Phase, shape = source))+
  geom_line(aes(x = sapply(z_mass, state_from_zmass)))

state_from_z11kt = function(z_kt){
  out = ((-0.6062494*z_kt) -8.246849e-17)
  return(out)
}

ggplot(coalesced_good, aes(x = SexState, y = z_11kt))+
  geom_point(aes(color = Phase, shape = source))+
  geom_line(aes(x = sapply(z_11kt, state_from_z11kt)))

state_from_zlogitTestis = function(z_logitTestis){
  # input z_score transformed logit transformed % testicular
  out = ((-0.7067754*z_logitTestis) -1.088467e-16)
  return(out)
}

ggplot(coalesced_good, aes(x = SexState, y = z_testis))+
  geom_point(aes(color = Phase, shape = source))+
  geom_line(aes(x = sapply(z_testis, state_from_zlogitTestis)))
#perf

### try to get a curve for behavior?
ggplot(multiome_good, aes(x = SexState, y = z_time))+
  geom_point(aes(color = Status)) # looks parabolic to me

ggplot(multiome_good, aes(x = SexState, y = z_time))+
  geom_point(aes(color = Status))+
  geom_line(aes(x = sapply(z_time, function(z_time){
    return((-0.2*z_time^2)-(1*z_time)+1)}))) # good start for nls

z_time_nls = nls(SexState~(a*z_time^2)-(b*z_time)+c,
                 start = list(
                   a = -0.2,
                   b = -1,
                   c = 1
                 ),
                 data = multiome_good)

ggplot(multiome_good, aes(x = SexState, y = z_time))+
  geom_point(aes(color = Status))+
  geom_line(aes(x = sapply(z_time, function(z_time){
    return((-0.1167*z_time^2)-(0.5699*z_time)+0.3517)}))) # good start for nls
# really ugly but that is the result ig

state_from_z_time = function(z_time){
   return((-0.1167*z_time^2)-(0.5699*z_time)+0.3517)
}

#behavior
ggplot(multiome_good, aes(x = SexState, y = z_behavior))+
  geom_point(aes(color = Status)) # looks sigmoid? to me

ggplot(multiome_good, aes(x = SexState, y = z_behavior))+
  geom_point(aes(color = Status))+
  geom_line(aes(x = sapply(z_behavior, function(z_behavior){
    return(
      1/(.02+exp(-z_behavior*(1-2)))-1
    )}))) # very good start for nls

ggplot(multiome_good, aes(x = SexState, y = z_behavior))+
  geom_point(aes(color = Status))+
  geom_line(aes(x = sapply(z_behavior, function(z_behavior){
    return(
      1/(.02+exp(z_behavior))-1
    )}))) 

z_behavior_nls = nls(SexState~  1/(a+exp(z_behavior))-c,
                 start = list(
                   a = 0.2,
                   c = -1
                 ),
                 data = multiome_good)

ggplot(multiome_good, aes(x = SexState, y = z_behavior))+
  geom_point(aes(color = Status))+
  geom_line(aes(x = sapply(z_behavior, function(z_behavior){
    return(
      1/(0.2378+exp(z_behavior))-0.7376
    )})))

sigmoid_func <- function(x) {
  lower <- -0.5
  upper <- 2
  growth_rate <- -6
  inflection <- -0.7
  
  return(lower + ((upper - lower) / (1 + exp(-growth_rate * (x - inflection)))))
}

ggplot(multiome_good, aes(x = z_behavior, y =SexState )) +
  geom_point(aes(color = Status)) +
  # Use stat_function to draw the curve across the range of x
  stat_function(fun = sigmoid_func, color = "blue", size = 1)

z_behavior_nls_2 = nls(SexState~  lower + ((upper - lower) / (1 + exp(-growth_rate * (z_behavior - inflection)))),
                 start = list(
            lower = -0.5,
            upper = 2,
            growth_rate = -6,
            inflection = -0.7
                 ),
                 data = multiome_good)


sigmoid_func2 <- function(x) {
  lower <-  -0.3484  
  upper <- 0.9224
  growth_rate <- -6.1528  
  inflection <-  -0.4308 
  
  return(lower + ((upper - lower) / (1 + exp(-growth_rate * (x - inflection)))))
}

ggplot(multiome_good, aes(x = z_behavior, y = SexState)) +
  geom_point(aes(color = Status)) +
  # Use stat_function to draw the curve across the range of x
  stat_function(fun = sigmoid_func2, color = "blue", size = 1)
# ok its trying but this is just wrong , mine is better

temp_state_from_z_behavior = function(z_behavior){
  #nls disagrees with me but if you look at the plots youll agree with me
    lower <- -0.5
  upper <- 2
  growth_rate <- -6
  inflection <- -0.7
  
  return(lower + ((upper - lower) / (1 + exp(-growth_rate * (z_behavior - inflection)))))

  
}
#e2
  e2_mod = lm(SexState~scale(Log_E2), data = coalesced_good)
summary(e2_mod)

state_from_z_e2 = function(z_e2){
  out = (z_e2*0.45342)+0.18467
  return(out)
}

#back transform testis
state_from_zlogitTestis_inv <- function(SexState){
  a <- -0.7067754
  b <- -1.088467e-16
  z_logit <- (SexState - b) / a
  return(z_logit)
}

####other way around


z_kt_from_state = function(SexState){
    #out = ((-0.6062494*z_kt) -8.246849e-17)
  z_kt = (SexState/-0.6062494)+8.246849e-17
return(z_kt)
  
}

z_logitTestis_from_state = function(SexState){
#out = ((-0.7067754*z_logitTestis) -1.088467e-16)  
  z_logitTestis_from_state = (SexState/-0.7067754)+1.088467e-16
  return(z_logitTestis_from_state)
}

z_e2_from_state = function(SexState){
#  out = (z_e2*0.45342)+0.18467

  z_e2_from_state = (SexState/0.45342)-0.18467
  return(z_e2_from_state)
}

z_mass_from_state = function(SexState){
#    out = ((0.6016658*z_mass) -3.988931e-16)
  z_mass_from_state = (SexState/0.6016658)+3.988931e-16
  return(z_mass_from_state)
}

z_mass_from_state = function(SexState){
#    out = ((0.6016658*z_mass) -3.988931e-16)
  z_mass_from_state = (SexState/0.6016658)+3.988931e-16
  return(z_mass_from_state)
}

# Get Mean and SD of the Logit Variable (to un-scale z_testis)
logit_mu <- mean(coalesced_good$logit_testis, na.rm = TRUE)
logit_sd <- sd(coalesced_good$logit_testis, na.rm = TRUE)

# Get Mean and SD of the Original Percent Variable (to create the final z_score)
pct_mu <- mean(coalesced_good$Percent_Testicular, na.rm = TRUE)
pct_sd <- sd(coalesced_good$Percent_Testicular, na.rm = TRUE)

z_percent_testicular_from_state <- function(SexState) {
  # --- Step 1: SexState -> z_logit ---
  # Using the coefficients you derived earlier
  slope <- -0.7067754
  intercept <- -1.088467e-16 
  z_logit <- (SexState - intercept) / slope
  
  # --- Step 2: z_logit -> raw logit value (Un-scale) ---
  raw_logit <- (z_logit * logit_sd) + logit_mu
  
  # --- Step 3: raw logit -> raw percent (Inverse Logit) ---
  # Calculate standard logistic probability
  p_adj <- 1 / (1 + exp(-raw_logit))
  
  # Reverse the car::logit adjustment (adjust = 0.01)
  # car::logit formula is: logit = log(p_adj / (1 - p_adj))
  # where p_adj = p * (1 - 2*a) + a
  # So: p = (p_adj - a) / (1 - 2*a)
  adj <- 0.01
  raw_pct <- (p_adj - adj) / (1 - 2 * adj)
  
  # Optional: Clip values to 0-1 range to handle math overflow extremes
  raw_pct <- pmin(pmax(raw_pct, 0), 1)

  # --- Step 4: raw percent -> z_percent (Re-scale) ---
  z_final <- (raw_pct - pct_mu) / pct_sd
  
  return(z_final)
}

ggplot(coalesced_good, aes(x = SexState, y = scale(Percent_Testicular)))+
  geom_point()+
  geom_line(aes(x = SexState, y=z_percent_testicular_from_state(SexState)))

z_behavior_from_state <- function(SexState) {
  lower <- -0.5
  upper <- 2
  k <- -6
  x0 <- -0.7

  x0 - (1 / k) * log((upper - SexState) / (SexState - lower))
}

xbeh = seq(-2,2, by = 0.1)
ybeh = temp_state_from_z_behavior(seq(-2,2, by = 0.1))

xtim = seq(-2,2, by = 0.1)
ytim = state_from_z_time(xtim)


plot_data = data.frame(SexState = seq(min(multiome_good$SexState, na.rm = T),
                                    max(multiome_good$SexState, na.rm =T),by = 0.1)
                       )

plot_data$z_kt= sapply(plot_data$SexState,z_kt_from_state )
plot_data$z_logit_testis= sapply(plot_data$SexState,z_logitTestis_from_state )
plot_data$z_e2= sapply(plot_data$SexState,z_e2_from_state )
plot_data$z_mass= sapply(plot_data$SexState,z_mass_from_state )
#plot_data$z_time= sapply(plot_data$SexState,z_time_from_state )
plot_data$z_behavior= sapply(plot_data$SexState,z_behavior_from_state )
plot_data$z_testis= sapply(plot_data$SexState,z_percent_testicular_from_state )



a=ggplot(plot_data, aes(x = SexState, y = z_kt))+
  geom_line(aes(color = '11KT'))+
  geom_line(aes(color = 'Testicular %', y = z_testis))+
    geom_line(aes(color = 'E2', y = z_e2))+
  geom_line(aes(color= 'Mass', y =z_mass))+
      geom_line(data = data.frame(), aes(color= 'Time in Nest',x =xtim, y =ytim))+
      geom_line(data = data.frame(), aes(color= 'Parental Behavior',x =xbeh, y =ybeh))+
    xlim(min(multiome_good$SexState, na.rm = T),
                                    max(multiome_good$SexState, na.rm =T))+
  labs(y = 'Z Score', colour = 'Measure')


coalesced_good$Phase = factor(coalesced_good$Phase, levels = c('M','D','F'))
multiome_good$Status = factor(multiome_good$Status, levels = c('M','D','E','NF','F'))

b = ggplot(data = multiome_good, aes(x = SexState, y = Status))+
  geom_boxplot()+
  geom_point()+
  xlim(min(multiome_good$SexState, na.rm = T),
                                    max(multiome_good$SexState, na.rm =T))

(a+b)+
  plot_layout(ncol=1)


### #save work
#saveRDS(state_from_z_time, 'Functions/Theory/SexState_from_z_time.rds')
#saveRDS(state_from_z11kt, 'Functions/Theory/SexState_from_z_11kt.rds')
#saveRDS(state_from_zlogitTestis, 'Functions/Theory/SexState_from_z_logit_tesits.rds')
#saveRDS(z_percent_testicular_from_state, 'Functions/Theory/SexState_from_z_percent_tesits.rds')
#saveRDS(state_from_zmass, 'Functions/Theory/SexState_from_z_mass_tesits.rds')
#saveRDS(temp_state_from_z_behavior, 'Functions/Theory/SexState_temp_from_z_behavior.rds')
#saveRDS(state_from_z_e2, 'Functions/Theory/SexState_from_z_e2.rds')

multiome = read.csv('Measures/all_data.csv')
multiome=multiome %>%
  dplyr::rename(mass_final=mass_final_cm)

multiome$Percent_Testicular[multiome$Status %in%c('NF','F')]=0
multiome$z_11kt = scale(multiome$Log_11KT)[,1]
multiome$z_mass = scale(multiome$mass_final)[,1]
multiome$z_time = scale(multiome$Time_Day_2)[,1]
multiome$z_behavior = scale(multiome$Behaviors_Day_2)[,1]
multiome$z_percent_testicular = scale(multiome$Percent_Testicular)[,1]
multiome$z_testis = scale(car::logit(multiome$Percent_Testicular, adjust = 0.01))[,1]
multiome$z_volume = scale(multiome$Log10_Volume)[,1]
multiome$SexState = predict(fit, multiome)*-1


status_to_phase = list('D'='I',
                       'E' = 'LI',
                       'EP' = 'LIP',
                       'S' = 'IP',
                       'M'='M',
                       'F'='F',
                       'NF'='NF',
                       'NM'='NM')
multiome$Phase = unlist(status_to_phase[multiome$Status])

#write.csv(multiome, 'Measures/2025_12_26 all_data.csv')


coalesced =read.csv('Measures/coalesced_data.csv')

coalesced$Percent_Testicular[coalesced$Phase %in%c('NF','F')]=0
coalesced$z_11kt = scale(coalesced$Log_11KT)[,1]
coalesced$z_mass = scale(coalesced$mass_final)[,1]
#coalesced$z_time = scale(coalesced$Time_Day_2)[,1]
#coalesced$z_behavior = scale(coalesced$Behaviors_Day_2)[,1]
coalesced$z_percent_testicular = scale(coalesced$Percent_Testicular)[,1]
coalesced$z_testis = scale(car::logit(coalesced$Percent_Testicular, adjust = 0.01))[,1]
#coalesced$z_volume = scale(coalesced$Log10_Volume)[,1]
coalesced$SexState = predict(fit, coalesced)*-1

#write.csv(coalesced, 'Measures/2026_01_08 Coalesced sex state.csv')

            

##### variance by states #####
vart =multiome%>%
  subset(Status %in% c('M','D','E',
                       'F'))%>%
  group_by(Status)%>%
  summarize(testis_sample_var = var(Percent_Testicular)/(n()-1),
            mean_status = mean(SexState))

vart$Status = factor(vart$Status, levels = c('M','D','E',
                       'F'))
ggplot(vart, aes(y= testis_sample_var, x = Status))+
  geom_bar(stat = 'identity')

ggplot(vart, aes(y= testis_sample_var, x = mean_status))+
  geom_point(aes(color = Status))

## 
multiome_complete = multiome[,c('Log_11KT',
            'Percent_Testicular',
            'Time_Day_2',
            'Behaviors_Day_2',
            'mass_final', 'Status')]%>%
  na.omit()

multiome_complete[,1:5]%>%
  scale()%>%
  prcomp()%>%
  fviz_pca_biplot(habillage =multiome_complete$Status)
multiome_complete$Status= factor(multiome_complete$Status, levels = c('S',
                                                                      'EP',
                                                                      'NM',
                                                                      'M',
                                                                      'D',
                                                                      'E',
                                                                      'NF',
                                                                      'F'))
multiome_long_scaled =multiome_complete %>%
  mutate(across(c(Log_11KT, Percent_Testicular, Time_Day_2, 
                  Behaviors_Day_2, mass_final), scale)) %>%
  pivot_longer(values_to = 'value',
               names_to = 'measure',
               cols = c(Log_11KT, Percent_Testicular, Time_Day_2, 
                       Behaviors_Day_2, mass_final))%>%
  subset(Status %in% c('M',
                       'D',
                       'E',
                       'NF',
                       "F"))

ggplot(multiome_long_scaled, aes(x = Status, y = value))+
  geom_boxplot()+
  geom_point()+
  facet_wrap(~measure)
            
##### actual plots ####

status_to_phase = list('D'='I',
                       'E' = 'LI',
                       'EP' = 'LIP',
                       'S' = 'IP',
                       'M'='M',
                       'F'='F',
                       'NF'='NF',
                       'NM'='NM')
measures$Phase = unlist(status_to_phase[measures$Status])
measures$Phase = factor(measures$Phase, levels =c('F','M','I','IP','LI','LIP','NF','NM'))
measures$Dominance = ifelse(measures$Phase %in%c('F','I','LI','NF'), 'Dominant', 'Subordinate')
measures$Behaviors_Day_2 = as.numeric(measures$Behaviors_Day_2)
### Sum Behaviors ----
beh_model = lm(Behaviors_Day_2~Phase, data = measures)
anova(beh_model, test = 'Chisq')
p1 = pairs(emmeans(beh_model, 'Phase'), adjust = 'none')
cld_beh = cld(emmeans(beh_model, 'Phase'), Letters = letters, adjust = "none", alpha = 0.05)
cld_beh_df = as.data.frame(cld_beh)

mean(measures$Behaviors_Day_2[measures$Phase=='M'])/mean(measures$Behaviors_Day_2[measures$Phase=='F'])
mean(measures$Behaviors_Day_2[measures$Phase=='I'])/mean(measures$Behaviors_Day_2[measures$Phase=='IP'])
mean(measures$Behaviors_Day_2[measures$Phase=='LIP'])/mean(measures$Behaviors_Day_2[measures$Phase=='LI'])
mean(measures$Behaviors_Day_2[measures$Phase=='NM'])/mean(measures$Behaviors_Day_2[measures$Phase=='NF'])

beh_mean_sd = measures%>%
  group_by(Phase, Dominance)%>%
  summarize(mean  =mean(Behaviors_Day_2),
            se = sd(Behaviors_Day_2)/sqrt(n()))

beh_plot = ggplot(measures, aes(y = Behaviors_Day_2, x = Phase, color = Dominance, shape = Dominance))+
  stat_summary(data = subset(measures, !Phase %in% c('NM', 'NF')),
               geom = 'crossbar', fun = 'mean') +
  geom_errorbar(data = beh_mean_sd, aes(x = Phase, y = mean, ymin = mean-se, ymax = mean+se),
                color ='black',
                width = 0.4, linewidth = 0.5)+
  labs(x = 'Phase', y = 'Parental Behaviors') +
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(linewidth=1, color = 'black'),
        legend.key = element_blank())+
  geom_vline(xintercept = 2.5)+
    geom_vline(xintercept = 4.5,)+
      geom_vline(xintercept = 6.5)+
  annotate("segment", 
           x = 1.5, xend = 1.5, 
           y = mean(measures$Behaviors_Day_2[measures$Phase=='F']), 
           yend = mean(measures$Behaviors_Day_2[measures$Phase=='M']),
           arrow = arrow(ends = "both", length = unit(0.2, "cm")))+
  geom_text(x = 2, y = mean(c(mean(measures$Behaviors_Day_2[measures$Phase=='F']),
                              mean(measures$Behaviors_Day_2[measures$Phase=='M']))-20),
            label = '2.21x', 
            color = 'black',
            size = 2.5,
             legend.key = element_blank())+
    annotate("segment", 
           x = 3.5, xend = 3.5, 
           y = mean(measures$Behaviors_Day_2[measures$Phase=='I']), 
           yend = mean(measures$Behaviors_Day_2[measures$Phase=='IP']),
           arrow = arrow(ends = "both", length = unit(0.2, "cm")))+
    geom_text(x = 4, y = mean(c(mean(measures$Behaviors_Day_2[measures$Phase=='I']),
                              mean(measures$Behaviors_Day_2[measures$Phase=='IP']))),
            label = '3.88x', 
            color = 'black',
            size = 2.5,
             legend.key = element_blank())+
      annotate("segment", 
           x = 5.5, xend = 5.5, 
           y = mean(measures$Behaviors_Day_2[measures$Phase=='LI']), 
           yend = mean(measures$Behaviors_Day_2[measures$Phase=='LIP']),
           arrow = arrow(ends = "both", length = unit(0.2, "cm")))+
      geom_text(x = 6, y = mean(c(mean(measures$Behaviors_Day_2[measures$Phase=='LI']),
                              mean(measures$Behaviors_Day_2[measures$Phase=='LIP']))),
            label = '7.64x', 
            color = 'black',
            size = 2.5,
             legend.key = element_blank())+
        annotate("segment", 
           x = 7.5, xend = 7.5, 
           y = mean(measures$Behaviors_Day_2[measures$Phase=='NF']), 
           yend = mean(measures$Behaviors_Day_2[measures$Phase=='NM']),
           arrow = arrow(ends = "both", length = unit(0.2, "cm")))+
        geom_text(x = 8, y = mean(c(mean(measures$Behaviors_Day_2[measures$Phase=='NF']),
                              mean(measures$Behaviors_Day_2[measures$Phase=='NM']))),
            label = '11.08x', 
            color = 'black',
            size = 2.5,
             legend.key = element_blank())+
    geom_point( size = 2, position = position_jitterdodge(0.25)) +
  geom_text(data = cld_beh_df, aes(x = Phase, y = 1.1*max(measures$Behaviors_Day_2), label = .group), 
          size = 3, inherit.aes = FALSE)+
  theme(legend.position = 'none')
beh_plot

ggsave(plot = beh_plot,
       file = "beh_plot.svg",
       device = "svg",
       units = "in",
       width = 3,
       height = 2,
       path = "Manuscript/Plots/Fig.1/Behaviors")




### Time ----
tim_mean_sd = measures%>%
  group_by(Phase, Dominance)%>%
  summarize(mean  =mean(Time_Day_2),
            se = sd(Time_Day_2)/sqrt(n()))

time_model = lm(Time_Day_2~Phase, data = measures)
anova(time_model, test = 'Chisq')
p2 = pairs(emmeans(time_model, 'Phase'), adjust = 'none')
cld_results = cld(emmeans(time_model, 'Phase'), Letters = letters, adjust = "none", alpha = 0.05)
cld_df = as.data.frame(cld_results)

mean(measures$Time_Day_2[measures$Phase=='M'])/mean(measures$Time_Day_2[measures$Phase=='F'])
mean(measures$Time_Day_2[measures$Phase=='I'])/mean(measures$Time_Day_2[measures$Phase=='IP'])
mean(measures$Time_Day_2[measures$Phase=='LIP'])/mean(measures$Time_Day_2[measures$Phase=='LI'])
mean(measures$Time_Day_2[measures$Phase=='NM'])/mean(measures$Time_Day_2[measures$Phase=='NF'])

time_plot = ggplot(measures, aes(y = Time_Day_2, x = Phase, color = Dominance, shape = Dominance))+
  stat_summary(data = subset(measures, !Phase %in% c('NM', 'NF')),
               geom = 'crossbar', fun = 'mean') +
    geom_errorbar(data = tim_mean_sd, aes(x = Phase, y = mean, ymin = mean-se, ymax = mean+se),
                color ='black',
                width = 0.4, linewidth = 0.5)+
  labs(x = 'Phase', y = 'Time in nest (s)') +
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(linewidth=1, color = 'black'),
        legend.key = element_blank())+
  geom_vline(xintercept = 2.5)+
    geom_vline(xintercept = 4.5,)+
      geom_vline(xintercept = 6.5)+
  annotate("segment", 
           x = 1.5, xend = 1.5, 
           y = mean(measures$Time_Day_2[measures$Phase=='F']), 
           yend = mean(measures$Time_Day_2[measures$Phase=='M']),
           arrow = arrow(ends = "both", length = unit(0.2, "cm")))+
  geom_text(x = 2, y = mean(c(mean(measures$Time_Day_2[measures$Phase=='F']),
                              mean(measures$Time_Day_2[measures$Phase=='M']))),
            label = '1.69x', 
            color = 'black',
            size = 2.5,
             legend.key = element_blank())+
    annotate("segment", 
           x = 3.5, xend = 3.5, 
           y = mean(measures$Time_Day_2[measures$Phase=='I']), 
           yend = mean(measures$Time_Day_2[measures$Phase=='IP']),
           arrow = arrow(ends = "both", length = unit(0.2, "cm")))+
    geom_text(x = 4, y = mean(c(mean(measures$Time_Day_2[measures$Phase=='I']),
                              mean(measures$Time_Day_2[measures$Phase=='IP']))),
            label = '2.14x', 
            color = 'black',
            size = 2.5,
             legend.key = element_blank())+
      annotate("segment", 
           x = 5.5, xend = 5.5, 
           y = mean(measures$Time_Day_2[measures$Phase=='LI']), 
           yend = mean(measures$Time_Day_2[measures$Phase=='LIP']),
           arrow = arrow(ends = "both", length = unit(0.2, "cm")))+
      geom_text(x = 6, y = mean(c(mean(measures$Time_Day_2[measures$Phase=='LI']),
                              mean(measures$Time_Day_2[measures$Phase=='LIP']))),
            label = '1.84x', 
            color = 'black',
            size = 2.5,
             legend.key = element_blank())+
        annotate("segment", 
           x = 7.5, xend = 7.5, 
           y = mean(measures$Time_Day_2[measures$Phase=='NF']), 
           yend = mean(measures$Time_Day_2[measures$Phase=='NM']),
           arrow = arrow(ends = "both", length = unit(0.2, "cm")))+
        geom_text(x = 8, y = mean(c(mean(measures$Time_Day_2[measures$Phase=='NF']),
                              mean(measures$Time_Day_2[measures$Phase=='NM']))),
            label = '10.67x', 
            color = 'black',
            size = 2.5,
             legend.key = element_blank())+
    geom_point( size = 2, position = position_jitterdodge(0.25)) +
  geom_text(data = cld_df, aes(x = Phase, y = 1.1*max(measures$Time_Day_2), label = .group), 
          size = 3, inherit.aes = FALSE)+
  
  theme(legend.position = 'none')
time_plot

ggsave(plot = time_plot,
       file = "time_plot_nolegend.svg",
       device = "svg",
       units = "in",
       width = 3,
       height = 2,
       path = "Manuscript/Plots/Fig.1/Time")


### 11KT ----
kt = measures%>%
  group_by(Phase, Dominance)%>%
  subset(!is.na(Log_11KT))%>%
  summarize(mean  =mean(Log_11KT,),
            se = sd(Log_11KT)/sqrt(n()))

kt_model = lm(Log_11KT~Phase, data = measures)
anova(kt_model, test = 'Chisq')
p3 = pairs(emmeans(kt_model, 'Phase'), adjust = 'none')
cld_kt = cld(emmeans(kt_model, 'Phase'), Letters = letters, adjust = "none", alpha = 0.05)
cld_df_kt = as.data.frame(cld_kt)

  mean(measures$Log_11KT[measures$Phase=='M']%>%na.omit())/mean(measures$Log_11KT[measures$Phase=='F']%>%na.omit())
mean(measures$Log_11KT[measures$Phase=='I'])/mean(measures$Log_11KT[measures$Phase=='IP'])
mean(measures$Log_11KT[measures$Phase=='LIP'])/mean(measures$Log_11KT[measures$Phase=='LI'])
mean(measures$Log_11KT[measures$Phase=='NM'])/mean(measures$Log_11KT[measures$Phase=='NF'])

kt_plot = ggplot(measures, aes(y = Log_11KT, x = Phase, color = Dominance, shape = Dominance))+
  stat_summary(data = subset(measures, !Phase %in% c('NM', 'NF')),
               geom = 'crossbar', fun = 'mean') +
    geom_errorbar(data = kt, aes(x = Phase, y = mean, ymin = mean-se, ymax = mean+se),
                color ='black',
                width = 0.4, linewidth = 0.5)+
  labs(x = 'Phase', y = 'Log10 11KT (pg/ml)') +
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(linewidth=1, color = 'black'),
        legend.key = element_blank())+
  geom_vline(xintercept = 2.5)+
    geom_vline(xintercept = 4.5,)+
      geom_vline(xintercept = 6.5)+
  annotate("segment", 
           x = 1.5, xend = 1.5, 
           y = mean(measures$Log_11KT[measures$Phase=='F']%>%na.omit()), 
           yend = mean(measures$Log_11KT[measures$Phase=='M']%>%na.omit()),
           arrow = arrow(ends = "both", length = unit(0.2, "cm")))+
  geom_text(x = 2, y = mean(c(mean(measures$Log_11KT[measures$Phase=='F' & !is.na(measures$Log_11KT)]),
                              mean(measures$Log_11KT[measures$Phase=='M'& !is.na(measures$Log_11KT)]))),
            label = '1.85x', 
            color = 'black',
            size = 2.5,
             legend.key = element_blank())+
    annotate("segment", 
           x = 3.5, xend = 3.5, 
           y = mean(measures$Log_11KT[measures$Phase=='I']), 
           yend = mean(measures$Log_11KT[measures$Phase=='IP']),
           arrow = arrow(ends = "both", length = unit(0.2, "cm")))+
    geom_text(x = 4, y = 0.05+mean(c(mean(measures$Log_11KT[measures$Phase=='I']),
                              mean(measures$Log_11KT[measures$Phase=='IP']))),
            label = '1.11x', 
            color = 'black',
            size = 2.5,
             legend.key = element_blank())+
      annotate("segment", 
           x = 5.5, xend = 5.5, 
           y = mean(measures$Log_11KT[measures$Phase=='LI']), 
           yend = mean(measures$Log_11KT[measures$Phase=='LIP']),
           arrow = arrow(ends = "both", length = unit(0.2, "cm")))+
      geom_text(x = 6, y = mean(c(mean(measures$Log_11KT[measures$Phase=='LI']),
                              mean(measures$Log_11KT[measures$Phase=='LIP']))),
            label = '1.69x', 
            color = 'black',
            size = 2.5,
             legend.key = element_blank())+
        annotate("segment", 
           x = 7.5, xend = 7.5, 
           y = mean(measures$Log_11KT[measures$Phase=='NF']), 
           yend = mean(measures$Log_11KT[measures$Phase=='NM']),
           arrow = arrow(ends = "both", length = unit(0.2, "cm")))+
        geom_text(x = 8, y = mean(c(mean(measures$Log_11KT[measures$Phase=='NF']),
                              mean(measures$Log_11KT[measures$Phase=='NM']))),
            label = '1.35x', 
            color = 'black',
            size = 2.5,
             legend.key = element_blank())+
    geom_point( size = 2, position = position_jitterdodge(0.25)) +
  geom_text(data = cld_df_kt, aes(x = Phase, y = 1.1*max(measures$Log_11KT%>%na.omit()), label = .group), 
          size = 3, inherit.aes = FALSE)+
  
  theme(legend.position = 'none')
kt_plot

ggsave(plot = kt_plot,
       file = "kt_plot.svg",
       device = "svg",
       units = "in",
       width = 3,
       height = 2,
       path = "Manuscript/Plots/Fig.1/kt")


###---- Percent Testicular -----

pct_test = measures%>%
  group_by(Phase, Dominance)%>%
  subset(!is.na(Percent_Testicular))%>%
  summarize(mean  =mean(Percent_Testicular,),
            se = sd(Percent_Testicular)/sqrt(n()))

test_model = lm(Percent_Testicular~Phase, data = measures)
anova(test_model, test = 'Chisq')
p4 = pairs(emmeans(test_model, 'Phase'), adjust = 'none')
cld_pct_test = cld(emmeans(test_model, 'Phase'), Letters = letters, adjust = "none", alpha = 0.05)
cld_pct_test_df = as.data.frame(cld_pct_test)

  #mean(measures$Percent_Testicular[measures$Phase=='M']%>%na.omit())/mean(measures$Log_11KT[measures$Phase=='F']%>%na.omit())
mean(measures$Percent_Testicular[measures$Phase=='I'])/mean(measures$Percent_Testicular[measures$Phase=='IP'])
mean(measures$Percent_Testicular[measures$Phase=='LIP'])/mean(measures$Percent_Testicular[measures$Phase=='LI'])
#mean(measures$Percent_Testicular[measures$Phase=='NM'])/mean(measures$Percent_Testicular[measures$Phase=='NF'])

test_plot = ggplot(subset(measures, !Phase %in%c('F','NF')), aes(y = Percent_Testicular, x = as.character(Phase)))+
  stat_summary(data = subset(measures, !Phase %in% c('NM', 'NF')),
               geom = 'crossbar', fun = 'mean') +
    geom_errorbar(data = pct_test, aes(x = as.character(Phase), y = mean, ymin = mean-se, ymax = mean+se),
                color ='black',
                width = 0.4, linewidth = 0.5)+
  labs(x = 'Phase', y = 'Testicular %') +
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(linewidth=1, color = 'black'),
        legend.key = element_blank())+
    geom_point( size = 2, position = position_jitterdodge(0.25), shape = 1) +
  geom_text(data = cld_pct_test_df, aes(x = Phase, y = 1.1*max(measures$Percent_Testicular%>%na.omit()), label = .group), 
          size = 3, inherit.aes = FALSE)+
  theme(legend.position = 'none')+
    scale_x_discrete(limits = c('M','I','IP','LI','LIP','NM')) +
        scale_y_continuous(labels = scales::percent)
test_plot

ggsave(plot = test_plot,
       file = "test_plot.svg",
       device = "svg",
       units = "in",
       width = 2.5,
       height = 2,
       path = "Manuscript/Plots/Fig.1/Testicular")

#### ovarian ----
pct_ov = measures%>%
  group_by(Phase, Dominance)%>%
  subset(!is.na(Percent_Ovarian))%>%
  summarize(mean  =mean(Percent_Ovarian,),
            se = sd(Percent_Ovarian)/sqrt(n()))

ov_model = lm(Percent_Ovarian~Phase, data = measures)
anova(ov_model, test = 'Chisq')
p5 = pairs(emmeans(ov_model, 'Phase'), adjust = 'none')
cld_pct_ov = cld(emmeans(ov_model, 'Phase'), Letters = letters, adjust = "none", alpha = 0.05)
cld_pct_ov_df = as.data.frame(cld_pct_ov)

  #mean(measures$Percent_Testicular[measures$Phase=='M']%>%na.omit())/mean(measures$Log_11KT[measures$Phase=='F']%>%na.omit())
mean(measures$Percent_Ovarian[measures$Phase=='I'])/mean(measures$Percent_Ovarian[measures$Phase=='IP'])
mean(measures$Percent_Ovarian[measures$Phase=='LIP'])/mean(measures$Percent_Ovarian[measures$Phase=='LI'])
#mean(measures$Percent_Testicular[measures$Phase=='NM'])/mean(measures$Percent_Testicular[measures$Phase=='NF'])

ov_plot = ggplot(subset(measures, !Phase %in%c('F','NF')), aes(y = Percent_Ovarian, x = as.character(Phase)))+
  stat_summary(data = subset(measures, !Phase %in% c('NM', 'NF')),
               geom = 'crossbar', fun = 'mean') +
    geom_errorbar(data = pct_ov, aes(x = as.character(Phase), y = mean, ymin = mean-se, ymax = mean+se),
                color ='black',
                width = 0.4, linewidth = 0.5)+
  labs(x = 'Phase', y = 'Ovarian %') +
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(linewidth=1, color = 'black'),
        legend.key = element_blank())+
    geom_point( size = 2, position = position_jitterdodge(0.25), shape = 1) +
  geom_text(data = cld_pct_ov_df, aes(x = Phase, y = 1.1*max(measures$Percent_Ovarian%>%na.omit()), label = .group), 
          size = 3, inherit.aes = FALSE)+
  theme(legend.position = 'none')+
    scale_x_discrete(limits = c('M','I','IP','LI','LIP','NM')) +
        scale_y_continuous(labels = scales::percent)
ov_plot

ggsave(plot = ov_plot,
       file = "ov_plot.svg",
       device = "svg",
       units = "in",
       width = 2.5,
       height = 2,
       path = "Manuscript/Plots/Fig.1/Ovarian")

### volume ----
measures$Log10_Volume <- as.numeric(measures$Log10_Volume)
ev_2.5x_model <- lm(Log10_Volume~Phase+length_final_cm, data = measures)
ev_2.5x_pairs <- as.data.frame(pairs(emmeans(ev_2.5x_model, 'Phase'), adjust = 'none'))

measures_gonad= subset(measures, Phase %notin% c('NF','F'))
  
residual_model <- lm(Log10_Volume~length_final_cm, data = measures_gonad)

measures_gonad$corrected_length =residual_model$residuals
anova(ev_2.5x_model,test = 'Chisq')

corrected_vol = measures_gonad%>%
  group_by(Phase)%>%
  subset(!is.na(corrected_length))%>%
  summarize(mean  =mean(corrected_length,),
            se = sd(corrected_length)/sqrt(n()))

cld_pct_vol = cld(emmeans(ev_2.5x_model, 'Phase'), Letters = letters, adjust = "none", alpha = 0.05)
cld_pct_vol_df = as.data.frame(cld_pct_vol)

vol_plot = ggplot(subset(measures_gonad, !Phase %in%c('F','NF')), aes(y = corrected_length, x = as.character(Phase)))+
  stat_summary(data = subset(measures_gonad, !Phase %in% c('NM', 'NF')),
               geom = 'crossbar', fun = 'mean') +
    geom_errorbar(data = corrected_vol, aes(x = as.character(Phase), y = mean, ymin = mean-se, ymax = mean+se),
                color ='black',
                width = 0.4, linewidth = 0.5)+
  labs(x = 'Phase', y = 'Corrected Volume') +
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(linewidth=1, color = 'black'),
        legend.key = element_blank())+
    geom_point( size = 2, position = position_jitterdodge(0.25), shape = 1) +
  geom_text(data = cld_pct_vol_df, aes(x = Phase,
                                      y = 1.15*max(measures_gonad$corrected_length%>%na.omit()),
                                      label = .group), 
          size = 3, inherit.aes = FALSE)+
  theme(legend.position = 'none')+
    scale_x_discrete(limits = c('M','I','IP','LI','LIP','NM')) 
vol_plot

ggsave(plot = vol_plot,
       file = "vol_plot.svg",
       device = "svg",
       units = "in",
       width = 2.5,
       height = 2,
       path = "Manuscript/Plots/Fig.1/Volume")



###
len = measures%>%
  group_by(Phase, Dominance)%>%
  subset(!is.na(length_final_cm))%>%
  summarize(mean  =mean(length_final_cm,),
            se = sd(length_final_cm)/sqrt(n()))

len_plot =ggplot(measures, aes(y = length_final_cm, x = Phase, color = Dominance, shape = Dominance))+
  stat_summary(data = subset(measures, !Phase %in% c('NM', 'NF')),
               geom = 'crossbar', fun = 'mean') +
    geom_errorbar(data = len, aes(x = Phase, y = mean, ymin = mean-se, ymax = mean+se),
                color ='black',
                width = 0.4, linewidth = 0.5)+
  labs(x = 'Phase', y = 'Length (cm)') +
  theme_bw()+
    geom_point( size = 2, position = position_jitterdodge(0.25)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(linewidth=1, color = 'black'),
        legend.key = element_blank())+
    geom_vline(xintercept = 2.5)+
    geom_vline(xintercept = 4.5,)+
      geom_vline(xintercept = 6.5)+
      theme(legend.position = 'none')


ggsave(plot = len_plot,
       file = "len_plot.svg",
       device = "svg",
       units = "in",
       width = 2.5,
       height = 2,
       path = "Manuscript/Plots/Fig.1/Length")



mass = measures%>%
  group_by(Phase, Dominance)%>%
  subset(!is.na(mass_final_cm))%>%
  summarize(mean  =mean(mass_final_cm,),
            se = sd(mass_final_cm)/sqrt(n()))

massplot = ggplot(measures, aes(y = mass_final_cm, x = Phase, color = Dominance, shape = Dominance))+
  stat_summary(data = subset(measures, !Phase %in% c('NM', 'NF')),
               geom = 'crossbar', fun = 'mean') +
    geom_errorbar(data = mass, aes(x = Phase, y = mean, ymin = mean-se, ymax = mean+se),
                color ='black',
                width = 0.4, linewidth = 0.5)+
  labs(x = 'Phase', y = 'Mass (g)') +
  theme_bw()+
    geom_point( size = 2, position = position_jitterdodge(0.25)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(linewidth=1, color = 'black'),
        legend.key = element_blank())+
    geom_vline(xintercept = 2.5)+
    geom_vline(xintercept = 4.5,)+
      geom_vline(xintercept = 6.5)+
    theme(legend.position = 'none')


ggsave(plot = massplot,
       file = "mass_plot.svg",
       device = "svg",
       units = "in",
       width = 2.5,
       height = 2,
       path = "Manuscript/Plots/Fig.1/mass")

##### latent variable ##### 
state_from_z_time = readRDS('Functions/Theory/SexState_from_z_time.rds')
state_from_z11kt=readRDS( 'Functions/Theory/SexState_from_z_11kt.rds')
state_from_zlogitTestis=readRDS( 'Functions/Theory/SexState_from_z_logit_tesits.rds')
z_percent_testicular_from_state =readRDS('Functions/Theory/SexState_from_z_percent_tesits.rds')
state_from_zmass =readRDS('Functions/Theory/SexState_from_z_mass_tesits.rds')
temp_state_from_z_behavior=readRDS( 'Functions/Theory/SexState_temp_from_z_behavior.rds')
state_from_z_e2=readRDS( 'Functions/Theory/SexState_from_z_e2.rds')



plot_data = data.frame(SexState = seq(min(multiome_good$SexState, na.rm = T),
                                    max(multiome_good$SexState, na.rm =T),by = 0.1)
                       )


plot_data$z_kt= sapply(plot_data$SexState,z_kt_from_state )
plot_data$z_logit_testis= sapply(plot_data$SexState,z_logitTestis_from_state )
plot_data$z_e2= sapply(plot_data$SexState,z_e2_from_state )
plot_data$z_mass= sapply(plot_data$SexState,z_mass_from_state )
#plot_data$z_time= sapply(plot_data$SexState,z_time_from_state )
plot_data$z_behavior= sapply(plot_data$SexState,z_behavior_from_state )
plot_data$z_testis= sapply(plot_data$SexState,z_percent_testicular_from_state )


xbeh = seq(-2,2, by = 0.1)
ybeh = temp_state_from_z_behavior(seq(-2,2, by = 0.1))

xtim = seq(-2,2, by = 0.1)
ytim = state_from_z_time(xtim)


neuro = function(SexState){
 y = exp(1)^(SexState)
 return(y-1.5)
}

z_neuro =sapply(xtim,neuro )

plot_data$z_neuro= sapply(plot_data$SexState,neuro )


a=ggplot(plot_data, aes(x = SexState, y = z_kt))+
  geom_line(aes(color = '11KT'), linewidth =2)+
  geom_line(aes(color = '% Testicular', y = z_testis), linewidth =2)+
    geom_line(aes(color = 'E2', y = z_e2), linewidth =2)+
  geom_line(aes(color= 'Mass', y =z_mass), linewidth =2)+
      geom_line(data = data.frame(), aes(color= 'Time in Nest',x =xtim, y =ytim), linewidth =2)+
        geom_line(data = data.frame(), aes(color= 'Neurogenesis',x =xtim, y =z_neuro), linewidth =2)+
      geom_line(data = data.frame(), aes(color= 'Parental Behavior',x =xbeh, y =ybeh), linewidth =2)+
    xlim(min(multiome_good$SexState, na.rm = T),
                                    max(multiome_good$SexState, na.rm =T))+
  labs(y = 'Arbitrary Units', x = '', color = 'Measurement')+
  theme_minimal()+
  ylim(-2, 4)+
  
 theme(axis.text.y = element_blank(), axis.text.x = element_blank())
a

ggsave(plot = a,
       file = "schematic.svg",
       device = "svg",
       units = "in",
       width = 6.7,
       height = 1.5,
       path = '/Users/ggraham/Desktop/multiome_poa/Manuscript/Plots/Manuscript v1.2')







