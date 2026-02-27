k1 = 0.4
k2 = 0.4
k3 = 0.4

simulate_model <- function(k1, k2, k3, steps = 20, P1_init = 100) {
  
  P1 <- P1_init
  P2 <- 0
  P3 <- 0
  N  <- 0
  
  for (t in 1:steps) {
    
    # transitions
    move12 <- P1 * k1
    move23 <- P2 * k2
    move3N <- P3 * k3
    
    # update states
    P1 <- P1 - move12
    P2 <- P2 + move12 - move23
    P3 <- P3 + move23 - move3N
    N  <- N  + move3N
  }
  
  return(c(P1=P1, P2=P2, P3=P3, N=N))
}

male <- simulate_model(k1=0.4, k2=0.4, k3=0.4)
male

intermediate <- simulate_model(k1=0.6, k2=0.3, k3=0.4)
intermediate

female <- simulate_model(k1=0.4, k2=0.3, k3=0.7)
female

to_percent <- function(x) {
  x / sum(x) * 100
}

simulate_model <- function(k1, k2, k3, steps = 20, P1_init = 100) {
  
  P1 <- P1_init
  P2 <- 0
  P3 <- 0
  N  <- 0
  
  for (t in 1:steps) {
    move12 <- P1 * k1
    move23 <- P2 * k2
    move3N <- P3 * k3
    
    P1 <- P1 - move12
    P2 <- P2 + move12 - move23
    P3 <- P3 + move23 - move3N
    N  <- N  + move3N
  }
  
  return(c(P1=P1, P2=P2, P3=P3, N=N))
}

to_percent <- function(x) {
  x / sum(x) * 100
}

# lets estimate actual means 
library(Seurat)
obj  = readRDS("~/Desktop/optimal_clustering_rna_only.rds")

sub_1_ = FindSubCluster(obj,1, graph.name='harmony.wsnn')
Idents(sub_1_) <- 'sub.cluster'
sub_1_ = subset(sub_1, sub.cluster %in%c('1_3', '1_2','1_0'))
sub_1_$Status = factor(sub_1_$Status, levels = c('NRM','M',"D",'E','NF','F'))

cells_ind = sub_1_@meta.data%>%
  group_by(individual)%>%
  summarize(n_cells = n())

cells_sub_ind = sub_1_@meta.data%>%
  group_by(individual, Status, sub.cluster)%>%
  summarize(n_cells_in = n())

cells_total = cells_ind%>%
  right_join(cells_sub_ind, by = 'individual')
cells_total$prop = cells_total$n_cells_in/cells_total$n_cells

mean_data =cells_total%>%group_by(Status, sub.cluster, individual)%>%
  summarize(mean = mean(prop))%>%subset(Status %in% c('M',
                                                      'D',
                                                      'F'))
mean_data$sub.cluster =factor(mean_data$sub.cluster, levels = c('1_3', '1_2', '1_0'))
ggplot(mean_data, aes(x=sub.cluster, y=mean, fill=Status)) +
  stat_summary(geom = 'bar', fun = 'mean', position = 'dodge') +
  geom_point(position = position_dodge(0.9))+
  theme_classic()


# Simulate conditions
male <- to_percent(simulate_model(0.2, 0.2, 0.4))
intermediate <- to_percent(simulate_model(0.25, 0.2, 0.4))
female <- to_percent(simulate_model(0.2, 0.2, 0.7))

# Combine into tidy dataframe
df <- rbind(
  data.frame(Condition="Male",         State=names(male),         Percent=male),
  data.frame(Condition="Intermediate", State=names(intermediate), Percent=intermediate),
  data.frame(Condition="Female",       State=names(female),       Percent=female)
)

df

df$Condition = factor(df$Condition, levels = c("Male", 'Intermediate', 'Female'))

ggplot(subset(df, State != 'N'), aes(x=State, y=Percent, fill=Condition)) +
  geom_bar(stat = 'identity', position = 'dodge') +
  theme_classic()


real_data = cells_total%>%group_by(Status, sub.cluster)%>%
  summarize(mean = mean(prop))%>%subset(Status %in% c('M',
                                                      'D',
                                                      'F'))%>%
  mutate(mean = mean*100, 
         half_mean = mean/2)

target_male <- c(P1=14.9, P2=13.8, P3= 21.4, N=(100)-(14.9+ 13.8+21.4))
target_female <- c(P1=14.5, P2=18.6, P3= 16.9, N=(100)-(14.5+ 18.6+16.9))
target_intermediate  <- c(P1=10.4, P2=15.9, P3= 23.6, N=(100)-(10.4+ 15.9+23.6))

loss_function <- function(par, target) {
  
  k1 <- par[1]
  k2 <- par[2]
  k3 <- par[3]
  
  # Penalize invalid k values
  if(any(par < 0) || any(par > 1)) return(1e6)
  
  sim <- to_percent(simulate_model(k1, k2, k3))
  
  sum((sim - target)^2)
}

fit_male <- optim(
  par = c(0.4, 0.4, 0.4),   # starting guesses
  fn  = loss_function,
  target = target_male,
  method = "L-BFGS-B",
  lower = c(0,0,0),
  upper = c(1,1,1)
)

fit_male$par

best_sim_m <- to_percent(simulate_model(
  fit_male$par[1],
  fit_male$par[2],
  fit_male$par[3]
))

rbind(Target=target_male, Simulated=best_sim)

# ok so we can fuck around with it to make these values true, so what?
# ok we are saying males have this baseline chance of each transition,
# so I should do with the other sexes now

fit_female <- optim(
  par = c(0.4, 0.4, 0.4),   # starting guesses
  fn  = loss_function,
  target = target_female,
  method = "L-BFGS-B",
  lower = c(0,0,0),
  upper = c(1,1,1)
)

fit_female$par

best_sim_f <- to_percent(simulate_model(
  fit_female$par[1],
  fit_female$par[2],
  fit_female$par[3]
))

rbind(Target=target_female, Simulated=best_sim)

fit_intermediate <- optim(
  par = c(0.4, 0.4, 0.4),   # starting guesses
  fn  = loss_function,
  target = target_intermediate,
  method = "L-BFGS-B",
  lower = c(0,0,0),
  upper = c(1,1,1)
)

fit_intermediate$par

best_sim_i <- to_percent(simulate_model(
  fit_intermediate$par[1],
  fit_intermediate$par[2],
  fit_intermediate$par[3]
))

rbind(Target=target_intermediate, Simulated=best_sim)

k_values = matrix(NA, 3, 3)
rownames(k_values) = c('M','I','F')
colnames(k_values) = c('k1','k2','k3')

k_values[,1] = fit_male$par
k_values[,2] = fit_intermediate$par
k_values[,3] = fit_female$par

k_values

for_plot = k_values%>%
  as.data.frame()%>%
  mutate(sex = rownames(k_values))%>%
  pivot_longer(cols=  c('k1','k2','k3'))
for_plot$sex = factor(for_plot$sex, levels = c('M','I','F'))

ggplot(for_plot, aes(x = name, y = value, fill = sex))+
  stat_summary(geom ='bar', position = 'dodge', fun = 'mean')+
  labs(x = 'Transition', y = 'Probability', fill = 'Phase')+
  scale_x_discrete(labels = c('1_3->1_2', '1_2->1_0', '1_0->N'))+
  theme_minimal()

# ok so here is the summary, the pattern of proportions 
#> is consistent with a model where intermediates and females 
#> have elevated 1_3 -> 1_2 transitions and 1_2->1_0 transitions, 
#> and there is a stepwise increase in the probability of a 1_0 -> n transition

for_plot$names = ifelse(for_plot$name == 'k1', '1_3->1_2',for_plot$name )
for_plot$names = ifelse(for_plot$name == 'k2', '1_2->1_0',for_plot$names )
for_plot$names = ifelse(for_plot$name == 'k3', '1_0->N',for_plot$names )
for_plot$names  = factor(for_plot$names , levels = c('1_3->1_2',
                                                     '1_2->1_0',
                                                     '1_0->N'))
probabilities =ggplot(for_plot, aes(x = sex, y = value, fill = sex)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ names) +
  theme_minimal() +
  labs(x = "Phase", y = "Probability")+
  theme(legend.position = 'none')

#ggsave(plot = probabilities,
 #      file = "2026_02_25 simulating neurogenesis transition_probabilities.svg",
  #     device = "svg",
   #    units = "in",
    #   width = 2,
     #  height = 1.5,
      # path = "Manuscript/Plots/RGC supplementary")

for_plot_2 = rbind(best_sim_m, 
                   best_sim_i, 
                   best_sim_f)%>%as.data.frame()
for_plot_2$sex = c('M', 'I', 'F')

for_plot_2_2 = for_plot_2%>%
  as.data.frame()%>%
  pivot_longer(cols=  c('P1','P2','P3', 'N'))%>%
  subset(name != 'N')
for_plot_2_2$sex = factor(for_plot_2_2$sex, levels = c('M','I','F'))
for_plot_2_2$name = factor(for_plot_2_2$name, levels = c('P1', 'P2', 'P3'))
for_plot_2_2$names = ifelse(for_plot_2_2$name == 'P1', '1_3', for_plot_2_2$name)
for_plot_2_2$names = ifelse(for_plot_2_2$name == 'P2', '1_2', for_plot_2_2$names)
for_plot_2_2$names = ifelse(for_plot_2_2$name == 'P3', '1_0', for_plot_2_2$names)
for_plot_2_2$names = factor(for_plot_2_2$names, levels = c('1_3', '1_2', '1_0'))

proportions =ggplot(for_plot_2_2, aes(x = sex, y = value/100, fill = sex)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ names, nrow= 1, ncol = 3) +
  theme_minimal() +
  labs(x = "Phase", y = "Simulated Proportion")+
  theme(legend.position = 'none')+
    scale_y_continuous(labels = scales::label_percent())
proportions

ggsave(plot = proportions,
       file = "2026_02_25 simulating neurogenesis proportions.svg",
       device = "svg",
       units = "in",
       width = 2,
       height = 1.5,
       path = "Manuscript/Plots/RGC supplementary")




               