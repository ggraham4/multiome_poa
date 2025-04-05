
### what is the range of gamma distributions that would be necessary for p != 0.05?

#first, simulate normal
normal_distribution <- rnorm(10000)
hist(normal_distribution)

empirical_0.05 <- quantile(normal_distribution, .95)


data = data.frame()
for(i in seq(2,11, by =0.001)){
  print(i)

    gamma_distribution <- rgamma(10000, i)
  
    gamma_p_value_pos = sum(gamma_distribution>=empirical_0.05 )/length(gamma_distribution)
    gamma_p_value_neg = sum(gamma_distribution<=empirical_0.05 )/length(gamma_distribution)
    
    newd <- data.frame(shape = i,
                       gamma_p_value_pos =gamma_p_value_pos,
                       gamma_p_value_neg=gamma_p_value_neg)
    data = rbind(newd, data)
}

data_subset = data[(data$gamma_p_value_pos<=0.05 | data$gamma_p_value_neg<=0.05) &
       (data$gamma_p_value_pos !=1 &data$gamma_p_value_pos !=0)&
       (data$gamma_p_value_neg !=1 &data$gamma_p_value_neg !=0 ),]

hist(data_subset$shape,20)
