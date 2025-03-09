
library(MASS)  # For mvrnorm
library(ggplot2)  # For plotting

data1 <- data.frame(PC1 = rpois(100, 2),
                   PC2 = rpois(100, 2),
                   Status = sample(c('M','D','F'), replace =T, size = 100))


data = data1[,1:2]
# Estimate parameters from data (two columns indicating coordinates)
mu_hat <- colMeans(data)  # Estimated mean
Sigma_hat <- cov(data)  # Estimated covariance matrix

# Step 2: Simulate data from the estimated bivariate normal distribution
simulated_data <- mvrnorm(n = 10000, mu = mu_hat, Sigma = Sigma_hat) ### here, this needs to be males and females only

# Step 3: Plot the simulated data
df_sim <- as.data.frame(simulated_data)
colnames(df_sim) <- c("X1", "X2")

ggplot(df_sim, aes(x = X1, y = X2)) +
  geom_point(color = "blue", alpha = 0.25) +
  labs(title = "Simulated Bivariate Normal Data") +
  theme_minimal()+
  geom_point(aes(x = 2, y = 5), color = 'red')

# Step 4: Calculate p-value for a given coordinate (a, b)
a <- data1[1,1]  # Example coordinate
b <- data1[1,2]

# Compute Mahalanobis distance
point <- c(a, b)
mahalanobis_dist <- mahalanobis(point, mu_hat, Sigma_hat)

# Compute p-value using chi-square distribution with 2 degrees of freedom
p_value <- 1 - pchisq(mahalanobis_dist, df = 2)
p_value


#so I think here I input the whole set of males and females and calcuate if the mean dominant is significant?



