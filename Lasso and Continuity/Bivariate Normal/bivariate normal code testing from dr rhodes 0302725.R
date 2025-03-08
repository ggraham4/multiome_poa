
library(MASS)  # For mvrnorm
library(ggplot2)  # For plotting


# Estimate parameters from data (two columns indicating coordinates)
mu_hat <- colMeans(data)  # Estimated mean
Sigma_hat <- cov(data)  # Estimated covariance matrix

# Step 2: Simulate data from the estimated bivariate normal distribution
simulated_data <- mvrnorm(n = 10000, mu = mu_hat, Sigma = Sigma_hat)

# Step 3: Plot the simulated data
df_sim <- as.data.frame(simulated_data)
colnames(df_sim) <- c("X1", "X2")

ggplot(df_sim, aes(x = X1, y = X2)) +
  geom_point(color = "blue", alpha = 0.5) +
  labs(title = "Simulated Bivariate Normal Data") +
  theme_minimal()

# Step 4: Calculate p-value for a given coordinate (a, b)
a <- 2  # Example coordinate
b <- 5

# Compute Mahalanobis distance
point <- c(a, b)
mahalanobis_dist <- mahalanobis(point, mu_hat, Sigma_hat)

# Compute p-value using chi-square distribution with 2 degrees of freedom
p_value <- 1 - pchisq(mahalanobis_dist, df = 2)
p_value