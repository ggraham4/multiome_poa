library(ggplot2)
library(dplyr)
library(zoo) # For the rolling window function

# 1. Sort the data by your Latent Variable (the X-axis of your hill)
measures <- measures %>% arrange(SexState)

measures2 = read.csv('Measures/coalesced_data.csv')

# 2. Calculate a rolling Standard Deviation
# 'width' is the number of fish to include in each window. 
# Adjust this based on your sample size (e.g., 5-10% of total N).
window_size <- 0.1*22

measures <- measures %>%
  mutate(
    # Rolling SD of the Latent Variable itself
    rolling_sd = rollapply(SexState, width = window_size, FUN = sd, fill = NA),
    # Or Rolling SD of a specific driver like Log_11KT
    rolling_sd_11KT = rollapply(Log_11KT, width = window_size, FUN = sd, fill = NA),
    
        rolling_sd_test = rollapply(Percent_Testicular, width = window_size, FUN = sd, fill = NA),
        rolling_sd_mass = rollapply(mass_final, width = window_size, FUN = sd, fill = NA)

  )

# 3. Create the Plot
ggplot(subset(measures, Status %in% c('M','D','E','F')), aes(x = SexState, y = rolling_sd)) +
  # Add a smooth line to see the "shape" of the variance
  geom_smooth(method = "loess", color = "firebrick", se = FALSE) + 
  # Add the raw points to see the density
  geom_point(alpha = 0.3, aes(color =Status)) +
  # Label based on your status to see where M, D, and F sit
  geom_rug(aes(color = Status), sides = "b") + 
  labs(
    title = "Variance across the Sex Change Landscape",
    x = "Latent Sex State (Male -> Female)",
    y = paste0("Rolling SD (Window Size = ", window_size, ")")
  ) +
  theme_minimal()

ggplot(subset(measures, Status %in% c('M','D','E','F')), aes(x = SexState, y = rolling_sd_11KT)) +
  # Add a smooth line to see the "shape" of the variance
  geom_smooth(method = "loess", color = "firebrick", se = FALSE) + 
  # Add the raw points to see the density
  geom_point(alpha = 0.3, aes(color =Status)) +
  # Label based on your status to see where M, D, and F sit
  geom_rug(aes(color = Status), sides = "b") + 
  labs(
    title = "11KT Variance across the Sex Change Landscape",
    x = "Latent Sex State (Male -> Female)",
    y = paste0("Rolling SD (Window Size = ", window_size, ")")
  ) +
  theme_minimal()

ggplot(subset(measures, Status %in% c('M','D','E','F')), aes(x = SexState, y = rolling_sd_test)) +
  # Add a smooth line to see the "shape" of the variance
  geom_smooth(method = "loess", color = "firebrick", se = FALSE) + 
  # Add the raw points to see the density
  geom_point(alpha = 0.3, aes(color =Status)) +
  # Label based on your status to see where M, D, and F sit
  geom_rug(aes(color = Status), sides = "b") + 
  labs(
    title = "Testicular Variance across the Sex Change Landscape",
    x = "Latent Sex State (Male -> Female)",
    y = paste0("Rolling SD (Window Size = ", window_size, ")")
  ) +
  theme_minimal()

ggplot(subset(measures, Status %in% c('M','D','E','F')), aes(x = SexState, y = rolling_sd_mass)) +
  # Add a smooth line to see the "shape" of the variance
  geom_smooth(method = "loess", color = "firebrick", se = FALSE) + 
  # Add the raw points to see the density
  geom_point(alpha = 0.3, aes(color =Status)) +
  # Label based on your status to see where M, D, and F sit
  geom_rug(aes(color = Status), sides = "b") + 
  labs(
    title = "Mass Variance across the Sex Change Landscape",
    x = "Latent Sex State (Male -> Female)",
    y = paste0("Rolling SD (Window Size = ", window_size, ")")
  ) +
  theme_minimal()



#####
measures2 = read.csv('Measures/2026_01_08 Coalesced sex state.csv')

window_size = 10

measures2 <- measures2 %>%
  mutate(
    rolling_sd = rollapply(SexState, width = window_size, FUN = sd, fill = NA),
    rolling_sd_11KT = rollapply(Log_11KT, width = window_size, FUN = sd, fill = NA),
    rolling_sd_test = rollapply(Percent_Testicular, width = window_size, FUN = sd, fill = NA),
    rolling_sd_mass = rollapply(mass_final, width = window_size, FUN = sd, fill = NA)
  )

measures2$rolling_sd_test[measures2$Phase=='F']=0


ggplot(subset(measures2, Phase %in% c('M','D','E','NF','F')), aes(x = SexState, y = rolling_sd)) +
  geom_smooth(method = "loess", color = "firebrick", se = FALSE) + 
  geom_point(alpha = 0.3, aes(color =Phase)) +
  geom_rug(aes(color = Phase), sides = "b") + 
  theme_minimal()

ggplot(subset(measures2, Phase %in% c('M','D','E','NF','F')), aes(x = SexState, y = rolling_sd_11KT)) +
  geom_smooth(method = "loess", color = "firebrick", se = FALSE) + 
  geom_point(alpha = 0.3, aes(color =Phase)) +
  geom_rug(aes(color = Phase), sides = "b") + 
  theme_minimal()

ggplot(subset(measures2, Phase %in% c('M','D','E','NF','F')), aes(x = SexState, y = rolling_sd_test)) +
  geom_smooth(method = "loess", color = "firebrick", se = FALSE) + 
  geom_point(alpha = 0.3, aes(color =Phase)) +
  geom_rug(aes(color = Phase), sides = "b") + 
  theme_minimal()

ggplot(subset(measures2, Phase %in% c('M','D','E','NF','F')), aes(x = SexState, y = rolling_sd_mass)) +
  geom_smooth(method = "loess", color = "firebrick", se = FALSE) + 
  geom_point(alpha = 0.3, aes(color =Phase)) +
  geom_rug(aes(color = Phase), sides = "b") + 
  theme_minimal()



#does not agree

