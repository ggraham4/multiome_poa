library(Polychrome)
# Create a 40-color palette
set.seed(1) # For reproducibility
p40 <- createPalette(8, c("#FF0000", "#00FF00", "#0000FF"), range = c(30, 80))

# Display the swatches
swatch(p40)
names(p40) <- NULL

newd = data.frame(Sex = seq(0,1, by =0.01))

newd$BodyMass = newd$Sex
newd$PercentTesticular <- 1 - (log(1 + 800 * newd$Sex) / log(801))
newd$Aggression <-  (log(1 + 8000 * newd$Sex) / log(8001))
newd$Parenting <- 1 / (1 + exp(30 * (newd$Sex - 0.65)))
newd$KT <-  1 / (1 + exp(30 * (newd$Sex - 0.65)))
newd$Vitellogenesis <- (newd$Sex * 4) -3
newd$E2 <- 1 / (1 + exp(-30 * (newd$Sex - 0.75)))
x <- pmin(newd$Sex, 0.75) / 0.75
newd$PercentTesticular <- 1 - (log(1 + 8 * x) / log(9))

newd$Vitellogenesis = ifelse(newd$Vitellogenesis < 0, NA, newd$Vitellogenesis )
newd$Vitellogenesis = ifelse(newd$Vitellogenesis < 0, NA, newd$Vitellogenesis )

newd$Neurogenesis <- ifelse(
  newd$Sex <= 0.75,
  (newd$Sex / 0.75)^2 * 0.5,
  0.5 + (1 - exp(-12 * (newd$Sex - 0.75))) * 0.5
)



p = ggplot(newd, aes(x = Sex, y = BodyMass))+
  ylim(0,1)+
  geom_line( linewidth = 2, aes(color = 'Body Mass'))+
  geom_line(aes(y = PercentTesticular,color = 'Percent Testicular'),linewidth = 2)+
  geom_line(aes(y = Aggression,color = 'Aggression'),linewidth = 2)+
    geom_line(aes(y = Parenting,color = 'Parental Care'),linewidth = 2)+
    geom_line(aes(y = KT,color = '11KT'),linewidth = 2, linetype = 2)+
        geom_line(aes(y = Neurogenesis,color = 'Neurogenesis'),linewidth = 2)+
      geom_line(aes(y = E2,color = 'E2'),linewidth = 2)+
      geom_line(aes(y = Vitellogenesis,color = 'Vitellogenesis'),linewidth = 2)+
  scale_color_manual(values = c(
    "#000000",  # black
"#D55E00",  # vermillion
"#7F7F7F",  # dark gray
"#0072B2",  # blue
"#CC79A7",  # purple-pink
 "#009E73",  # green
"#E69F00",# orange
"#F0E442"   # yellow
))+
  theme_classic()+
  labs(x = 'Phase', y = 'Arbitrary Value', color = 'Measurement')+
  scale_y_continuous(
  breaks = c(0, 1),
  labels = c("Min", "Max")
) +
scale_x_continuous(
  breaks = c(0, 0.5, 0.75, 1),
  labels = c("M", "I", "LI", "F")
)+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

p

ggsave(plot = p,
       file = "smooth_transition.svg",
       device = "svg",
       units = "in",
       width = 8,
       height = 3,
       path = "Manuscript/Plots/Fig.1")

