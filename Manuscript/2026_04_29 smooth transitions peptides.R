library(ggplot2)

N <- 500
xs <- seq(0, 1, length.out = N)

smooth_pts <- function(pts_x, pts_y, xs, bw = 25) {
  sapply(xs, function(x) {
    w <- exp(-bw * (x - pts_x)^2)
    sum(w * pts_y) / sum(w)
  })
}

sigmoid <- function(x, center, slope) 1 / (1 + exp(-slope * (x - center)))

phase_labels <- c(0, 1/3, 2/3, 1)

newd <- data.frame(x = xs)
newd$cyp19a1b  <- smooth_pts(phase_labels, c(1, 0, 0, 1), xs)
newd$tacr3a    <- smooth_pts(phase_labels, c(0.5, 1, 0, 0), xs)
newd$npy7r     <- smooth_pts(phase_labels, c(1, 0.25, 0, 0), xs)
newd$drd3      <- smooth_pts(phase_labels, c(0, 0, 0, 1), xs)
newd$cckb      <- smooth_pts(phase_labels, c(0, 0, 1, 1), xs)
newd$gnrh1     <- 0.5
newd$E2        <-smooth_pts(phase_labels, c(0, 0, 1, 1), xs)
newd$KT11      <- pmax(0, 1 - sigmoid(xs, 0.45, 20))
newd$LH        <- pmax(0, 1 - sigmoid(xs, 0.45, 20))
newd$FSH       <- pmin(1, sigmoid(xs, 0.65, 18))

colors <- c(
  "cyp19a1b" = "#882E72", #9
  "tacr3a"   = "#5289C7", #12
  "npy7r"    = "#4EB265",#15
  "drd3"     = "#F7F056",#18
  "cckb"     = "#EE8026", #24
  "gnrh1"    = "#1965B0", # 26
  "E2"       = "black",#21
  "11KT"     = "#CAE0AB",#17
  "LH"       = "#7BAFDE",# 14
  "FSH"      = "#E8601C" # 10 
)

linetypes <- c(
  "cyp19a1b" = "solid", "tacr3a" = "solid", "npy7r" = "solid",
  "drd3" = "dashed", "cckb" = "solid",
  "gnrh1" = "solid",
  "E2" = "dashed", "11KT" = "solid", "LH" = "dotted", "FSH" = "dotted"
)

p <- ggplot(newd, aes(x = x)) +
  geom_line(aes(y = cyp19a1b, color = "cyp19a1b", linetype = "cyp19a1b"), linewidth = 1.8) +
  geom_line(aes(y = tacr3a,   color = "tacr3a",   linetype = "tacr3a"),   linewidth = 1.8) +
  geom_line(aes(y = npy7r,    color = "npy7r",    linetype = "npy7r"),    linewidth = 1.8) +
  geom_line(aes(y = drd3,     color = "drd3",     linetype = "drd3"),     linewidth = 1.8) +
  geom_line(aes(y = cckb,     color = "cckb",     linetype = "cckb"),     linewidth = 1.8) +
  geom_line(aes(y = gnrh1,    color = "gnrh1",    linetype = "gnrh1"),    linewidth = 1.8) +
  geom_line(aes(y = E2,       color = "E2",       linetype = "E2"),       linewidth = 1.8) +
  geom_line(aes(y = KT11,     color = "11KT",     linetype = "11KT"),     linewidth = 1.8) +
  geom_line(aes(y = LH,       color = "LH",       linetype = "LH"),       linewidth = 1.8) +
  geom_line(aes(y = FSH,      color = "FSH",      linetype = "FSH"),      linewidth = 1.8) +
  scale_color_manual(values = colors) +
  scale_linetype_manual(values = linetypes) +
  scale_y_continuous(breaks = c(0, 1), labels = c("Min", "Max"), limits = c(0, 1)) +
  scale_x_continuous(
    breaks = c(0, 1/3, 2/3, 1),
    labels = c("M", "I", "LI", "F")
  ) +
  theme_classic() +
  labs(color = "Measurement", linetype = "Measurement") +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x  = element_text(size = 10),
    legend.text  = element_text(face = "italic")
  )

p

ggsave(plot = p,
       file = "smooth_transitions_post_neuroendocrine_phases.svg",
       device = "svg",
       units = "in",
       width = 6.7,
       height = 2.2, 
       path = "Manuscript/Plots/Manuscript v1.2")
