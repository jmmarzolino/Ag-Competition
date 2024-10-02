library(tidyverse)
library(readr)
library(dplyr)
library(ggpubr)
library(ggplot2)

library(tidyr)
library(car)
library(gridExtra)
library(dunn.test)

### Correlation For TW

g1 <- ggplot(Rep_2021_2022_Mixed, aes(total_seed_mass_g, total_seed_mass_g_2, add = "reg.line")) +
  geom_jitter(alpha = .5) +
  stat_cor(label.y = 175) +
  geom_abline(slope = 1, intercept = 0, color = 'red') +
  geom_smooth(method = 'lm') +
  labs(x = "Rep 1",
       y = "Rep 2",
       title = "Total Weight 2021-2022")

### Correlation for FT

g3 <- ggplot(Rep_2021_2022_Mixed, aes(Flowering_Date, Flowering_Date_2, add = "reg.line")) +
  geom_jitter(alpha = .5) +
  stat_cor(label.y = 120) +
  geom_abline(slope = 1, intercept = 0, color = 'red') +
  geom_smooth(method = 'lm') +
  labs(x = "Rep 1",
       y = "Rep 2",
       title = "Flowering Time 2021-2022")

### Correlation for Fec

g5 <- ggplot(Rep_2021_2022_Mixed, aes(Fecundity, Fecundity_2, add = "reg.line")) +
  geom_jitter(alpha = .5) +
  stat_cor(label.y = 3500, label.x = 1000) +
  geom_abline(slope = 1, intercept = 0, color = 'red') +
  geom_smooth(method = 'lm') +
  labs(x = "Rep 1",
       y = "Rep 2",
       title = "Fecundity 2021-2022")

### Correlation for Fit

g7 <- ggplot(Rep_2021_2022_Mixed, aes(Fitness, Fitness_2, add = "reg.line")) +
  geom_jitter(alpha = .5) +
  stat_cor(label.y = 35000) +
  geom_abline(slope = 1, intercept = 0, color = 'red') +
  geom_smooth(method = 'lm') +
  labs(x = "Rep 1",
       y = "Rep 2",
       title = "Fitness 2021-2022")

h <- arrangeGrob(g1, g2, g3, g4, g5, g6, g7, g8, nrow = 4, ncol = 2, top = "Correlation of Mixed Replicates For All Phenotypes")
ggsave("scripts/plotting/Quality_Check_Mixed_Corr_All_Years.png", h, width = 14, height = 16)

