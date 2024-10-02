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

g1 <- ggplot(Rep_2021_2022_Single, aes(total_seed_mass_g, total_seed_mass_g_2, add = "reg.line")) +
  geom_jitter(alpha = .5) +
  stat_cor(label.y = 175) +
  geom_abline(slope = 1, intercept = 0, color = 'red') +
  geom_smooth(method = 'lm') +
  labs(x = "Rep 1",
       y = "Rep 2",
       title = "Total Weight 2021-2022")

### Correlation for FT

g3 <- ggplot(Rep_2021_2022_Single, aes(Flowering_Date, Flowering_Date_2, add = "reg.line")) +
  geom_jitter(alpha = .5) +
  stat_cor(label.y = 120) +
  geom_abline(slope = 1, intercept = 0, color = 'red') +
  geom_smooth(method = 'lm') +
  labs(x = "Rep 1",
       y = "Rep 2",
       title = "Flowering Time 2021-2022")

### Correlation for Fec

g5 <- ggplot(Rep_2021_2022_Single, aes(Fecundity, Fecundity_2, add = "reg.line")) +
  geom_jitter(alpha = .5) +
  stat_cor(label.y = 2000, label.x = 6000) +
  geom_abline(slope = 1, intercept = 0, color = 'red') +
  geom_smooth(method = 'lm') +
  labs(x = "Rep 1",
       y = "Rep 2",
       title = 'Fecundity 2021-2022')

### Correlation for Fit

g7 <- ggplot(Rep_2021_2022_Single, aes(Fitness, Fitness_2, add = "reg.line")) +
  geom_jitter(alpha = .5) +
  stat_cor(label.y = 25000, label.x = 50000) +
  geom_abline(slope = 1, intercept = 0, color = 'red') +
  geom_smooth(method = 'lm') +
  labs(x = "Rep 1",
       y = "Rep 2",
       title = 'Fitness 2021-2022')

h <- arrangeGrob(g1, g2, g3, g4, g5, g6, g7, g8, ncol = 2, nrow = 4, top = "Correlation of Single Replicates For All Phenotypes")
ggsave("scripts/plotting/Quality_Check_Single_Corr_All_Years.pdf", h, width = 12, height = 14)
