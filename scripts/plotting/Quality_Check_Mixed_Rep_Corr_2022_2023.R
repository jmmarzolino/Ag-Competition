library(tidyverse)
library(readr)
library(dplyr)
library(ggpubr)
library(ggplot2)

library(tidyr)
library(car)
library(gridExtra)
library(dunn.test)

Replicate_corr_tbl_Mixed <- Replicate_corr_tbl %>% filter(Condition == "mixed")

### Correlation of Replicates (Total Seed Weight) 

g2 <- ggplot(Replicate_corr_tbl_Mixed, aes(x=`Brown Bag Weight`, y=`Brown Bag Weight_2`, add = "reg.line")) +
  geom_point(alpha = .5) +
  geom_abline(slope =1, intercept= 0, color = "red") +
  stat_cor(label.y = 235, label.x = 20) +
  geom_smooth(method = "lm") +
  labs(x = "Rep 1",
       y= "Rep 2",
       title = "Total Weight 2022-20223") 
ggsave("scripts/plotting/Correlation_of_Mixed_Rep_2022_2023_TW.png")

res <- Replicate_corr_tbl %>% summarise(residuals = resid(lm(`Brown Bag Weight_2` ~ `Brown Bag Weight`)))
Replicate_corr_tbl$Residuals <- res$residuals

### Correlation of Replicates (Flowering Time)

g4 <- ggplot(Replicate_corr_tbl_Mixed, aes(FT_DAYS, FT_DAYS_2)) +
  geom_jitter(alpha = .5) +
  geom_smooth(method = 'lm') +
  stat_cor(label.x = 100, label.y = 135) +
  geom_abline(slope = 1, intercept = 0, col = "red") +
  labs(x = "Rep 1",
       y = "Rep 2",
       title = "Flowering Time 2022-2023")
ggsave("scripts/plotting/Correlation_of_Mixed_Rep_2022_2023_FT.png")

### Correlation of Replicates (Fecundity)

g6 <- ggplot(Replicate_corr_tbl_Mixed, aes(x=Fecundity, y= Fecundity_2)) +
  geom_point(alpha = .5) +
  stat_cor(label.y = 4200, label.x = 20) +
  geom_abline(slope =1, intercept= 0, color = "red") +
  geom_smooth(method = "lm") +
  labs(x = "Rep 1",
       y= "Rep 2",
       title = "Fecundity 2022-2023") 
ggsave("scripts/plotting/Correlation_of_Mixed_Rep_2022_2023_Fecundity.png")

res <- Replicate_corr_tbl %>% summarise(residuals = resid(lm(Fecundity_2 ~ Fecundity)))
Replicate_corr_tbl$Residuals <- res$residuals

### Correlation of Replicates (Fitness)

g8 <- ggplot(Replicate_corr_tbl_Mixed, aes(Fitness, Fitness_2)) +
  geom_jitter(alpha = .5) +
  geom_smooth(method = 'lm') +
  stat_cor(label.x = 100, label.y = 40000) +
  geom_abline(slope = 1, intercept = 0, col = "red") +
  labs(x = "Rep 1",
       y = "Rep 2",
       title = "Fitness 2022-2023")
ggsave("scripts/plotting/Correlation_of_Mixed_Rep_2022_2023_Fitness.png")

res <- Replicate_corr_tbl %>% summarise(residuals = resid(lm(Fitness_2 ~ Fitness)))
Replicate_corr_tbl$Residuals <- res$residuals
