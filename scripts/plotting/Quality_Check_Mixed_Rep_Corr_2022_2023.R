#!/usr/bin/env Rscript

#SBATCH --ntasks=1
#SBATCH --mem=30G
#SBATCH --time=02:00:00
#SBATCH --output=/rhome/jmarz001/bigdata/Ag-Competition/scripts/CCII_Ag_Comp_Graphs_5.stdout
#SBATCH -p koeniglab

library(tidyverse)
library(ggpubr)
library(car)
library(gridExtra)
library(dunn.test)

Replicate_corr_tbl_Mixed <- Replicate_corr_tbl %>% filter(Condition == "mixed")

### Correlation of Replicates (Total Seed Weight) 

g2 <- ggplot(Replicate_corr_tbl_Mixed, aes(x=TOTAL_WEIGHT, y=`TOTAL_WEIGHT_2`, add = "reg.line")) +
  geom_point(alpha = .5) +
  geom_abline(slope =1, intercept= 0, color = "red") +
  stat_cor(label.y = 235, label.x = 20) +
  geom_smooth(method = "lm") +
  labs(x = "Rep 1",
       y= "Rep 2",
       title = "Total Weight 2022-20223") 
ggsave("scripts/plotting/Correlation_of_Mixed_Rep_2022_2023_TW.png")

res <- Replicate_corr_tbl %>% summarise(residuals = resid(lm(`TOTAL_WEIGHT_2` ~ TOTAL_WEIGHT)))
Replicate_corr_tbl$Residuals <- res$residuals

### Correlation of Replicates (Flowering Time)

g4 <- ggplot(Replicate_corr_tbl_Mixed, aes(FT, FT_2)) +
  geom_jitter(alpha = .5) +
  geom_smooth(method = "lm") +
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
  geom_smooth(method = "lm") +
  stat_cor(label.x = 100, label.y = 40000) +
  geom_abline(slope = 1, intercept = 0, col = "red") +
  labs(x = "Rep 1",
       y = "Rep 2",
       title = "Fitness 2022-2023")
ggsave("scripts/plotting/Correlation_of_Mixed_Rep_2022_2023_Fitness.png")

res <- Replicate_corr_tbl %>% summarise(residuals = resid(lm(Fitness_2 ~ Fitness)))
Replicate_corr_tbl$Residuals <- res$residuals
