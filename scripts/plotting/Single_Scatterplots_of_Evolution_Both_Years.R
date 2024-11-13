#!/usr/bin/env Rscript
#SBATCH --mem=30G
#SBATCH --time=02:00:00
#SBATCH --output=/rhome/jmarz001/bigdata/Ag-Competition/scripts/Single_Scatterplots_of_Evolution_Both_Years.stdout
#SBATCH -p koeniglab

library(tidyverse)
library(ggpubr)
library(car)
library(gridExtra)

setwd("/bigdata/koeniglab/jmarz001/Ag-Competition/data")

df <- PHENO_FULL_AVERAGE %>% filter(Condition == "single")

#fitness = (Plot Germination * Fecundity)
### 02a_Single_Fitness_over_Generation.R

df$Generation <- as.numeric(df$Generation)
a <- ggplot(df, aes(x= Generation, y = ABS_FITNESS)) +
  geom_jitter(alpha = .5) +
  geom_smooth(method = "lm") +
  geom_hline(aes(yintercept = 0), color = "red") +
  geom_boxplot(aes(Generation, ABS_FITNESS, group = Generation), width = 1.5, alpha = .5) +
  stat_regline_equation(label.y = 2.5) +
  labs(x = "Generation",
       y = "Absolute Fitness",
       title = "Evolution of Centered Absolute Fitness")
ggsave("/bigdata/koeniglab/jmarz001/Ag-Competition/results/02a_Generational_Change_in_Single_Fitness.png")

### 02ai_Single_Fecundity_over_Generations.R

b <- ggplot(df, aes(Generation, FECUNDITY, add = "reg.line")) +
  geom_jitter(alpha = .5) +
  geom_smooth(method = lm) +
  geom_hline(aes(yintercept = 0), color = "red") +
  geom_boxplot(aes(Generation, FECUNDITY, group = Generation), width = 1.5, alpha = .5) +
  stat_regline_equation(label.y = 3) +
  labs(x = "Generation",
       y = "Average Fecundity",
       title = "Evolution of Centered Fecundity")
ggsave("/bigdata/koeniglab/jmarz001/Ag-Competition/results/02ai_Generational_Change_in_Single_Fecundity.png")

### 02aii_Single_FT_over_Generations.R

c <- ggplot(df, aes(Generation, FT, add = "reg.line")) +
  geom_jitter(alpha = .5) +
  geom_smooth(method = lm) +
  stat_regline_equation(label.x = 40) +
  geom_boxplot(aes(Generation, FT, group = Generation), width = 1.5, alpha = .5)+
  labs(x = "Generation",
       y = "Average Flowering Time (Days)",
       title = "Evolution of Average Flowering Time")
ggsave("/bigdata/koeniglab/jmarz001/Ag-Competition/results/02aii_Generational_Change_in_Single_FT.png")

### 02aiii_Single_TW_over_Generations.R

d <- ggplot(df, aes(Generation, TOTAL_MASS, add = "reg.line")) +
  geom_jitter(alpha = .5) +
  geom_smooth(method = lm) +
  geom_boxplot(aes(Generation,TOTAL_MASS, group = Generation), width = 1.5, alpha = .5)+
  stat_regline_equation(label.y = 135) +
  labs(x = "Generation",
       y = "Average Total Seed Weight (g)",
       title = "Evolution of Average Total Weight (grams)")
ggsave("/bigdata/koeniglab/jmarz001/Ag-Competition/results/02aiii_Generational_Change_in_Single_TW.png")

### 02aiiii_Single_100SW_over_Generations.R

e <- ggplot(df, aes(Generation, SEED_WEIGHT_100, add = "reg.line")) +
  geom_jitter(alpha = .5) +
  geom_smooth(method = "lm") +
  geom_boxplot(aes(Generation, SEED_WEIGHT_100, group = Generation), width = 1.5, alpha = .5) +
  stat_regline_equation(label.y = 5.9) +
  labs(y = "Average 100 Seed Weight (grams)",
       title = "Evolution of Average 100 Seed Weight (grams)")
ggsave("/bigdata/koeniglab/jmarz001/Ag-Competition/results/02aiiii_Generational_Change_in_Single_100SW.png")







### 02ci_Mixed_Relative_Fitness_to_Atlas.R

mf <- ggplot(Mixed, aes(Generation, Fitness, add = "reg.line")) +
  geom_jitter(alpha = .4) +
  ylim(0, 40000) +
  stat_regline_equation() +
  geom_hline(yintercept = 21347.22, color = "red") +
  geom_smooth(method = "lm") +
  labs(title = "Relative Fitness of Atlas Compared to Mixed Condition")
ggsave("scripts/plotting/02ci_Mixed_Relative_Fitness_to_Atlas.png")

jf <- ggplot(Averaged_Full, aes(Generation, Fitness, add = "reg.line", color = Condition)) +
  geom_jitter(alpha = .6) +
  stat_regline_equation() +
  geom_smooth(method = "lm" ) +
  geom_hline(yintercept = 21347.22, color = "red") +
  labs(title = "Relative Fitness of Atlas Compared to Both Conditions")
ggsave("scripts/plotting/02cii_Relative_Fitness_to_Atlas_Combined.png", width = 12, height = 7)

### 02cii_Relative_Fitness_to_Atlas_Combined.R

g <- arrangeGrob(sf, mf, jf, nrow = 2)
ggsave("scripts/plotting/02cii_Relative_Fitness_to_Atlas_Combined.png", g, width = 12, height = 7)
