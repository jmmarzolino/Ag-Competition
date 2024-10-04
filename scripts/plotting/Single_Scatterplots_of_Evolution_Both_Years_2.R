#!/usr/bin/env Rscript

#SBATCH --ntasks=1
#SBATCH --mem=30G
#SBATCH --time=02:00:00
#SBATCH --output=/rhome/jmarz001/bigdata/Ag-Competition/competition1.stdout
#SBATCH -p koeniglab

library(tidyverse)
library(readr)
library(dplyr)
library(ggpubr)
library(ggplot2)

library(tidyr)
library(car)
library(gridExtra)

setwd("/bigdata/koeniglab/jmarz001/Ag-Competition/data")

# Single Scatterplots 2021-2022

Single_2021_2022 <- PHENO_FULL_AVERAGE %>% filter(Condition == "single" & Exp_year == 2022)

#fitness = (Plot Germination * Fecundity)
### 02a_Single_Fitness_over_Generation.R

Single_2021_2022$Generation <- as.numeric(Single_2021_2022$Generation)
a <- ggplot(Single_2021_2022, aes(x= Generation, y = ABS_FITNESS)) +
  geom_jitter(alpha = .5) +
  geom_smooth(method = "lm") +
  geom_hline(aes(yintercept = 0), color = "red") +
  geom_boxplot(aes(Generation, ABS_FITNESS, group = Generation), width = 1.5, alpha = .5) +
  stat_regline_equation(label.y = 2.5) +
  labs(x = "Generation",
       y = "Absolute Fitness",
       title = "Evolution of Centered Absolute Fitness 2021-2022")
ggsave("/bigdata/koeniglab/jmarz001/Ag-Competition/results/02a_Generational_Change_in_Single_Fitness_2021_2022.png")

### 02ai_Single_Fecundity_over_Generations.R

b <- ggplot(Single_2021_2022, aes(Generation, FECUNDITY, add = "reg.line")) +
  geom_jitter(alpha = .5) +
  geom_smooth(method = lm) +
  geom_hline(aes(yintercept = 0), color = "red") +
  geom_boxplot(aes(Generation, FECUNDITY, group = Generation), width = 1.5, alpha = .5) +
  stat_regline_equation(label.y = 3) +
  labs(x = "Generation",
       y = "Average Fecundity",
       title = "Evolution of Centered Fecundity 2021-2022")
ggsave("/bigdata/koeniglab/jmarz001/Ag-Competition/results/02ai_Generational_Change_in_Single_Fecundity_2021_2022.png")

### 02aii_Single_FT_over_Generations.R

c <- ggplot(Single_2021_2022, aes(Generation, FT, add = "reg.line")) +
  geom_jitter(alpha = .5) +
  geom_smooth(method = lm) +
  stat_regline_equation(label.x = 40) +
  geom_boxplot(aes(Generation, FT, group = Generation), width = 1.5, alpha = .5)+
  labs(x = "Generation",
       y = "Average Flowering Time (Days)",
       title = "Evolution of Average Flowering Time 2021-2022")
ggsave("/bigdata/koeniglab/jmarz001/Ag-Competition/results/02aii_Generational_Change_in_Single_FT_2021_2022.png")

### 02aiii_Single_TW_over_Generations.R

d <- ggplot(Single_2021_2022, aes(Generation, TOTAL_MASS, add = "reg.line")) +
  geom_jitter(alpha = .5) +
  geom_smooth(method = lm) +
  geom_boxplot(aes(Generation,TOTAL_MASS, group = Generation), width = 1.5, alpha = .5)+
  stat_regline_equation(label.y = 135) +
  labs(x = "Generation",
       y = "Average Total Seed Weight (g)",
       title = "Evolution of Average Total Weight (grams) 2021-2022")
ggsave("/bigdata/koeniglab/jmarz001/Ag-Competition/results/02aiii_Generational_Change_in_Single_TW_2021_2022.png")

### 02aiiii_Single_100SW_over_Generations.R

e <- ggplot(Single_2021_2022, aes(Generation, SEED_WEIGHT_100, add = "reg.line")) +
  geom_jitter(alpha = .5) +
  geom_smooth(method = "lm") +
  geom_boxplot(aes(Generation, SEED_WEIGHT_100, group = Generation), width = 1.5, alpha = .5) +
  stat_regline_equation(label.y = 5.9) +
  labs(y = "Average 100 Seed Weight (grams)",
       title = "Evolution of Average 100 Seed Weight (grams) 2021-2022")
ggsave("/bigdata/koeniglab/jmarz001/Ag-Competition/results/02aiiii_Generational_Change_in_Single_100SW_2021_2022.png")

###################################################

# Single Scatterplots 2022-2023

Single_2022_2023 <- PHENO_FULL_AVERAGE %>% filter(Condition == "single" & Exp_year == 2023)

#fitness = (Plot Germination * Fecundity)
### 02b_Single_Fitness_over_Generation.R

Single_2022_2023$Generation <- as.numeric(Single_2022_2023$Generation)
f <- ggplot(Single_2022_2023, aes(x= Generation, y = ABS_FITNESS)) +
  geom_jitter(alpha = .5) +
  geom_smooth(method = "lm") +
  geom_hline(aes(yintercept = 0), color = "red") +
  geom_boxplot(aes(Generation, ABS_FITNESS, group = Generation), width = 1.5, alpha = .5) +
  stat_regline_equation(label.y = 2.1) +
  labs(x = "Generation",
       y = "Absolute Fitness",
       title = "Evolution of Centered Absolute Fitness 2022-2023")
ggsave("/bigdata/koeniglab/jmarz001/Ag-Competition/results/02b_Generational_Change_in_Single_Fitness_2022_2023.png")

### 02bi_Single_Fecundity_over_Generations.R

g <- ggplot(Single_2022_2023, aes(Generation, FECUNDITY, add = "reg.line")) +
  geom_jitter(alpha = .5) +
  geom_smooth(method = lm) +
  geom_hline(aes(yintercept = 0), color = "red") +
  geom_boxplot(aes(Generation, FECUNDITY, group = Generation), width = 1.5, alpha = .5) +
  stat_regline_equation(label.y = 4) +
  labs(x = "Generation",
       y = "Average Fecundity",
       title = "Evolution of Centered Fecundity 2022-2023")
ggsave("/bigdata/koeniglab/jmarz001/Ag-Competition/results/02bi_Generational_Change_in_Single_Fecundity_2022_2023.png")

### 02bii_Single_FT_over_Generations.R

h <- ggplot(Single_2022_2023, aes(Generation, FT, add = "reg.line")) +
  geom_jitter(alpha = .5) +
  geom_smooth(method = lm) +
  stat_regline_equation(label.x = 40) +
  geom_boxplot(aes(Generation, FT, group = Generation), width = 1.5, alpha = .5)+
  labs(x = "Generation",
       y = "Average Flowering Time (Days)",
       title = "Evolution of Average Flowering Time 2022-2023")
ggsave("/bigdata/koeniglab/jmarz001/Ag-Competition/results/02bii_Generational_Change_in_Single_FT_2022_2023.png")

### 02biii_Single_TW_over_Generations.R

i <- ggplot(Single_2022_2023, aes(Generation, TOTAL_MASS, add = "reg.line")) +
  geom_jitter(alpha = .5) +
  geom_smooth(method = lm) +
  geom_boxplot(aes(Generation,TOTAL_MASS, group = Generation), width = 1.5, alpha = .5)+
  stat_regline_equation(label.y = 155) +
  labs(x = "Generation",
       y = "Average Total Seed Weight (g)",
       title = "Evolution of Average Total Weight (grams) 2022-2023")
ggsave("/bigdata/koeniglab/jmarz001/Ag-Competition/results/02biii_Generational_Change_in_Single_TW_2022_2023.png")

### 02biiii_Single_100SW_over_Generations.R

j <- ggplot(Single_2022_2023, aes(Generation, SEED_WEIGHT_100, add = "reg.line")) +
  geom_jitter(alpha = .5) +
  geom_smooth(method = "lm") +
  geom_boxplot(aes(Generation, SEED_WEIGHT_100, group = Generation), width = 1.5, alpha = .5) +
  stat_regline_equation(label.y = 6.1) +
  labs(y = "Average 100 Seed Weight (grams)",
       title = "Evolution of Average 100 Seed Weight (grams) 2022-2023")
ggsave("/bigdata/koeniglab/jmarz001/Ag-Competition/results/02biiii_Generational_Change_in_Single_100SW_2022_2023.png")

### 02biiiii_Combined_Single_Evolution_Scatterplots.R

k <- arrangeGrob(a, b, c, d, e, f, g, h, i, j, top = "Evolution of Our Four Measured Phenotypes (Single)", nrow = 2, ncol = 5)
ggsave("/bigdata/koeniglab/jmarz001/Ag-Competition/results/02biiiii_Combined_Single_Evolution_Scatterplots.png",k, width = 30, height = 14)












### 02ci_Mixed_Relative_Fitness_to_Atlas.R

mf <- ggplot(Mixed_2021_2022, aes(Generation, Fitness, add = "reg.line")) +
  geom_jitter(alpha = .4) +
  ylim(0, 40000) +
  stat_regline_equation() +
  geom_hline(yintercept = 21347.22, color = "red") +
  geom_smooth(method = "lm") +
  labs(title = "Relative Fitness of Atlas Compared to Mixed Condition")
ggsave("scripts/plotting/02ci_Mixed_Relative_Fitness_to_Atlas_2021_2022.png")

jf <- ggplot(Averaged_Full_2021_2022, aes(Generation, Fitness, add = "reg.line", color = Condition)) +
  geom_jitter(alpha = .6) +
  stat_regline_equation() +
  geom_smooth(method = "lm" ) +
  geom_hline(yintercept = 21347.22, color = "red") +
  labs(title = "Relative Fitness of Atlas Compared to Both Conditions")
ggsave("scripts/plotting/02cii_Relative_Fitness_to_Atlas_Combined_2021_2022.png", width = 12, height = 7)

### 02cii_Relative_Fitness_to_Atlas_Combined.R

g <- arrangeGrob(sf, mf, jf, nrow = 2)
ggsave("scripts/plotting/02cii_Relative_Fitness_to_Atlas_Combined_2021_2022.png", g, width = 12, height = 7)
