#!/usr/bin/env Rscript


#SBATCH --mem=30G
#SBATCH --time=02:00:00
#SBATCH --output=/rhome/jmarz001/bigdata/Ag-Competition/scripts/competition1.stdout
#SBATCH -p koeniglab

library(tidyverse)


library(ggpubr)



library(car)
library(gridExtra)

# Mixed Scatterplots 2021-2022

Mixed_2021_2022 <- PHENO_FULL_AVERAGE %>% filter(Condition == "mixed" & Exp_year == 2022)

#fitness = (Plot Germination * Fecundity)
### 02a_Mixed_Fitness_over_Generation.R

Mixed_2021_2022$Generation <- as.numeric(Mixed_2021_2022$Generation)
a <- ggplot(Mixed_2021_2022, aes(x= Generation, y = ABS_FITNESS)) +
  geom_jitter(alpha = .5) +
  geom_smooth(method = "lm") +
  geom_hline(aes(yintercept = 0), color = "red") +
  geom_boxplot(aes(Generation, ABS_FITNESS, group = Generation), width = 1.5, alpha = .5) +
  stat_regline_equation(label.y = 2) +
  labs(x = "Generation",
       y = "Absolute Fitness",
       title = "Evolution of Centered Absolute Fitness 2021-2022")
ggsave("scripts/plotting/02a_Generational_Change_in_Mixed_Fitness_2021_2022.png")

### 02ai_Mixed_Fecundity_over_Generations.R

b <- ggplot(Mixed_2021_2022, aes(Generation, FECUNDITY, add = "reg.line")) +
  geom_jitter(alpha = .5) +
  geom_smooth(method = lm) +
  geom_hline(aes(yintercept = 0), color = "red") +
  geom_boxplot(aes(Generation, FECUNDITY, group = Generation), width = 1.5, alpha = .5) +
  stat_regline_equation(label.y = 3) +
  labs(x = "Generation",
       y = "Average Fecundity",
       title = "Evolution of Centered Fecundity 2021-2022")
ggsave("scripts/plotting/02ai_Generational_Change_in_Mixed_Fecundity_2021_2022.png")

### 02aii_Mixed_FT_over_Generations.R

c <- ggplot(Mixed_2021_2022, aes(Generation, FT, add = "reg.line")) +
  geom_jitter(alpha = .5) +
  geom_smooth(method = lm) +
  stat_regline_equation(label.x = 40) +
  geom_boxplot(aes(Generation, FT, group = Generation), width = 1.5, alpha = .5)+
  labs(x = "Generation",
       y = "Average Flowering Time (Days)",
       title = "Evolution of Average Flowering Time 2021-2022")
ggsave("scripts/plotting/02aii_Generational_Change_in_Mixed_FT_2021_2022.png")

### 02aiii_Mixed_TW_over_Generations.R

d <- ggplot(Mixed_2021_2022, aes(Generation, TOTAL_MASS, add = "reg.line")) +
  geom_jitter(alpha = .5) +
  geom_smooth(method = lm) +
  geom_boxplot(aes(Generation,TOTAL_MASS, group = Generation), width = 1.5, alpha = .5)+
  stat_regline_equation(label.y = 135) +
  labs(x = "Generation",
       y = "Average Total Seed Weight (g)",
       title = "Evolution of Average Total Weight (grams) 2021-2022")
ggsave("scripts/plotting/02aiii_Generational_Change_in_Mixed_TW_2021_2022.png")

### 02aiiii_Mixed_100SW_over_Generations.R

e <- ggplot(Single_2021_2022, aes(Generation, SEED_WEIGHT_100, add = "reg.line")) +
  geom_jitter(alpha = .5) +
  geom_smooth(method = "lm") +
  geom_boxplot(aes(Generation, SEED_WEIGHT_100, group = Generation), width = 1.5, alpha = .5) +
  stat_regline_equation(label.x = 5) +
  labs(y = "Average 100 Seed Weight (grams)",
       title = "Evolution of Average 100 Seed Weight (grams) 2021-2022")
ggsave("scripts/plotting/02aiiii_Generational_Change_in_Mixed_100SW_2021_2022.png")

###################################################

# Mixed Scatterplots 2022-2023

Mixed_2022_2023 <- PHENO_FULL_AVERAGE %>% filter(Condition == "mixed" & Exp_year == 2023)

#fitness = (Plot Germination * Fecundity)
### 02b_Mixed_Fitness_over_Generation.R

Mixed_2022_2023$Generation <- as.numeric(Mixed_2022_2023$Generation)
f <- ggplot(Mixed_2022_2023, aes(x= Generation, y = ABS_FITNESS)) +
  geom_jitter(alpha = .5) +
  geom_smooth(method = "lm") +
  geom_hline(aes(yintercept = 0), color = "red") +
  geom_boxplot(aes(Generation, ABS_FITNESS, group = Generation), width = 1.5, alpha = .5) +
  stat_regline_equation(label.y = 2.1) +
  labs(x = "Generation",
       y = "Absolute Fitness",
       title = "Evolution of Centered Absolute Fitness 2022-2023")
ggsave("scripts/plotting/02b_Generational_Change_in_Mixed_Fitness_2022_2023.png")

### 02bi_Mixed_Fecundity_over_Generations.R

g <- ggplot(Mixed_2022_2023, aes(Generation, FECUNDITY, add = "reg.line")) +
  geom_jitter(alpha = .5) +
  geom_smooth(method = lm) +
  geom_hline(aes(yintercept = 0), color = "red") +
  geom_boxplot(aes(Generation, FECUNDITY, group = Generation), width = 1.5, alpha = .5) +
  stat_regline_equation(label.y = 4) +
  labs(x = "Generation",
       y = "Average Fecundity",
       title = "Evolution of Centered Fecundity 2022-2023")
ggsave("scripts/plotting/02bi_Generational_Change_in_Mixed_Fecundity_2022_2023.png")

### 02bii_Mixed_FT_over_Generations.R

h <- ggplot(Mixed_2022_2023, aes(Generation, FT, add = "reg.line")) +
  geom_jitter(alpha = .5) +
  geom_smooth(method = lm) +
  stat_regline_equation(label.x = 40) +
  geom_boxplot(aes(Generation, FT, group = Generation), width = 1.5, alpha = .5)+
  labs(x = "Generation",
       y = "Average Flowering Time (Days)",
       title = "Evolution of Average Flowering Time 2022-2023")
ggsave("scripts/plotting/02bii_Generational_Change_in_Mixed_FT_2022_2023.png")

### 02biii_Mixed_TW_over_Generations.R

i <- ggplot(Mixed_2022_2023, aes(Generation, TOTAL_MASS, add = "reg.line")) +
  geom_jitter(alpha = .5) +
  geom_smooth(method = lm) +
  geom_boxplot(aes(Generation,TOTAL_MASS, group = Generation), width = 1.5, alpha = .5)+
  stat_regline_equation(label.y = 155) +
  labs(x = "Generation",
       y = "Average Total Seed Weight (g)",
       title = "Evolution of Average Total Weight (grams) 2022-2023")
ggsave("scripts/plotting/02biii_Generational_Change_in_Mixed_TW_2022_2023.png")

### 02biiii_Mixed_100SW_over_Generations.R

j <- ggplot(Mixed_2022_2023, aes(Generation, SEED_WEIGHT_100, add = "reg.line")) +
  geom_jitter(alpha = .5) +
  geom_smooth(method = "lm") +
  geom_boxplot(aes(Generation, SEED_WEIGHT_100, group = Generation), width = 1.5, alpha = .5) +
  stat_regline_equation(label.y = 6.1) +
  labs(y = "Average 100 Seed Weight (grams)",
       title = "Evolution of Average 100 Seed Weight (grams) 2022-2023")
ggsave("scripts/plotting/02biiii_Generational_Change_in_Mixed_100SW_2022_2023.png")

### 02biiiii_Combined_Mixed_Evolution_Scatterplots.R

k <- arrangeGrob(a, b, c, d, e, f, g, h, i, j, nrow = 2, ncol = 5,top = "Evolution of Our Measured Phenotypes (Mixed)")
ggsave("scripts/plotting/02biiiii_Combined_Mixed_Evolution_Scatterplots.png",k, width = 30, height = 14)












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
