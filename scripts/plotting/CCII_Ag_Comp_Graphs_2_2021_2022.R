library(tidyverse)
library(readr)
library(dplyr)
library(ggpubr)
library(ggplot2)
library(tidyr)
library(car)
library(gridExtra)

#fitness = (Plot Germination * Fecundity)
### 02_Single_Fitness_over_Generation.R

fa <- ggplot(Single_2021_2022, aes(Generation, Centered_Fit, add = "reg.line")) +
  geom_jitter(alpha = .5) +
  geom_smooth(method = lm) +
  geom_hline(aes(yintercept = 0), color = 'red') +
  geom_boxplot(aes(Generation, Centered_Fit, group = Generation), width = 1.5, alpha = .5) + 
  stat_regline_equation(label.y = 35000) +
  scale_y_continuous(breaks = seq(-10000, 50000, 5000)) +
  labs(x = "Generation",
       y = "Average Fitness") 
ggsave("scripts/plotting/02_Generational_Change_in_Single_Fitness_2021_2022.png")

### 02a_Single_Fecundity_over_Generations.R

fb <- ggplot(Single_2021_2022, aes(Generation, Centered_Fec, add = "reg.line")) +
  geom_jitter(alpha = .5) +
  geom_smooth(method = lm) +
  geom_hline(aes(yintercept = 0), color = 'red') +
  geom_boxplot(aes(Generation, Centered_Fec, group = Generation), width = 1.5, alpha = .5) +
  stat_regline_equation(label.y = 3000) +
  scale_y_continuous(breaks = seq(-2000, 5000, 1000)) +
  labs(x = "Generation",
       y = "Average Fecundity") 
ggsave("scripts/plotting/02a_Generational_Change_in_Single_Fecundity_2021_2022.png")

### 02ai_Single_FT_over_Generations.R

fc <- ggplot(Single_2021_2022, aes(Generation, Centered_FT, add = "reg.line")) +
  geom_jitter(alpha = .5) +
  geom_smooth(method = lm) +
  geom_hline(aes(yintercept = 0), color = 'red') +
  geom_boxplot(aes(Generation, Centered_FT, group = Generation), width = 1.5, alpha = .5)+
  stat_regline_equation(label.y = 30) +
  scale_y_continuous(breaks = seq(-40, 30, 10)) +
  labs(x = "Generation",
       y = "Average Flowering Time (Days)")
ggsave("scripts/plotting/02ai_Generational_Change_in_Single_FT_2021_2022.png")

### 02aii_Single_TW_over_Generations.R

fd <- ggplot(Single_2021_2022, aes(Generation, Centered_TW, add = "reg.line")) +
  geom_jitter(alpha = .5) +
  geom_smooth(method = lm) +
  geom_hline(aes(yintercept = 0), color = 'red') +
  geom_boxplot(aes(Generation,Centered_TW, group = Generation), width = 1.5, alpha = .5)+
  stat_regline_equation(label.y = 90) +
  labs(x = "Generation",
       y = "Average Total Seed Weight (g)") 
ggsave("scripts/plotting/02aii_Generational_Change_in_Single_TW_2021_2022.png")

### 02aiii_Combined_Single_Evolution_Scatterplots.R

g <- grid.arrange(fd, fc, fb, fa, top = "Evolution of Our Four Measured Phenotypes")
ggsave("scripts/plotting/02aiii_Combined_Single_Evolution_Scatterplots_2021_2022.png",g, width = 10, height = 8)

### 02b_Fecundity_Distributions_Over_Generations.R

ggplot(Averaged_Full_2021_2022, aes(x = Fecundity, color = Condition, fill = Condition)) +
  geom_histogram(position = "identity", bins = 45, alpha = .3) +
  labs(x = "Fecundity",
       y = "Frequency",
       title = "Fecundity Over Generations") +
  facet_wrap(~Generation, scales = "free_x")
ggsave("scripts/plotting/02b_Fecundity_Over_Generations_Distributions_2021_2022.png")

Averaged_Full_2021_2022$Generation <- as.factor(Averaged_Full_2021_2022$Generation)
ggplot(Averaged_Full_2021_2022, aes(x = Fecundity, fill = Generation, group = Generation)) +
  geom_histogram(alpha = .5, position = 'identity', binwidth = 70) +
  scale_fill_brewer(palette = "Blues")
ggsave("scripts/plotting/02b_Overlapping_Distributions_Fecundity_Over_Generations_2021_2022.png")

### 02bi_100_SW_Distributions.R

ggplot(Averaged_Full_2021_2022, aes(x = "100_seed_weight", group = Generation, fill = Generation)) +
  geom_histogram(alpha =.5, position = 'identity') +
  scale_fill_brewer(palette = "Blues")
ggsave("scripts/plotting/02bi_Overlapping_Distributions_for_Generational_Change_100SW_Distributions_2021_2022.png")

### 02c_Single_Relative_Fitness_to_Atlas.R

sf <- ggplot(Single_2021_2022, aes(Generation, Fitness, add = "reg.line")) +
  geom_jitter(alpha = .4) +
  stat_regline_equation(label.x = 40) +
  ylim(0, 40000) +
  geom_smooth(method = 'lm') +
  geom_hline(yintercept = 21347.22, color = "red") +
  labs(title = "Relative Fitness of Atlas Compared to Single Condition") 
ggsave("scripts/plotting/02c_Single_Relative_Fitness_to_Atlas_2021_2022.png")

### 02ci_Mixed_Relative_Fitness_to_Atlas.R

mf <- ggplot(Mixed_2021_2022, aes(Generation, Fitness, add = "reg.line")) +
  geom_jitter(alpha = .4) +
  ylim(0, 40000) +
  stat_regline_equation() +
  geom_hline(yintercept = 21347.22, color = "red") +
  geom_smooth(method = 'lm') +
  labs(title = "Relative Fitness of Atlas Compared to Mixed Condition")
ggsave("scripts/plotting/02ci_Mixed_Relative_Fitness_to_Atlas_2021_2022.png")

jf <- ggplot(Averaged_Full_2021_2022, aes(Generation, Fitness, add = "reg.line", color = Condition)) +
  geom_jitter(alpha = .6) +
  stat_regline_equation() +
  geom_smooth(method = 'lm' ) +
  geom_hline(yintercept = 21347.22, color = "red") +
  labs(title = "Relative Fitness of Atlas Compared to Both Conditions")
ggsave("scripts/plotting/02cii_Relative_Fitness_to_Atlas_Combined_2021_2022.png", width = 12, height = 7)

### 02cii_Relative_Fitness_to_Atlas_Combined.R

g <- arrangeGrob(sf, mf, jf, nrow = 2)
ggsave("scripts/plotting/02cii_Relative_Fitness_to_Atlas_Combined_2021_2022.png", g, width = 12, height = 7)
