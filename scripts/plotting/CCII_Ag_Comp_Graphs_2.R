#!/usr/bin/env Rscript
#SBATCH --ntasks=1
#SBATCH --mem=30G
#SBATCH --time=02:00:00
#SBATCH --output=/rhome/jmarz001/bigdata/Ag-Competition/scripts/CCII_Ag_Comp_Graphs_2.stdout
#SBATCH -p koeniglab

# This script plots scatterplots with linear regressions for each trait in the single subpopulation, as well as distributions for these traits

library(tidyverse)
library(ggpubr)
library(ggplot2)

setwd("/bigdata/koeniglab/jmarz001/Ag-Competition/data")
fitness_df <- read_delim("FITNESS.tsv")

# Making single dataframe

s_fit <- fitness_df[fitness_df$Condition == "single", ]

# Adding columns for standardized data

s_fit <- s_fit %>% mutate(stand_fit = scale(FECUNDITY),
                          stand_fec = scale(FITNESS),
                          stand_tw = scale(TOTAL_MASS),
                          stand_100 = scale(SEED_WEIGHT_100),
                          stand_ft = scale(FT))

### 02_standard_fit_over_gen.R

a <- ggplot(s_fit, aes(Generation, stand_fit, add = "reg.line")) +
geom_jitter(alpha = .5) +
  geom_smooth(method = lm) +
  geom_hline(aes(yintercept = 0), color = "red") +
  geom_boxplot(aes(Generation, stand_fit, group = Generation), width = 1.5, alpha = .5) +
  stat_regline_equation() +
  scale_y_continuous(breaks = seq(-4, 4, 2), limits = c(-4,4)) +
  labs(x = "Generation",
       y = "Average Fitness") +
  theme_bw()

### 02_standard_fec_over_gen.R

b <- ggplot(s_fit, aes(Generation, stand_fec, add = "reg.line")) +
  geom_jitter(alpha = .5) +
  geom_smooth(method = lm) +
  geom_hline(aes(yintercept = 0), color = "red") +
  geom_boxplot(aes(Generation, stand_fec, group = Generation), width = 1.5, alpha = .5) +
  stat_regline_equation() +
  labs(x = "Generation",
       y = "Average Fecundity") +
  scale_y_continuous(breaks = seq(-4, 4, 2), limits = c(-4,4)) +
  theme_bw()

### 02_standard_fit_over_gen.R

c <- ggplot(s_fit, aes(Generation, stand_ft, add = "reg.line")) +
  geom_jitter(alpha = .5) +
  geom_smooth(method = lm) +
  geom_hline(aes(yintercept = 0), color = "red") +
  stat_regline_equation() +
  geom_boxplot(aes(Generation, stand_ft, group = Generation), width = 1.5, alpha = .5) +
  labs(x = "Generation",
       y = "Average Flowering Time (Days)") +
  scale_y_continuous(breaks = seq(-4, 4,2), limits = c(-4,4)) +
  theme_bw()

### 02_standard_tw_over_gen.R

d <- ggplot(s_fit, aes(Generation, stand_tw, add = "reg.line")) +
  geom_jitter(alpha = .5) +
  geom_smooth(method = lm) +
  geom_hline(aes(yintercept = 0), color = "red") +
  geom_boxplot(aes(Generation, stand_tw, group = Generation), width = 1.5, alpha = .5)+
  stat_regline_equation() +
  labs(x = "Generation",
       y = "Average Total Seed Weight (g)") +
  scale_y_continuous(breaks = seq(-4, 4,2), limits = c(-4,4)) +
  theme_bw()

### 02_standard_100sw_over_gen.R

e <- ggplot(s_fit, aes(Generation, stand_100, add = "reg.line")) +
  geom_jitter(alpha = .5) +
  geom_smooth(method = lm) +
  geom_hline(aes(yintercept = 0), color = "red") +
  geom_boxplot(aes(Generation, stand_100, group = Generation), width = 1.5, alpha = .5)+
  stat_regline_equation() +
  scale_y_continuous(breaks = seq(-4, 4, 2), limits = c(-4,4)) +
  labs(x = "Generation",
       y = "Average 100 Seed Weight (g)") +
  theme_bw()

### 02_combined_single_evolution_scatterplots.R

y <- grid.arrange(a, b, c, d, e, top = "Evolution of Our Four Measured Phenotypes", nrow = 2, ncol = 3)
ggsave("/bigdata/koeniglab/jmarz001/Ag-Competition/results/combined_single_evolution_scatterplots.png",y, width = 14, height = 10)



### Unscaled scatterplots

s_fit <- fitness_df[fitness_df$Condition == "single", ]

### 02a_fit_over_gen.R

a <- ggplot(s_fit, aes(Generation, FITNESS, add = "reg.line")) +
  geom_jitter(alpha = .5) +
  geom_smooth(method = lm) +
  geom_hline(aes(yintercept = mean(FITNESS)), color = "red") +
  geom_boxplot(aes(Generation, FITNESS, group = Generation), width = 1.5, alpha = .5) +
  stat_regline_equation() +
  labs(x = "Generation",
       y = "Average Fitness") +
  theme_bw()

### 02a_fec_over_gen.R

b <- ggplot(s_fit, aes(Generation, FECUNDITY, add = "reg.line")) +
  geom_jitter(alpha = .5) +
  geom_smooth(method = lm) +
  geom_hline(aes(yintercept = mean(FECUNDITY)), color = "red") +
  geom_boxplot(aes(Generation, FECUNDITY, group = Generation), width = 1.5, alpha = .5) +
  stat_regline_equation() +
  labs(x = "Generation",
       y = "Average Fecundity") +
  theme_bw()

### 02a_ft_over_gen.R

c <- ggplot(s_fit, aes(Generation, FT, add = "reg.line")) +
  geom_jitter(alpha = .5) +
  geom_smooth(method = lm) +
  geom_hline(aes(yintercept = mean(FT)), color = "red") +
  stat_regline_equation() +
  geom_boxplot(aes(Generation, FT, group = Generation), width = 1.5, alpha = .5) +
  labs(x = "Generation",
       y = "Average Flowering Time (Days)") +
  theme_bw()

### 02a_tw_over_gen.R

d <- ggplot(s_fit, aes(Generation, TOTAL_MASS, add = "reg.line")) +
  geom_jitter(alpha = .5) +
  geom_smooth(method = lm) +
  geom_hline(aes(yintercept = mean(TOTAL_MASS)), color = "red") +
  geom_boxplot(aes(Generation, TOTAL_MASS, group = Generation), width = 1.5, alpha = .5)+
  stat_regline_equation() +
  labs(x = "Generation",
       y = "Average Total Seed Weight (g)") +
  theme_bw()

### 02a_100sw_over_gen.R

e <- ggplot(s_fit, aes(Generation, SEED_WEIGHT_100, add = "reg.line")) +
  geom_jitter(alpha = .5) +
  geom_smooth(method = lm) +
  geom_hline(aes(yintercept = mean(SEED_WEIGHT_100)), color = "red") +
  geom_boxplot(aes(Generation, SEED_WEIGHT_100, group = Generation), width = 1.5, alpha = .5)+
  stat_regline_equation() +
  labs(x = "Generation",
       y = "Average 100 Seed Weight (g)") +
  theme_bw()

### 02_combined_single_evolution_scatterplots.R

y <- grid.arrange(a, b, c, d, e, top = "Evolution of Our Four Measured Phenotypes", nrow = 2, ncol = 3)
ggsave("/bigdata/koeniglab/jmarz001/Ag-Competition/results/combined_single_evolution_scatterplots_unscaled.png",y, width = 14, height = 10)







### 02b_Fecundity_Distributions_Over_Generations.R

ggplot(Average_Haplo_rep, aes(x = Fecundity, color = Condition, fill = Condition)) +
  geom_histogram(position = "identity", bins = 45, alpha = .3) +
  labs(x = "Fecundity",
       y = "Frequency",
       title = "Fecundity Over Generations") +
  facet_wrap(~Generation, scales = "free_x")
ggsave("scripts/plotting/02b_Fecundity_Over_Generations_Distributions_2022_2023.png")

Average_Haplo_rep$Generation <- as.factor(Average_Haplo_rep$Generation)
ggplot(Average_Haplo_rep, aes(x = Fecundity, fill = Generation, group = Generation)) +
  geom_histogram(alpha = .5, position = "identity", binwidth = 70) +
  scale_fill_brewer(palette = "Blues")
ggsave("scripts/plotting/02b_Overlapping_Distributions_Fecundity_Over_Generations_2022_2023.png")

### 02bi_100_SW_Distributions.R

ggplot(Average_Haplo_rep, aes(x = `100 seed weight`, group = Generation, fill = Generation)) +
  geom_histogram(alpha =.5, position = "identity") +
  scale_fill_brewer(palette = "Blues")
ggsave("scripts/plotting/02bi_Overlapping_Distributions_for_Generational_Change_100SW_Distributions_2022_2023.png")

### 02c_Single_Relative_Fitness_to_Atlas.R

sf <- ggplot(Rep_Single, aes(Generation, Fitness, add = "reg.line")) +
  geom_jitter(alpha = .4) +
  stat_regline_equation(label.x = 40) +
  ylim(0, 40000) +
  geom_smooth(method = "lm") +
  geom_hline(yintercept = 21347.22, color = "red") +
  labs(title = "Relative Fitness of Atlas Compared to Single Condition") 
ggsave("scripts/plotting/02c_Single_Relative_Fitness_to_Atlas_2022_2023.png")

### 02ci_Mixed_Relative_Fitness_to_Atlas.R

mf <- ggplot(Rep_Mixed, aes(Generation, Fitness, add = "reg.line")) +
  geom_jitter(alpha = .4) +
  ylim(0, 40000) +
  stat_regline_equation() +
  geom_hline(yintercept = 21347.22, color = "red") +
  geom_smooth(method = "lm") +
  labs(title = "Relative Fitness of Atlas Compared to Mixed Condition")
ggsave("scripts/plotting/02ci_Mixed_Relative_Fitness_to_Atlas_2022_2023.png")

jf <- ggplot(Average_Haplo_rep, aes(Generation, Fitness, add = "reg.line", color = Condition)) +
  geom_jitter(alpha = .6) +
  stat_regline_equation() +
  geom_smooth(method = "lm" ) +
  geom_hline(yintercept = 21347.22, color = "red") +
  labs(title = "Relative Fitness of Atlas Compared to Both Conditions")
ggsave("scripts/plotting/02cii_Relative_Fitness_to_Atlas_Combined_2022_2023.png", width = 12, height = 7)

### 02cii_Relative_Fitness_to_Atlas_Combined.R

g <- arrangeGrob(sf, mf, jf, nrow = 2)
ggsave("scripts/plotting/02cii_Relative_Fitness_to_Atlas_Combined_2022_2023.png", g, width = 12, height = 7)

write_delim(Rep_Single, "Rep_Single")