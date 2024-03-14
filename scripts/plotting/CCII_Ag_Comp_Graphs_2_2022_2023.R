library(tidyverse)
library(readr)
library(dplyr)
library(ggpubr)
library(ggplot2)
library(tidyr)
library(car)
library(gridExtra)

### Load Data

Full_Data <- read_delim("~/Documents/GitHub/Ag-Competition/Full_Data")
Average_Haplo_rep <- read_delim("~/Documents/GitHub/Ag-Competition/Average_Haplo_rep")
Rep_Mixed <- read_delim("~/Documents/GitHub/Ag-Competition/Rep_Mixed")
Rep_Single <- read_delim("~/Documents/GitHub/Ag-Competition/Rep_Single")

#fitness = (Plot Germination * Fecundity
### 02_Single_Fitness_over_Generation.R

fa <- ggplot(Rep_Single, aes(Generation, Fitness, add = "reg.line")) +
geom_jitter(alpha = .5) +
  geom_smooth(method = lm) +
  geom_boxplot(aes(Generation, Fitness, group = Generation), width = 1.5, alpha = .5) +
  stat_regline_equation(label.y = 35000) +
  scale_y_continuous(breaks = seq(0, 40000, 5000)) +
  labs(x = "Generation",
       y = "Average Fitness") 
ggsave("scripts/plotting/02_Generational_Change_in_Single_Fitness.png")

### 02a_Single_Fecundity_over_Generations.R

fb <- ggplot(Rep_Single, aes(Generation, Fecundity, add = "reg.line")) +
  geom_jitter(alpha = .5) +
  geom_smooth(method = lm) +
  geom_boxplot(aes(Generation, Fecundity, group = Generation), width = 1.5, alpha = .5) +
  stat_regline_equation(label.y = 4000) +
  labs(x = "Generation",
       y = "Average Fecundity") 
ggsave("scripts/plotting/02a_Generational_Change_in_Single_Fecundity.png")

### 02ai_Single_FT_over_Generations.R

fc <- ggplot(Rep_Single, aes(Generation, FT_DAYS, add = "reg.line")) +
  geom_jitter(alpha = .5) +
  geom_smooth(method = lm) +
  geom_boxplot(aes(Generation, FT_DAYS, group = Generation), width = 1.5, alpha = .5)+
  stat_regline_equation(label.y = 140) +
  labs(x = "Generation",
       y = "Average Flowering Time (Days)")
ggsave("scripts/plotting/02ai_Generational_Change_in_Single_FT.png")

### 02aii_Single_Yield_over_Generations.R

fd <- ggplot(Rep_Single, aes(Generation, `Brown Bag Weight`, add = "reg.line")) +
  geom_jitter(alpha = .5) +
  geom_smooth(method = lm) +
  geom_boxplot(aes(Generation,`Brown Bag Weight`, group = Generation), width = 1.5, alpha = .5)+
  stat_regline_equation(label.y = 160) +
  labs(x = "Generation",
       y = "Average Yield (g)")
ggsave("scripts/plotting/02aii_Generational_Change_in_Single_Yield.png")

### 02aiii_Combined_Single_Evolution_Scatterplots.R

g <- grid.arrange(fd, fc, fb, fa, top = "Evolution of Our Four Measured Phenotypes")
ggsave("scripts/plotting/02aiii_Combined_Single_Evolution_Scatterplots.png",g, width = 10, height = 8)

### 02b_Fecundity_Distributions_Over_Generations.R

ggplot(Average_Haplo_rep, aes(x = Fecundity, color = Condition, fill = Condition)) +
  geom_histogram(position = "identity", bins = 45, alpha = .3) +
  labs(x = "Fecundity",
       y = "Frequency",
       title = "Fecundity Over Generations") +
  facet_wrap(~Generation, scales = "free_x")
ggsave("scripts/plotting/02b_Fecundity_Over_Generations_Distributions.png")

Average_Haplo_rep$Generation <- as.factor(Average_Haplo_rep$Generation)
ggplot(Average_Haplo_rep, aes(x = Fecundity, fill = Generation, group = Generation)) +
  geom_histogram(alpha = .5, position = 'identity', binwidth = 70) +
  scale_fill_brewer(palette = "Blues")
ggsave("scripts/plotting/02b_Overlapping_Distributions_Fecundity_Over_Generations.png")

### 02bi_100_SW_Distributions.R

ggplot(Side_By_Side_Replicates, aes(x = Avg_100_SW)) +
  geom_histogram() +
  labs(x = "100 Seed Weight",
       y = "Frequency",
       title = "Genertaional Changes in Average 100 Seed Weight") +
  facet_wrap(~Generation)
ggsave("scripts/plotting/02bi_Distributions_for_Generational_Change_100SW_Distributions.png")

ggplot(Average_Haplo_rep, aes(x = `100 seed weight`, group = Generation, fill = Generation)) +
  geom_histogram(alpha =.5, position = 'identity') +
  scale_fill_brewer(palette = "Blues")
ggsave("scripts/plotting/02bi_Overlapping_Distributions_for_Generational_Change_100SW_Distributions.png")

### 02c_Single_Relative_Fitness_to_Atlas.R

sf <- ggplot(Rep_Single, aes(Generation, Fitness, add = "reg.line")) +
  geom_jitter(alpha = .4) +
  stat_regline_equation(label.x = 40) +
  ylim(0, 40000) +
  geom_smooth(method = 'lm') +
  geom_hline(yintercept = 21347.22, color = "red") +
  labs(title = "Relative Fitness of Atlas Compared to Single Condition") 
ggsave("scripts/plotting/02c_Single_Relative_Fitness_to_Atlas.png")

### 02ci_Mixed_Relative_Fitness_to_Atlas.R

mf <- ggplot(Rep_Mixed, aes(Generation, Fitness, add = "reg.line")) +
  geom_jitter(alpha = .4) +
  ylim(0, 40000) +
  stat_regline_equation() +
  geom_hline(yintercept = 21347.22, color = "red") +
  geom_smooth(method = 'lm') +
  labs(title = "Relative Fitness of Atlas Compared to Mixed Condition")
ggsave("scripts/plotting/02ci_Mixed_Relative_Fitness_to_Atlas.png")

jf <- ggplot(Average_Haplo_rep, aes(Generation, Fitness, add = "reg.line", color = Condition)) +
  geom_jitter(alpha = .6) +
  stat_regline_equation() +
  geom_smooth(method = 'lm' ) +
  geom_hline(yintercept = 21347.22, color = "red") +
  labs(title = "Relative Fitness of Atlas Compared to Both Conditions")
ggsave("scripts/plotting/02cii_Relative_Fitness_to_Atlas_Combined.png", width = 12, height = 7)

### 02cii_Relative_Fitness_to_Atlas_Combined.R

g <- arrangeGrob(sf, mf, jf, nrow = 2)
ggsave("scripts/plotting/02cii_Relative_Fitness_to_Atlas_Combined.png", g, width = 12, height = 7)
