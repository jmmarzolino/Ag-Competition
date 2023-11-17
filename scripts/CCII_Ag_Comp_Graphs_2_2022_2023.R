library(tidyverse)
library(readr)
library(dplyr)
library(ggpubr)
library(ggplot2)
library(tidyr)

### Load Data

Full_Data <- read_delim("~/Documents/GitHub/Ag-Competition/Full_Data")
Side_By_Side_Replicates <- read_delim("~/Documents/GitHub/Ag-Competition/Side_By_Side_Replicates")
Rep_Mixed <- read_delim("~/Documents/GitHub/Ag-Competition/Rep_Mixed")
Rep_Single <- read_delim("~/Documents/GitHub/Ag-Competition/Rep_Single")

#### 02_Single_Fitness_Over_Generations.R
#fitness = (Plot Germination * Fecundity)

Rep_Single$Generation <- as.numeric(Rep_Single$Generation)
ggplot(Rep_Single, aes(Generation, Avg_Fit)) +
  geom_jitter(alpha =.3) + 
  geom_smooth(method = lm) +
  labs(x = "Generation",
       y = "Fitness (Plot Germination * Fecundity)",
       title = "Fitness over Generations")

ggsave("02a_Single_Fitness_Over_Generations.png")

### 2ai_FT_over_Generations.R

ggplot(Side_By_Side_Replicates, aes(Generation, FT_DAYS)) +
  geom_jitter(alpha = .3) +
  geom_smooth(method = lm) +
  labs(y = "Flowering Time (Days)",
       title = "Flowering Time over Generations") +
  facet_wrap(~Condition)

ggsave("2ai_FT_over_Generations.png")

### 2aii_Total_Seed_Weight_Over_Time

ggplot(Side_By_Side_Replicates, aes(Generation, Avg_Total_Weight)) +
  geom_jitter(alpha = .3) +
  geom_smooth(method = "lm") +
  labs(x = "Generation",
       y = "Average Total Seed Weight (grams)")

ggsave("2aii_Overall_Trend_Avg_Weight_Over_Generations.png")

ggplot(Side_By_Side_Replicates, aes(Generation, Avg_Total_Weight, add = "reg.line")) +
  geom_jitter(alpha = .3) +
  geom_smooth(method = "lm") +
  stat_regline_equation(label.y =200) +
  labs(x = "Generation",
       y = "Average Total Seed Weight (grams)") +
  facet_wrap(~Condition)

ggsave("2aii_Total_Avg_Weight_Over_Generations_Between_Conditions.png")

### 02b_Single_Fitness_over_Generation.R

ggplot(Rep_Single, aes(Generation, Avg_Fit, add = "reg.line")) +
geom_jitter(alpha = .5) +
  geom_smooth(method = lm) +
  stat_regline_equation(label.y = 35000) +
  labs(x = "Generation",
       y = "Average Fitness",
       title = "Average Fitness Over Generation") 

ggsave("02b_Single_Fitness_over_Generation.png")

### 02c_Fecundity_Distributions_Over_Generations.R

ggplot(Side_By_Side_Replicates, aes(x = Avg_Fec)) +
  geom_histogram() +
  labs(x = "Fecundity",
       y = "Frequency",
       title = "Fecundity Over Generations") +
  facet_wrap(~Generation)

ggsave("02c_Fecundity_Over_Generations_Distributions.png")

Side_By_Side_Replicates$Generation <- as.factor(Side_By_Side_Replicates$Generation)
ggplot(Side_By_Side_Replicates, aes(x = Avg_Fec, fill = Generation, group = Generation)) +
  geom_histogram(alpha = .5, position = 'identity', binwidth = 70) +
  scale_fill_brewer(palette = "Blues")

ggsave("02c_Overlapping_Histograms_Fecundity_Over_Generations.png")

ANOVA_Fedundity <- aov(Fecundity ~ as.factor(Generation), Full_Data)
summary(ANOVA_Fedundity)
TukeyHSD(ANOVA_Fedundity)

### 02ci_100_SW_Distributions

ggplot(Side_By_Side_Replicates, aes(x = `100 seed weight`)) +
  geom_histogram() +
  labs(x = "100 Seed Weight",
       y = "Frequency",
       title = "100 Seed Weight Over Time") +
  facet_wrap(~Generation)

ggsave("02ci_100SW_Over_Generations_Distributions.png")

ggplot(Side_By_Side_Replicates, aes(x = `100 seed weight`, group = Generation, fill = Generation)) +
  geom_histogram(alpha =.5, position = 'identity') +
  scale_fill_brewer(palette = "Blues")

ANOVA_100SW_Over_Generations <- aov(`100 seed weight`~as.factor(Generation), Side_By_Side_Replicates)
summary(ANOVA_100SW_Over_Generations)
TukeyHSD(ANOVA_100SW_Over_Generations)

### 02di_Relative_Fitness_to_Atlas.R

ggplot(Rep_Single, aes(Generation, Avg_Fit)) +
  geom_jitter() +
  ylim(0, 30000) +
  geom_hline(yintercept = 17834.67, color = "red") +
  labs(y = "Relative Fitness",
       title = "Relative Fitness Over Generations") +
  facet_wrap(~Generation)

ggsave("02di_Relative_Fitness_to_Atlas_Faceted_By_Generation.png")

ggplot(Rep_Single, aes(Generation, Avg_Fit)) +
  geom_jitter() +
  ylim(0, 30000) +
  geom_hline(yintercept = 17834.67, color = "red") +
  geom_smooth(method = lm) +
  labs(y = "Relative Fitness",
       title = "Relative Fitness Over Generations")

ggsave("Relative_Fitness_to_Atlas.png")
