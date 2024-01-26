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

#### 02_Single_Fitness_Over_Generations.R
#fitness = (Plot Germination * Fecundity)

Rep_Single$Generation <- as.numeric(Rep_Single$Generation)
ggplot(Rep_Single, aes(Generation, Fitness), add = "reg.line") +
  geom_jitter(alpha =.3) + 
  geom_smooth(method = lm) +
  stat_regline_equation(label.y = 35000) +
  labs(x = "Generation",
       y = "Plot Fitness (Plot Germination * Fecundity)",
       title = "Evolution of Average Plot Fitness") 

ggsave("02a_Single_Fitness_Over_Generations.png")

ggplot(Rep_Single, aes(Generation, Fecundity)) +
  geom_jitter(alpha =.3) + 
  geom_smooth(method = lm) +
  labs(x = "Generation",
       y = "Fecundity",
       title = "Average Fecundity over Generations")

### 2ai_Avg_FT_over_Generations.R

ggplot(Average_Haplo_rep, aes(Generation, FT_DAYS, add = "reg.line")) +
  geom_jitter(alpha = .3) +
  stat_regline_equation(label.y = 140) +
  geom_smooth(method = lm) +
  labs(y = "Flowering Time (Days)",
       title = " Average Flowering Time over Generations") +
  facet_wrap(~Condition)

ggsave("2ai_FT_over_Generations.png")

### 2aii_Avg_Total_Seed_Weight_Over_Time

fe <- ggplot(Average_Haplo_rep, aes(Generation, `Brown Bag Weight`, add = "reg.line")) +
  geom_jitter(alpha = .3) +
  geom_smooth(method = lm) +
  stat_regline_equation(label.y = 200) +
  labs(x = "Generation",
       y = "Average Total Seed Weight (grams)") +
  facet_wrap(~Condition)

ggsave("2aii_Overall_Trend_Avg_Weight_Over_Generations.png")

ff <- ggplot(Average_Haplo_rep, aes(Generation, `Brown Bag Weight`, col = Condition)) +
  geom_jitter(alpha = .5) +
  geom_smooth(method = lm) +
  labs(x = "Generation",
       y = "Average Total Seed Weight (grams)") 

ggsave("2aii_Total_Avg_Weight_Over_Generations_Between_Conditions.png")

Compare_TW <- grid.arrange(fe,ff)
ggsave("Average_Total_Weight_Comparisons.png", Compare_TW)

### 02aiii_Avg_FT_Over_Generations.R

fe <- ggplot(Average_Haplo_rep, aes(Generation, FT_DAYS, add = "reg.line")) +
  geom_jitter(alpha = .3) +
  geom_smooth(method = lm) +
  stat_regline_equation(label.y = 130) +
  labs(x = "Generation",
       y = "Average Flowering Time (Days)") +
  facet_wrap(~Condition)

ff <- ggplot(Average_Haplo_rep, aes(Generation, FT_DAYS, col = Condition)) +
  geom_jitter(alpha = .5) +
  geom_smooth(method = lm) +
  labs(x = "Generation",
       y = "Average Total Seed Weight (grams)") 


Compare_TW <- grid.arrange(fe,ff)
ggsave("Average_FT_Comparisons.png", Compare_TW)
       
### 02aiiii_Avg_Fec_Over_Generations.R

fe <- ggplot(Average_Haplo_rep, aes(Generation, Fecundity, add = "reg.line")) +
  geom_jitter(alpha = .3) +
  geom_smooth(method = lm) +
  stat_regline_equation(label.y = 3500) +
  labs(x = "Generation",
       y = "Average Fecundity") +
  facet_wrap(~Condition)

ff <- ggplot(Average_Haplo_rep, aes(Generation, Fecundity, col = Condition)) +
  geom_jitter(alpha = .5) +
  geom_smooth(method = lm) +
  labs(x = "Generation",
       y = "Average Fecundity") 

Compare_TW <- grid.arrange(fe,ff)
ggsave("Average_Fecundity_Comparisons.png", Compare_TW)

### 02aiiiii_Avg_Fit_Over_Generations.R

fe <- ggplot(Average_Haplo_rep, aes(Generation, Fitness, add = "reg.line")) +
  geom_jitter(alpha = .3) +
  geom_smooth(method = lm) +
  stat_regline_equation(label.y = 31000) +
  labs(x = "Generation",
       y = "Average Fitness") +
  facet_wrap(~Condition)

ff <- ggplot(Average_Haplo_rep, aes(Generation, Fitness, col = Condition)) +
  geom_jitter(alpha = .5) +
  geom_smooth(method = lm) +
  labs(x = "Generation",
       y = "Average Fitness") 

Compare_TW <- grid.arrange(fe,ff)
ggsave("Average_Fitness_Comparisons.png", Compare_TW)

### 02b_Single_Fitness_over_Generation.R

fa <- ggplot(Rep_Single, aes(Generation, Fitness, add = "reg.line")) +
geom_jitter(alpha = .5) +
  geom_smooth(method = lm) +
  stat_regline_equation(label.y = 35000) +
  labs(x = "Generation",
       y = "Average Fitness",
       title = "Average Fitness Over Generations") 

ggsave("02b_Single_Fitness_over_Generation.png")

### 02b_Single_Fecundity_over_Generations

fb <- ggplot(Rep_Single, aes(Generation, Fecundity, add = "reg.line")) +
  geom_jitter(alpha = .5) +
  geom_smooth(method = lm) +
  stat_regline_equation(label.y = 4000) +
  labs(x = "Generation",
       y = "Average Fecundity",
       title = "Average Fecundity Over Generations") 

### 02b_Single_FT_over_Generations

fc <- ggplot(Rep_Single, aes(Generation, FT_DAYS, add = "reg.line")) +
  geom_jitter(alpha = .5) +
  geom_smooth(method = lm) +
  stat_regline_equation(label.y = 145) +
  labs(x = "Generation",
       y = "Average Flowering Time (Days)",
       title = "Average Flowering Time Over Generations")

### 02b_Single_Total_Weight_over_Generations

fd <- ggplot(Rep_Single, aes(Generation, `Brown Bag Weight`, add = "reg.line")) +
  geom_jitter(alpha = .5) +
  geom_smooth(method = lm) +
  stat_regline_equation(label.y = 160) +
  labs(x = "Generation",
       y = "Average Total Weight (g)",
       title = "Average Total Weight Over Generations")

grid.arrange(fa, fb, fc, fd, ncol = 2, nrow = 2)

### 02c_Fecundity_Distributions_Over_Generations.R

ggplot(Average_Haplo_rep, aes(x = Fecundity)) +
  geom_histogram(bins = 45) +
  labs(x = "Fecundity",
       y = "Frequency",
       title = "Fecundity Over Generations") +
  facet_wrap(~Generation)

ggsave("02c_Fecundity_Over_Generations_Distributions.png")

Average_Haplo_rep$Generation <- as.factor(Average_Haplo_rep$Generation)
ggplot(Average_Haplo_rep, aes(x = Fecundity, fill = Generation, group = Generation)) +
  geom_histogram(alpha = .5, position = 'identity', binwidth = 70) +
  scale_fill_brewer(palette = "Blues")

ggsave("02c_Overlapping_Histograms_Fecundity_Over_Generations.png")

ANOVA_Fedundity <- aov(Fecundity ~ as.factor(Generation), Full_Data)
summary(ANOVA_Fedundity)
TukeyHSD(ANOVA_Fedundity)
levene.test(Fecundity ~ as.factor(Generation, Full_Data))

### 02ci_100_SW_Distributions

ggplot(Side_By_Side_Replicates, aes(x = Avg_100_SW)) +
  geom_histogram() +
  labs(x = "100 Seed Weight",
       y = "Frequency",
       title = "Genertaional Changes in Average 100 Seed Weight") +
  facet_wrap(~Generation)

ggsave("02ci_100SW_Over_Generations_Distributions.png")

ggplot(Average_Haplo_rep, aes(x = `100 seed weight`, group = Generation, fill = Generation)) +
  geom_histogram(alpha =.5, position = 'identity') +
  scale_fill_brewer(palette = "Blues")

ANOVA_100SW_Over_Generations <- aov(`100 seed weight`~as.factor(Generation), Average_Haplo_rep)
summary(ANOVA_100SW_Over_Generations)
TukeyHSD(ANOVA_100SW_Over_Generations)

### 02di_Relative_Single_Fitness_to_Atlas.R
## overall fitness

ggplot(Rep_Single, aes(Generation, Fitness, add = "reg.line")) +
  geom_jitter() +
  stat_regline_equation(label.x = 40) +
  ylim(0, 40000) +
  geom_smooth(method = 'lm') +
  geom_hline(yintercept = 21347.22, color = "red") +
  labs(y = "Relative Fitness",
       title = "Relative Fitness Over Generations") 

ggsave("02di_Single_Relative_Fitness_to_Atlas_Faceted_By_Generation.png")

## mixed fitness

ggplot(Rep_Mixed, aes(Generation, Fitness)) +
  geom_jitter() +
  ylim(0, 30000) +
  geom_hline(yintercept = 21347.22, color = "red") +
  geom_smooth(method = lm) +
  labs(y = "Relative Fitness",
       title = "Relative Fitness Over Generations (Mixed)") 