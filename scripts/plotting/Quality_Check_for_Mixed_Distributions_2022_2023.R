#!/usr/bin/env Rscript
#SBATCH --mem=30G
#SBATCH --time=02:00:00
#SBATCH --output=/rhome/jmarz001/bigdata/Ag-Competition/scripts/Quality_Check_for_Mixed_Distributions_2022_2023.stdout
#SBATCH -p koeniglab

library(tidyverse)
library(ggpubr)
library(car)
library(gridExtra)
library(dunn.test)

### Load Data

Rep_Mixed <- read_delim("~/Documents/GitHub/Ag-Competition/Rep_Mixed")

###### Flowering Time Distributions

### Histograms for Avg_FT Over Generations (Mixed)

ggplot(Rep_Mixed, aes(x = FT)) +
  geom_histogram(binwidth = .6) +
  facet_wrap(~Generation) +
  labs(x = "Average Flowering Time (Days after sowing)",
       y = "Frequency",
       title = "Generational Change in Average Flowering Time (Mixed)")

### ANOVA and Homogeneity of Variance (Mixed)

leveneTest(FT ~ as.factor(Generation), Rep_Mixed)
ANOVA_FT <- aov(FT ~ as.factor(Generation), Rep_Mixed)
summary(ANOVA_FT)
TukeyHSD(ANOVA_FT)

### QQ-plot for FT & Shapiro Normality (Mixed)

T <- Rep_Mixed %>% select(FT, Generation, Genotype) %>%  arrange(Generation)
T <- T %>% mutate(FT  = FT)
T <- T %>% select(c("Generation", "Genotype", "FT"))

par(mfrow = c(2, 3))
for (i in unique(T$Generation)) {
  p <- T %>% filter(Generation == i) %>% select(FT)
  qqnorm(p$FT, main = paste0("Generation ", i))
  qqline(p$FT)
}

for(i in unique(T$Generation)){
  p <- T %>% filter(Generation == i) %>% select(FT)
  s <- shapiro.test(p$FT)
  print(s)
}

table(T$Generation)

### Kruskal-Wallis Test

kruskal.test(FT ~ Generation, Rep_Mixed)
dunn.test(Rep_Mixed$FT, g = Rep_Mixed$Generation)

### Overlapping Histograms for Flowering Time

ggplot(Rep_Mixed, aes(x = FT, fill = Generation, group = Generation)) +
  geom_histogram(alpha = .5, position = "identity") +
  labs(x = "Average Flowering Time (Days after sowing)",
       y = "Frequency",
       title = "Average Flowering Time Over Generations")



####### Total Seed Weight Distributions and Statistical Analysis

### Histograms for TW (Mixed)

ggplot(Rep_Mixed, aes(x = TOTAL_WEIGHT)) +
  geom_histogram(binwidth = 2) +
  facet_wrap(~Generation, scales = "free_x") +
  labs(x = "Total Seed Weight (g)",
       y = "Frequency",
       title = "Generational Change in Total Seed Weight (g)") +
  scale_x_continuous(breaks = seq(0, 200, 10))

### ANOVA and Homogeneity of Variance for Total Weight

ANOVA_TW <- aov(TOTAL_WEIGHT ~ as.factor(Generation), Rep_Mixed)
summary(ANOVA_TW)
TukeyHSD(ANOVA_TW)
leveneTest(TOTAL_WEIGHT ~ as.factor(Generation), Rep_Mixed)

### QQ-Plot for Total Weight & Shaprio Normality (Mixed)

T <- Rep_Mixed %>% select(TOTAL_WEIGHT, Generation, Genotype) %>%  arrange(Generation)
T <- T %>% mutate(TW  = TOTAL_WEIGHT)
T <- T %>% select(c("Generation", "Genotype", "TW"))

par(mfrow = c(2, 3))
for (i in unique(T$Generation)) {
  p <- T %>% filter(Generation == i) %>% select(TW)
  qqnorm(p$TW, main = paste0("Generation ", i))
  qqline(p$TW)
}

for(i in unique(T$Generation)){
  p <- T %>% filter(Generation == i) %>% select(TW)
  s <- shapiro.test(p$TW)
  print(s)
}


##### Fecundity Distributions and Statistical Analysis

### Histograms for Average Fecundity Over Generations (Mixed)

ggplot(Rep_Mixed, aes(x = Fecundity)) +
  geom_histogram(binwidth = .7) +
  facet_grid(~Generation) +
  stat_bin(bins = 70) +
  labs(x = "Average Fecundity",
       y = "Frequency",
       title = "Generational Change in Average Fecundity")

### Testing for Homogenity of Variance and ANOVA for Average Fecundity

ANOVA_Avg_Fec <- aov(Fecundity ~ as.factor(Generation), Rep_Mixed)
summary(ANOVA_Avg_Fec)
TukeyHSD(ANOVA_Avg_Fec)
leveneTest(Fecundity ~ as.factor(Generation), Rep_Mixed)

### QQ plot Average Fecundity (Mixed)

T <- Rep_Mixed %>% select(Fecundity, Generation, Genotype) %>%  arrange(Generation)
T <- T %>% mutate(fec  = Fecundity)
T <- T %>% select(c("Generation", "Genotype", "fec"))

par(mfrow = c(2, 3))
for (i in unique(T$Generation)) {
  p <- T %>% filter(Generation == i) %>% select(fec)
  qqnorm(p$fec, main = paste0("Generation ", i))
  qqline(p$fec)
}

for(i in unique(T$Generation)){
  p <- T %>% filter(Generation == i) %>% select(fec)
  s <- shapiro.test(p$fec)
  print(s)
}



##### Fitness Distributions & Statistical Analysis

### Histograms for Average Fit Over Generations (Mixed)

ggplot(Rep_Mixed, aes(x = Fitness)) +
  geom_histogram(binwidth = .7) +
  facet_grid(~Generation) +
  stat_bin(bins = 70) +
  labs(x = "Average Fitness",
       y = "Frequency",
       title = "Average Fitness Over Generations")

### Testing for Homogenity of Variance and ANOVA for Average Fitness (Mixed)

ANOVA_Avg_Fit <- aov(Fitness ~ as.factor(Generation), Rep_Mixed)
summary(ANOVA_Avg_Fit)
TukeyHSD(ANOVA_Avg_Fit)
leveneTest(Fitness ~ as.factor(Generation), Rep_Mixed)

### QQ plot Average Fitness (Mixed)

T <- Rep_Mixed %>% select(Fitness, Generation, Genotype) %>%  arrange(Generation)
T <- T %>% mutate(fit  = Fitness)
T <- T %>% select(c("Generation", "Genotype", "fit"))

par(mfrow = c(2, 3))
for (i in unique(T$Generation)) {
  p <- T %>% filter(Generation == i) %>% select(fit)
  qqnorm(p$fit, main = paste0("Generation ", i))
  qqline(p$fit)
}

for(i in unique(T$Generation)){
  p <- T %>% filter(Generation == i) %>% select(fit)
  s <- shapiro.test(p$fit)
  print(s)
}
