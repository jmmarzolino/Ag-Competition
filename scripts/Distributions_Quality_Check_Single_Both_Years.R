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
library(dunn.test)

# SINGLE 2021-2022 DISTRIBUTIONS

# Creating dataframe for distributions

Single_2021_2022 <- PHENO_FULL_AVERAGE %>% filter(Condition == "single" & Exp_year == 2022)

### Average Total Weight Distributions
graph_tmp <- Single_2021_2022 %>% group_by(Generation, Genotypes, Exp_year) %>% summarise(Avg_TW = mean(total_seed_mass_g)) %>% ungroup()
tmp <- graph_tmp %>% group_by(Generation) %>% summarise(Generation_Avg = mean(Avg_TW, na.rm = T))
graph_tmp <- full_join(graph_tmp, tmp, by = 'Generation')

ggplot(graph_tmp, aes(x = Avg_TW)) +
  geom_histogram(binwidth = 4) +
  geom_vline(aes(xintercept = graph_tmp$Generation_Avg), color = 'red') +
  facet_grid(~Generation) +
  scale_x_continuous(breaks = seq(0, 150, 25)) +
  labs(x = 'Average Total Weight (grams)')


# Levene Test (checking for homogeneity of variance, assumption of ANOVA) - insignificant means we can assume homogeneity of variance

leveneTest(total_seed_mass_g ~ as.factor(Generation), Single_2021_2022)

### QQ-plot & Shapiro Normality for TW - (SOME GROUPS NOT NORMAL, need Kruskal Wallis and Dunn test)

h <- Single_2021_2022 %>% select(total_seed_mass_g, Generation, Genotypes) %>%  arrange(Generation)
h <- h %>% mutate(TW = total_seed_mass_g)
h <- h %>% select(c("Generation", "Genotypes", "TW"))
h$Genotypes <- factor(h$Genotypes, levels = h$Genotypes[order(h$TW)])

par(mfrow = c(2, 3))
for (i in unique(h$Generation)) {
  p <- h %>% filter(Generation == i) %>% select(TW)
  qqnorm(p$TW, main = paste0("Generation ", i))
  qqline(p$TW)
}

for(i in unique(h$Generation)){
  p <- h %>% filter(Generation == i) %>% select(TW)
  s <- shapiro.test(p$TW)
  print(s)
}

### Kruskal Wallis and Dunn Tests

kruskal.test(total_seed_mass_g ~ Generation, Single_2021_2022)
dunn.test(Single_2021_2022$total_seed_mass_g, Single_2021_2022$Generation)


### Average Flowering Time Distributions

graph_tmp <- Single_2021_2022 %>% group_by(Generation, Genotypes, Exp_year) %>% summarise(Avg_FT = mean(FT_DAYS)) %>% ungroup()
tmp <- graph_tmp %>% group_by(Generation) %>% summarise(Generation_Avg = mean(Avg_FT, na.rm = T))
graph_tmp <- full_join(graph_tmp, tmp, by = 'Generation')

ggplot(graph_tmp, aes(x = Avg_FT)) +
  geom_histogram(binwidth = 1.1) +
  geom_vline(aes(xintercept = Generation_Avg, color = 'red')) +
  facet_grid(~Generation) +
  labs(x = 'Average Flowering Time (Days After Sowing)')

# Levene test - unequal variances between groups

leveneTest(FT_DAYS ~ as.factor(Generation), Single_2021_2022)

### QQ-plot for Single FT & Shapiro Normality | (SOME GROUPS NOT NORMAL)

h <- Single_2021_2022 %>% select(FT_DAYS, Generation, Genotypes) %>%  arrange(Generation)
h <- h %>% mutate(FT = FT_DAYS)
h <- h %>% select(c("Generation", "Genotypes", "FT"))
h$Genotypes <- factor(h$Genotypes, levels = h$Genotypes[order(h$FT)])

par(mfrow = c(2, 3))
for (i in unique(h$Generation)) {
  p <- h %>% filter(Generation == i) %>% select(FT)
  qqnorm(p$FT, main = paste0("Generation ", i))
  qqline(p$FT)
}

for(i in unique(h$Generation)){
  p <- h %>% filter(Generation == i) %>% select(FT)
  s <- shapiro.test(p$FT)
  print(s)
}

### Kruskal Wallis and Dunn test

kruskal.test(FT_DAYS ~ Generation, Single_2021_2022)
dunn.test(Single_2021_2022$FT_DAYS, Single_2021_2022$Generation)


### Average Fecundity Distributions

graph_tmp <- Single_2021_2022 %>% group_by(Generation, Genotypes, Exp_year) %>% summarise(Avg_Fec = mean(FECUNDITY)) %>% ungroup()
tmp <- graph_tmp %>% group_by(Generation) %>% summarise(Generation_Avg = mean(Avg_Fec, na.rm = T))
graph_tmp <- full_join(graph_tmp, tmp, by = 'Generation')

ggplot(graph_tmp, aes(x = Avg_Fec)) +
  geom_histogram(binwidth = .09) +
  geom_vline(aes(xintercept = Generation_Avg, color = 'red')) +
  facet_grid(~Generation) +
  scale_y_continuous(breaks = seq(0, 15, 1)) +
  labs(x = "Centered Fecundity")

# Levene Test

leveneTest(FECUNDITY ~ as.factor(Generation), Single_2021_2022)

### Homogeneity of Variance and AVOVA/Tukey Post-hoc

ANOVA_Fec <- aov(Fecundity ~ as.factor(Generation), Single_2021_2022)
summary(ANOVA_Fec)
TukeyHSD(ANOVA_Fec)

### QQ-plot for Single Fec & Shapiro Normality (NOT NORMAL)

h <- Single_2021_2022 %>% select(FECUNDITY, Generation, Genotypes) %>%  arrange(Generation)
h <- h %>% mutate(Fec = FECUNDITY)
h <- h %>% select(c("Generation", "Genotypes", "Fec"))
h$Genotypes <- factor(h$Genotypes, levels = h$Genotypes[order(h$Fec)])

par(mfrow = c(2, 3))
for (i in unique(h$Generation)) {
  p <- h %>% filter(Generation == i) %>% select(Fec)
  qqnorm(p$Fec, main = paste0("Generation ", i))
  qqline(p$Fec)
}

for(i in unique(h$Generation)){
  p <- h %>% filter(Generation == i) %>% select(Fec)
  s <- shapiro.test(p$Fec)
  print(s)
}

### Kruskal Wallis and Dunn test

kruskal.test(FECUNDITY ~ Generation, Single_2021_2022)
dunn.test(Single_2021_2022$FECUNDITY, Single_2021_2022$Generation)


### Average Fitness Distributions

graph_tmp <- Single_2021_2022 %>% group_by(Generation, Genotypes, Exp_year) %>% summarise(Avg_Fit = mean(ABS_FITNESS, na.rm = T)) %>% ungroup()
tmp <- graph_tmp %>% group_by(Generation) %>% summarise(Generation_Avg = mean(Avg_Fit, na.rm = T))
graph_tmp <- full_join(graph_tmp, tmp, by = 'Generation')

ggplot(graph_tmp, aes(x = Avg_Fit)) +
  geom_histogram(binwidth = .08) +
  geom_vline(aes(xintercept = Generation_Avg, color = 'red')) +
  facet_grid(~Generation) +
  labs(x = 'Centered Fitness')

# Levene Test

leveneTest(ABS_FITNESS ~ as.factor(Generation), Single_2021_2022)

### QQ-plot for Fitness & Shapiro Normality (NOT ALL NORMAL)

h <- Single_2021_2022 %>% select(ABS_FITNESS, Generation, Genotypes) %>%  arrange(Generation)
h <- h %>% mutate(Fit = ABS_FITNESS)
h <- h %>% select(c("Generation", "Genotypes", "Fit"))
h$Genotypes <- factor(h$Genotypes, levels = h$Genotypes[order(h$Fit)])

par(mfrow = c(2, 3))
for (i in unique(h$Generation)) {
  p <- h %>% filter(Generation == i) %>% select(Fit)
  qqnorm(p$Fit, main = paste0("Generation ", i))
  qqline(p$Fit)
}

for(i in unique(h$Generation)){
  p <- h %>% filter(Generation == i) %>% select(Fit)
  s <- shapiro.test(p$Fit)
  print(s)
}

# Kruskal-Wallis Test and Dunn test

kruskal.test(ABS_FITNESS ~ Generation, Single_2021_2022)
dunn.test(Single_2021_2022$ABS_FITNESS, Single_2021_2022$Generation)

# 100 SW Distributions

graph_tmp <- Single_2021_2022 %>% group_by(Generation, Genotypes, Exp_year) %>% summarise(Avg_100_SW = mean(`100_seed_weight`)) %>% ungroup()
tmp <- graph_tmp %>% group_by(Generation) %>% summarise(Generation_Avg = mean(Avg_100_SW, na.rm = T))
graph_tmp <- full_join(graph_tmp, tmp, by = 'Generation')

ggplot(graph_tmp, aes(x = Avg_100_SW)) +
  geom_histogram(binwidth = .07) +
  geom_vline(aes(xintercept = Generation_Avg, color = 'red')) +
  facet_grid(~Generation) +
  scale_x_continuous(breaks = seq(0, 60000, 10000)) +
  theme(axis.text.x = element_text(angle = 45)) +
  labs(x = 'Average 100 Seed Weight (grams)')

# Levene Test (EQUAL VARIANCE)

leveneTest(`100_seed_weight` ~ as.factor(Generation), Single_2021_2022)

### QQ-plot & Shapiro Normality for 100 SW (NORMAL)

h <- Single_2021_2022 %>% select(`100_seed_weight`, Generation, Genotypes) %>%  arrange(Generation)
h <- h %>% mutate(SW = `100_seed_weight`)
h <- h %>% select(c("Generation", "Genotypes", "SW"))
h$Genotypes <- factor(h$Genotypes, levels = h$Genotypes[order(h$SW)])

par(mfrow = c(2, 3))
for (i in unique(h$Generation)) {
  p <- h %>% filter(Generation == i) %>% select(SW)
  qqnorm(p$SW, main = paste0("Generation ", i))
  qqline(p$SW)
}

for(i in unique(h$Generation)){
  p <- h %>% filter(Generation == i) %>% select(SW)
  s <- shapiro.test(p$SW)
  print(s)
}

# ANOVA & TUKEY for 100 SW

ANOVA_100 <- aov(`100_seed_weight` ~ as.factor(Generation), Single_2021_2022)
summary(ANOVA_100)
TukeyHSD(ANOVA_100)


############################



# SINGLE 2022-2023 DISTRIBUTIONS

# Creating dataframes for 2022-2023 distributions

Single_2022_2023 <- PHENO_FULL_AVERAGE %>% filter(Condition == "single" & Exp_year == 2023)

# Average TW Distributions

### Average Total Weight Distributions
graph_tmp <- Single_2022_2023 %>% group_by(Generation, Genotypes, Exp_year) %>% summarise(Avg_TW = mean(total_seed_mass_g)) %>% ungroup()
tmp <- graph_tmp %>% group_by(Generation) %>% summarise(Generation_Avg = mean(Avg_TW, na.rm = T))
graph_tmp <- full_join(graph_tmp, tmp, by = 'Generation')

ggplot(graph_tmp, aes(x = Avg_TW)) +
  geom_histogram(binwidth = 3.5) +
  geom_vline(aes(xintercept = graph_tmp$Generation_Avg), color = 'red') +
  facet_grid(~Generation) +
  scale_x_continuous(breaks = seq(0, 150, 25)) +
  labs(x = 'Average Total Weight (grams)')

# Levene Test (EQUAL VARIANCE)

leveneTest(total_seed_mass_g ~ as.factor(Generation), Single_2022_2023)

### QQ-plot & Shapiro Normality for TW - (NORMAL)

h <- Single_2022_2023 %>% select(total_seed_mass_g, Generation, Genotypes) %>%  arrange(Generation)
h <- h %>% mutate(TW = total_seed_mass_g)
h <- h %>% select(c("Generation", "Genotypes", "TW"))
h$Genotypes <- factor(h$Genotypes, levels = h$Genotypes[order(h$TW)])

par(mfrow = c(2, 3))
for (i in unique(h$Generation)) {
  p <- h %>% filter(Generation == i) %>% select(TW)
  qqnorm(p$TW, main = paste0("Generation ", i))
  qqline(p$TW)
}

for(i in unique(h$Generation)){
  p <- h %>% filter(Generation == i) %>% select(TW)
  s <- shapiro.test(p$TW)
  print(s)
}

### ANOVA and TUKEY TW

ANOVA_TW <- aov(total_seed_mass_g ~ as.factor(Generation), Single_2022_2023)
summary(ANOVA_TW)
TukeyHSD(ANOVA_TW)

### Average Flowering Time Distributions

graph_tmp <- Single_2022_2023 %>% group_by(Generation, Genotypes, Exp_year) %>% summarise(Avg_FT = mean(FT_DAYS)) %>% ungroup()
tmp <- graph_tmp %>% group_by(Generation) %>% summarise(Generation_Avg = mean(Avg_FT, na.rm = T))
graph_tmp <- full_join(graph_tmp, tmp, by = 'Generation')

ggplot(graph_tmp, aes(x = Avg_FT)) +
  geom_histogram(binwidth = 1.1) +
  geom_vline(aes(xintercept = Generation_Avg, color = 'red')) +
  facet_grid(~Generation) +
  labs(x = 'Average Flowering Time (Days After Sowing)')

# Levene test - unequal variances between groups

leveneTest(FT_DAYS ~ as.factor(Generation), Single_2022_2023)

### QQ-plot for Single FT & Shapiro Normality | (SOME GROUPS NOT NORMAL)

h <- Single_2022_2023 %>% select(FT_DAYS, Generation, Genotypes) %>%  arrange(Generation)
h <- h %>% mutate(FT = FT_DAYS)
h <- h %>% select(c("Generation", "Genotypes", "FT"))
h$Genotypes <- factor(h$Genotypes, levels = h$Genotypes[order(h$FT)])

par(mfrow = c(2, 3))
for (i in unique(h$Generation)) {
  p <- h %>% filter(Generation == i) %>% select(FT)
  qqnorm(p$FT, main = paste0("Generation ", i))
  qqline(p$FT)
}

for(i in unique(h$Generation)){
  p <- h %>% filter(Generation == i) %>% select(FT)
  s <- shapiro.test(p$FT)
  print(s)
}

### Kruskal Wallis and Dunn test

kruskal.test(FT_DAYS ~ Generation, Single_2022_2023)
dunn.test(Single_2022_2023$FT_DAYS, Single_2022_2023$Generation)


### Average Fecundity Distributions

graph_tmp <- Single_2022_2023 %>% group_by(Generation, Genotypes, Exp_year) %>% summarise(Avg_Fec = mean(FECUNDITY)) %>% ungroup()
tmp <- graph_tmp %>% group_by(Generation) %>% summarise(Generation_Avg = mean(Avg_Fec, na.rm = T))
graph_tmp <- full_join(graph_tmp, tmp, by = 'Generation')

ggplot(graph_tmp, aes(x = Avg_Fec)) +
  geom_histogram(binwidth = .09) +
  geom_vline(aes(xintercept = Generation_Avg, color = 'red')) +
  facet_grid(~Generation) +
  scale_y_continuous(breaks = seq(0, 15, 1)) +
  labs(x = "Centered Fecundity")

# Levene Test (EQUAL VARIANCE)

leveneTest(FECUNDITY ~ as.factor(Generation), Single_2022_2023)

### QQ-plot for Single Fec & Shapiro Normality (NOT NORMAL)

h <- Single_2022_2023 %>% select(FECUNDITY, Generation, Genotypes) %>%  arrange(Generation)
h <- h %>% mutate(Fec = FECUNDITY)
h <- h %>% select(c("Generation", "Genotypes", "Fec"))
h$Genotypes <- factor(h$Genotypes, levels = h$Genotypes[order(h$Fec)])

par(mfrow = c(2, 3))
for (i in unique(h$Generation)) {
  p <- h %>% filter(Generation == i) %>% select(Fec)
  qqnorm(p$Fec, main = paste0("Generation ", i))
  qqline(p$Fec)
}

for(i in unique(h$Generation)){
  p <- h %>% filter(Generation == i) %>% select(Fec)
  s <- shapiro.test(p$Fec)
  print(s)
}

### Kruskal Wallis and Dunn test

kruskal.test(FECUNDITY ~ Generation, Single_2022_2023)
dunn.test(Single_2022_2023$FECUNDITY, Single_2022_2023$Generation)

### Average Fitness Distributions

graph_tmp <- Single_2022_2023 %>% group_by(Generation, Genotypes, Exp_year) %>% summarise(Avg_Fit = mean(ABS_FITNESS, na.rm = T)) %>% ungroup()
tmp <- graph_tmp %>% group_by(Generation) %>% summarise(Generation_Avg = mean(Avg_Fit, na.rm = T))
graph_tmp <- full_join(graph_tmp, tmp, by = 'Generation')

ggplot(graph_tmp, aes(x = Avg_Fit)) +
  geom_histogram(binwidth = .08) +
  geom_vline(aes(xintercept = Generation_Avg, color = 'red')) +
  facet_grid(~Generation) +
  labs(x = 'Centered Fitness')

# Levene Test (EQUAL VAR)

leveneTest(ABS_FITNESS ~ as.factor(Generation), Single_2022_2023)

### QQ-plot for Fitness & Shapiro Normality (NORMAL)

h <- Single_2022_2023 %>% select(ABS_FITNESS, Generation, Genotypes) %>%  arrange(Generation)
h <- h %>% mutate(Fit = ABS_FITNESS)
h <- h %>% select(c("Generation", "Genotypes", "Fit"))
h$Genotypes <- factor(h$Genotypes, levels = h$Genotypes[order(h$Fit)])

par(mfrow = c(2, 3))
for (i in unique(h$Generation)) {
  p <- h %>% filter(Generation == i) %>% select(Fit)
  qqnorm(p$Fit, main = paste0("Generation ", i))
  qqline(p$Fit)
}

for(i in unique(h$Generation)){
  p <- h %>% filter(Generation == i) %>% select(Fit)
  s <- shapiro.test(p$Fit)
  print(s)
}

# ANOVA AND TUKEY (FITNESS)

ANOVA_Fit <- aov(ABS_FITNESS ~ as.factor(Generation), Single_2022_2023)
summary(ANOVA_Fit)
TukeyHSD(ANOVA_Fit)


# 100 SW Distributions

graph_tmp <- Single_2022_2023 %>% group_by(Generation, Genotypes, Exp_year) %>% summarise(Avg_100_SW = mean(`100_seed_weight`)) %>% ungroup()
tmp <- graph_tmp %>% group_by(Generation) %>% summarise(Generation_Avg = mean(Avg_100_SW, na.rm = T))
graph_tmp <- full_join(graph_tmp, tmp, by = 'Generation')

ggplot(graph_tmp, aes(x = Avg_100_SW)) +
  geom_histogram(binwidth = .07) +
  geom_vline(aes(xintercept = Generation_Avg, color = 'red')) +
  facet_grid(~Generation) +
  scale_x_continuous(breaks = seq(0, 60000, 10000)) +
  theme(axis.text.x = element_text(angle = 45)) +
  labs(x = 'Average 100 Seed Weight (grams)')

# Levene Test (EQUAL VARIANCE) - 100 SW

leveneTest(`100_seed_weight` ~ as.factor(Generation), Single_2022_2023)

### QQ-plot & Shapiro Normality for 100 SW (NORMAL, for the most part)

h <- Single_2022_2023 %>% select(`100_seed_weight`, Generation, Genotypes) %>%  arrange(Generation)
h <- h %>% mutate(SW = `100_seed_weight`)
h <- h %>% select(c("Generation", "Genotypes", "SW"))
h$Genotypes <- factor(h$Genotypes, levels = h$Genotypes[order(h$SW)])

par(mfrow = c(2, 3))
for (i in unique(h$Generation)) {
  p <- h %>% filter(Generation == i) %>% select(SW)
  qqnorm(p$SW, main = paste0("Generation ", i))
  qqline(p$SW)
}

for(i in unique(h$Generation)){
  p <- h %>% filter(Generation == i) %>% select(SW)
  s <- shapiro.test(p$SW)
  print(s)
}

# Checking the number of samples per group

table(Single_2022_2023$Generation)

# Kruskal Wallis and Dunn - 100 SW

kruskal.test(`100_seed_weight` ~ Generation, Single_2022_2023)
dunn.test(Single_2022_2023$`100_seed_weight`, Single_2022_2023$Generation)
