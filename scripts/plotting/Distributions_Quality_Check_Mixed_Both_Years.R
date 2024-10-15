#!/usr/bin/env Rscript
#SBATCH --mem=30G
#SBATCH --time=02:00:00
#SBATCH --output=/rhome/jmarz001/bigdata/Ag-Competition/scripts/Distributions_Quality_Check_Mixed_Both_Years.stdout
#SBATCH -p koeniglab

library(tidyverse)


library(ggpubr)



library(car)
library(gridExtra)
library(dunn.test)
library(data.table)

# MIXED 2021-2022 DISTRIBUTIONS
df <- fread("")
# Creating dataframe for distributions

Mixed_2021_2022 <- PHENO_FULL_AVERAGE %>% filter(Condition == "mixed" & Exp_year == 2022)

### Average Total Weight Distributions
graph_tmp <- Mixed_2021_2022 %>% group_by(Generation, Genotype, Exp_year) %>% summarise(Avg_TW = mean(TOTAL_MASS)) %>% ungroup()
tmp <- graph_tmp %>% group_by(Generation) %>% summarise(Generation_Avg = mean(Avg_TW, na.rm = T))
graph_tmp <- full_join(graph_tmp, tmp, by = "Generation")

ggplot(graph_tmp, aes(x = Avg_TW)) +
  geom_histogram(binwidth = 4) +
  geom_vline(aes(xintercept = graph_tmp$Generation_Avg), color = "red") +
  facet_grid(~Generation) +
  scale_x_continuous(breaks = seq(0, 150, 25)) +
  labs(x = "Average Total Weight (grams)")

# Levene Test (checking for homogeneity of variance, assumption of ANOVA) - insignificant means we can assume homogeneity of variance

leveneTest(TOTAL_MASS ~ as.factor(Generation), Mixed_2021_2022)

### QQ-plot & Shapiro Normality for TW - (SOME GROUPS NOT NORMAL, need Kruskal Wallis and Dunn test)

h <- Mixed_2021_2022 %>% select(TOTAL_MASS, Generation, Genotype) %>%  arrange(Generation)
h <- h %>% mutate(TW = TOTAL_MASS)
h <- h %>% select(c("Generation", "Genotype", "TW"))
h$Genotype <- factor(h$Genotype, levels = h$Genotype[order(h$TW)])

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

kruskal.test(TOTAL_MASS ~ Generation, Mixed_2021_2022)
dunn.test(Mixed_2021_2022$TOTAL_MASS, Mixed_2021_2022$Generation)


### Average Flowering Time Distributions

graph_tmp <- Mixed_2021_2022 %>% group_by(Generation, Genotype, Exp_year) %>% summarise(Avg_FT = mean(FT)) %>% ungroup()
tmp <- graph_tmp %>% group_by(Generation) %>% summarise(Generation_Avg = mean(Avg_FT, na.rm = T))
graph_tmp <- full_join(graph_tmp, tmp, by = "Generation")

ggplot(graph_tmp, aes(x = Avg_FT)) +
  geom_histogram(binwidth = 1.1) +
  geom_vline(aes(xintercept = Generation_Avg, color = "red")) +
  facet_grid(~Generation) +
  labs(x = "Average Flowering Time (Days After Sowing)")

# Levene test - unequal variances between groups

leveneTest(FT ~ as.factor(Generation), Mixed_2021_2022)

### QQ-plot for Single FT & Shapiro Normality | (SOME GROUPS NOT NORMAL)

h <- Mixed_2021_2022 %>% select(FT, Generation, Genotype) %>%  arrange(Generation)
h <- h %>% mutate(FT = FT)
h <- h %>% select(c("Generation", "Genotype", "FT"))
h$Genotype <- factor(h$Genotype, levels = h$Genotype[order(h$FT)])

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

kruskal.test(FT ~ Generation, Mixed_2021_2022)
dunn.test(Mixed_2021_2022$FT, Mixed_2021_2022$Generation)

### Average Fecundity Distributions

graph_tmp <- Mixed_2021_2022 %>% group_by(Generation, Genotype, Exp_year) %>% summarise(Avg_Fec = mean(FECUNDITY)) %>% ungroup()
tmp <- graph_tmp %>% group_by(Generation) %>% summarise(Generation_Avg = mean(Avg_Fec, na.rm = T))
graph_tmp <- full_join(graph_tmp, tmp, by = "Generation")

ggplot(graph_tmp, aes(x = Avg_Fec)) +
  geom_histogram(binwidth = .09) +
  geom_vline(aes(xintercept = Generation_Avg, color = "red")) +
  facet_grid(~Generation) +
  scale_y_continuous(breaks = seq(0, 15, 1)) +
  labs(x = "Centered Fecundity")

# Levene Test (EQUAL VAR)

leveneTest(FECUNDITY ~ as.factor(Generation), Mixed_2021_2022)

### QQ-plot for Single Fec & Shapiro Normality (NOT NORMAL)

h <- Mixed_2021_2022 %>% select(FECUNDITY, Generation, Genotype) %>%  arrange(Generation)
h <- h %>% mutate(Fec = FECUNDITY)
h <- h %>% select(c("Generation", "Genotype", "Fec"))
h$Genotype <- factor(h$Genotype, levels = h$Genotype[order(h$Fec)])

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

kruskal.test(FECUNDITY ~ Generation, Mixed_2021_2022)
dunn.test(Mixed_2021_2022$FEC, Mixed_2021_2022$Generation)


### Average Fitness Distributions

graph_tmp <- Mixed_2021_2022 %>% group_by(Generation, Genotype, Exp_year) %>% summarise(Avg_Fit = mean(ABS_FITNESS, na.rm = T)) %>% ungroup()
tmp <- graph_tmp %>% group_by(Generation) %>% summarise(Generation_Avg = mean(Avg_Fit, na.rm = T))
graph_tmp <- full_join(graph_tmp, tmp, by = "Generation")

ggplot(graph_tmp, aes(x = Avg_Fit)) +
  geom_histogram(binwidth = .08) +
  geom_vline(aes(xintercept = Generation_Avg, color = "red")) +
  facet_grid(~Generation) +
  labs(x = "Centered Fitness")

# Levene Test (EQUAL VARIANCE)

leveneTest(ABS_FITNESS ~ as.factor(Generation), Mixed_2021_2022)

### QQ-plot for Fitness & Shapiro Normality (NOT ALL NORMAL)

h <- Mixed_2021_2022 %>% select(ABS_FITNESS, Generation, Genotype) %>%  arrange(Generation)
h <- h %>% mutate(Fit = ABS_FITNESS)
h <- h %>% select(c("Generation", "Genotype", "Fit"))
h$Genotype <- factor(h$Genotype, levels = h$Genotype[order(h$Fit)])

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

kruskal.test(ABS_FITNESS ~ Generation, Mixed_2021_2022)
dunn.test(Mixed_2021_2022$ABS_FITNESS, Mixed_2021_2022$Generation)

# 100 SW Distributions

graph_tmp <- Mixed_2021_2022 %>% group_by(Generation, Genotype, Exp_year) %>% summarise(Avg_100_SW = mean(SEED_WEIGHT_100)) %>% ungroup()
tmp <- graph_tmp %>% group_by(Generation) %>% summarise(Generation_Avg = mean(Avg_100_SW, na.rm = T))
graph_tmp <- full_join(graph_tmp, tmp, by = "Generation")

ggplot(graph_tmp, aes(x = Avg_100_SW)) +
  geom_histogram(binwidth = .07) +
  geom_vline(aes(xintercept = Generation_Avg, color = "red")) +
  facet_grid(~Generation) +
  scale_x_continuous(breaks = seq(0, 60000, 10000)) +
  theme(axis.text.x = element_text(angle = 45)) +
  labs(x = "Average 100 Seed Weight (grams)")

# Levene Test (EQUAL VARIANCE)

leveneTest(SEED_WEIGHT_100 ~ as.factor(Generation), Mixed_2021_2022)

### QQ-plot & Shapiro Normality for 100 SW (NORMAL)

h <- Single_2021_2022 %>% select(SEED_WEIGHT_100, Generation, Genotype) %>%  arrange(Generation)
h <- h %>% mutate(SW = SEED_WEIGHT_100)
h <- h %>% select(c("Generation", "Genotype", "SW"))
h$Genotype <- factor(h$Genotype, levels = h$Genotype[order(h$SW)])

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

ANOVA_100 <- aov(SEED_WEIGHT_100 ~ as.factor(Generation), Mixed_2021_2022)
summary(ANOVA_100)
TukeyHSD(ANOVA_100)


############################



# MIXED 2022-2023 DISTRIBUTIONS

# Creating dataframes for 2022-2023 distributions

Mixed_2022_2023 <- PHENO_FULL_AVERAGE %>% filter(Condition == "mixed" & Exp_year == 2023)

# Average TW Distributions

### Average Total Weight Distributions
graph_tmp <- Mixed_2022_2023 %>% group_by(Generation, Genotype, Exp_year) %>% summarise(Avg_TW = mean(TOTAL_MASS)) %>% ungroup()
tmp <- graph_tmp %>% group_by(Generation) %>% summarise(Generation_Avg = mean(Avg_TW, na.rm = T))
graph_tmp <- full_join(graph_tmp, tmp, by = "Generation")

ggplot(graph_tmp, aes(x = Avg_TW)) +
  geom_histogram(binwidth = 3.5) +
  geom_vline(aes(xintercept = Generation_Avg), color = "red") +
  facet_grid(~Generation) +
  scale_x_continuous(breaks = seq(0, 150, 25)) +
  labs(x = "Average Total Weight (grams)")

# Levene Test (EQUAL VARIANCE)

leveneTest(TOTAL_MASS ~ as.factor(Generation), Single_2022_2023)

### QQ-plot & Shapiro Normality for TW - (NORMAL)

h <- Single_2022_2023 %>% select(TOTAL_MASS, Generation, Genotype) %>%  arrange(Generation)
h <- h %>% mutate(TW = TOTAL_MASS)
h <- h %>% select(c("Generation", "Genotype", "TW"))
h$Genotype <- factor(h$Genotype, levels = h$Genotype[order(h$TW)])

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

ANOVA_TW <- aov(TOTAL_MASS ~ as.factor(Generation), Mixed_2022_2023)
summary(ANOVA_TW)
TukeyHSD(ANOVA_TW)

### Average Flowering Time Distributions

graph_tmp <- Mixed_2022_2023 %>% group_by(Generation, Genotype, Exp_year) %>% summarise(Avg_FT = mean(FT)) %>% ungroup()
tmp <- graph_tmp %>% group_by(Generation) %>% summarise(Generation_Avg = mean(Avg_FT, na.rm = T))
graph_tmp <- full_join(graph_tmp, tmp, by = "Generation")

ggplot(graph_tmp, aes(x = Avg_FT)) +
  geom_histogram(binwidth = 1.1) +
  geom_vline(aes(xintercept = Generation_Avg, color = "red")) +
  facet_grid(~Generation) +
  labs(x = "Average Flowering Time (Days After Sowing)")

# Levene test - unequal variances between groups

leveneTest(FT ~ as.factor(Generation), Mixed_2022_2023)

### QQ-plot for Single FT & Shapiro Normality | (SOME GROUPS NOT NORMAL)

h <- Mixed_2022_2023 %>% select(FT, Generation, Genotype) %>%  arrange(Generation)
h <- h %>% mutate(FT = FT)
h <- h %>% select(c("Generation", "Genotype", "FT"))
h$Genotype <- factor(h$Genotype, levels = h$Genotype[order(h$FT)])

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

kruskal.test(FT ~ Generation, Mixed_2022_2023)
dunn.test(Mixed_2022_2023$FT, Mixed_2022_2023$Generation)


### Average Fecundity Distributions

graph_tmp <- Mixed_2022_2023 %>% group_by(Generation, Genotype, Exp_year) %>% summarise(Avg_Fec = mean(FECUNDITY)) %>% ungroup()
tmp <- graph_tmp %>% group_by(Generation) %>% summarise(Generation_Avg = mean(Avg_Fec, na.rm = T))
graph_tmp <- full_join(graph_tmp, tmp, by = "Generation")

ggplot(graph_tmp, aes(x = Avg_Fec)) +
  geom_histogram(binwidth = .09) +
  geom_vline(aes(xintercept = Generation_Avg, color = "red")) +
  facet_grid(~Generation) +
  scale_y_continuous(breaks = seq(0, 15, 1)) +
  labs(x = "Centered Fecundity")

# Levene Test (EQUAL VARIANCE)

leveneTest(FECUNDITY ~ as.factor(Generation), Mixed_2022_2023)

### QQ-plot for Single Fec & Shapiro Normality (NOT NORMAL)

h <- Mixed_2022_2023 %>% select(FECUNDITY, Generation, Genotype) %>%  arrange(Generation)
h <- h %>% mutate(Fec = FECUNDITY)
h <- h %>% select(c("Generation", "Genotype", "Fec"))
h$Genotype <- factor(h$Genotype, levels = h$Genotype[order(h$Fec)])

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

kruskal.test(FECUNDITY ~ Generation, Mixed_2022_2023)
dunn.test(Mixed_2022_2023$FEC, Mixed_2022_2023$Generation)

### Average Fitness Distributions

graph_tmp <- Mixed_2022_2023 %>% group_by(Generation, Genotype, Exp_year) %>% summarise(Avg_Fit = mean(ABS_FITNESS, na.rm = T)) %>% ungroup()
tmp <- graph_tmp %>% group_by(Generation) %>% summarise(Generation_Avg = mean(Avg_Fit, na.rm = T))
graph_tmp <- full_join(graph_tmp, tmp, by = "Generation")

ggplot(graph_tmp, aes(x = Avg_Fit)) +
  geom_histogram(binwidth = .08) +
  geom_vline(aes(xintercept = Generation_Avg, color = "red")) +
  facet_grid(~Generation) +
  labs(x = "Centered Fitness")

# Levene Test (EQUAL VAR)

leveneTest(ABS_FITNESS ~ as.factor(Generation), Mixed_2022_2023)

### QQ-plot for Fitness & Shapiro Normality (NORMAL)

h <- Mixed_2022_2023 %>% select(ABS_FITNESS, Generation, Genotype) %>%  arrange(Generation)
h <- h %>% mutate(Fit = ABS_FITNESS)
h <- h %>% select(c("Generation", "Genotype", "Fit"))
h$Genotype <- factor(h$Genotype, levels = h$Genotype[order(h$Fit)])

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

ANOVA_Fit <- aov(ABS_FITNESS ~ as.factor(Generation), Mixed_2022_2023)
summary(ANOVA_Fit)
TukeyHSD(ANOVA_Fit)


# 100 SW Distributions

graph_tmp <- Mixed_2022_2023 %>% group_by(Generation, Genotype, Exp_year) %>% summarise(Avg_100_SW = mean(SEED_WEIGHT_100)) %>% ungroup()
tmp <- graph_tmp %>% group_by(Generation) %>% summarise(Generation_Avg = mean(Avg_100_SW, na.rm = T))
graph_tmp <- full_join(graph_tmp, tmp, by = "Generation")

ggplot(graph_tmp, aes(x = Avg_100_SW)) +
  geom_histogram(binwidth = .07) +
  geom_vline(aes(xintercept = Generation_Avg, color = "red")) +
  facet_grid(~Generation) +
  scale_x_continuous(breaks = seq(0, 60000, 10000)) +
  theme(axis.text.x = element_text(angle = 45)) +
  labs(x = "Average 100 Seed Weight (grams)")

# Levene Test (EQUAL VARIANCE) - 100 SW

leveneTest(SEED_WEIGHT_100 ~ as.factor(Generation), Mixed_2022_2023)

### QQ-plot & Shapiro Normality for 100 SW (NOT ALL NORMAL)

h <- Mixed_2022_2023 %>% select(SEED_WEIGHT_100, Generation, Genotype) %>%  arrange(Generation)
h <- h %>% mutate(SW = SEED_WEIGHT_100)
h <- h %>% select(c("Generation", "Genotype", "SW"))
h$Genotype <- factor(h$Genotype, levels = h$Genotype[order(h$SW)])

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

kruskal.test(SEED_WEIGHT_100 ~ Generation, Mixed_2022_2023)
dunn.test(Mixed_2022_2023$SEED_WEIGHT_100, Mixed_2022_2023$Generation)
