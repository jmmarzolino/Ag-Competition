library(tidyverse)
library(readr)
library(dplyr)
library(ggpubr)
library(ggplot2)

library(tidyr)
library(car)
library(gridExtra)
library(dunn.test)

### Average Total Weight Distributions
tmp <- Mixed_2021_2022 %>% group_by(Generation) %>% summarise(Avg_TW = mean(TOTAL_MASS)) %>% ungroup() 
graph_tmp <- full_join(Mixed_2021_2022, tmp, by = "Generation")

ggplot(graph_tmp, aes(x = TOTAL_MASS)) +
  geom_histogram(binwidth = 4) +
  geom_vline(aes(xintercept = Avg_TW, color = "red")) +
  facet_grid(~Generation) +
  scale_x_continuous(breaks = seq(0, 150, 25))

### Homogeneity of Variance and AVOVA/Tukey Post-hoc TW | (Homogeneity of Variance)

ANOVA_TW <- aov(TOTAL_MASS ~ as.factor(Generation), Mixed_2021_2022)
summary(ANOVA_TW)
TukeyHSD(ANOVA_TW)
leveneTest(TOTAL_MASS ~ as.factor(Generation), Mixed_2021_2022)

### QQ-plot for TW & Shapiro Normality | (NOT NORMAL, need Kruskal Wallis and Dunn test)

T <- Mixed_2021_2022 %>% select(TOTAL_MASS, Generation, Genotype) %>%  arrange(Generation)
T <- T %>% mutate(TW = TOTAL_MASS)
T <- T %>% select(c("Generation", "Genotype", "TW"))
T$Genotype <- factor(T$Genotype, levels = T$Genotype[order(T$TW)])

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

### Kruskal Wallis and Dunn Tests

kruskal.test(TOTAL_MASS ~ Generation, Mixed_2021_2022)
dunn.test(Mixed_2021_2022$TOTAL_MASS, Mixed_2021_2022$Generation)





### Average Flowering Time Distributions

tmp <- Mixed_2021_2022 %>% group_by(Generation) %>% summarise(Avg_FT = mean(FT, na.rm = TRUE)) %>% ungroup() 
graph_tmp <- full_join(Mixed_2021_2022, tmp, by = "Generation")

ggplot(graph_tmp, aes(x = FT)) +
  geom_histogram(binwidth = 1.1) +
  geom_vline(aes(xintercept = Avg_FT, color = "red")) +
  facet_grid(~Generation) 

### Homogeneity of Variance and AVOVA/Tukey Post-hoc | (No Homoggeneity of Variance)

ANOVA_FT <- aov(FT ~ as.factor(Generation), Mixed_2021_2022)
summary(ANOVA_FT)
TukeyHSD(ANOVA_FT)
leveneTest(FT ~ as.factor(Generation), Mixed_2021_2022)

### QQ-plot for Single FT & Shapiro Normality | (Not Normal)

T <- Mixed_2021_2022 %>% select(FT, Generation, Genotype) %>%  arrange(Generation)
T <- T %>% mutate(FT = FT)
T <- T %>% select(c("Generation", "Genotype", "FT"))
T$Genotype <- factor(T$Genotype, levels = T$Genotype[order(T$FT)])

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

### Kruskal Wallis and Dunn test

kruskal.test(FT ~ Generation, Mixed_2021_2022)
dunn.test(Mixed_2021_2022$FT, Mixed_2021_2022$Generation)




### Average Fecundity Distributions

tmp <- Mixed_2021_2022 %>% select(c(Generation, Genotype, Fecundity)) %>% group_by(Generation) %>% summarise(Avg_Fec = mean(Fecundity, na.rm = TRUE)) %>% ungroup()
graph_tmp <- full_join(Mixed_2021_2022, tmp, by = "Generation")

ggplot(graph_tmp, aes(x = Fecundity)) +
  geom_histogram() +
  geom_vline(aes(xintercept = Avg_Fec, color = "red")) +
  facet_grid(~Generation) +
  scale_y_continuous(breaks = seq(0, 15, 1))

### Homogeneity of Variance and AVOVA/Tukey Post-hoc

ANOVA_Fec <- aov(Fecundity ~ as.factor(Generation), Mixed_2021_2022)
summary(ANOVA_Fec)
TukeyHSD(ANOVA_Fec)
leveneTest(Fecundity ~ as.factor(Generation), Mixed_2021_2022)

### QQ-plot for Single Fec & Shapiro Normality

T <- Mixed_2021_2022 %>% select(Fecundity, Generation, Genotype) %>%  arrange(Generation)
T <- T %>% mutate(Fec = Fecundity)
T <- T %>% select(c("Generation", "Genotype", "Fec"))
T$Genotype <- factor(T$Genotype, levels = T$Genotype[order(T$Fec)])

par(mfrow = c(2, 3))
for (i in unique(T$Generation)) {
  p <- T %>% filter(Generation == i) %>% select(Fec)
  qqnorm(p$Fec, main = paste0("Generation ", i))
  qqline(p$Fec)
}

for(i in unique(T$Generation)){
  p <- T %>% filter(Generation == i) %>% select(Fec)
  s <- shapiro.test(p$Fec)
  print(s)
}

### Kruskal Wallis and Dunn test

kruskal.test(Fecundity ~ Generation, Mixed_2021_2022)
dunn.test(Mixed_2021_2022$FEC, g = Mixed_2021_2022$Generation)






### Average Fitness Distributions

tmp <- Mixed_2021_2022 %>% group_by(Generation) %>% summarise(Avg_Fit = mean(Fitness, na.rm = TRUE)) %>% ungroup() 
graph_tmp <- full_join(Mixed_2021_2022, tmp, by = "Generation")

ggplot(graph_tmp, aes(x = Fitness)) +
  geom_histogram() +
  geom_vline(aes(xintercept = Avg_Fit, color = "red")) +
  facet_grid(~Generation) +
  scale_x_continuous(breaks = seq(0, 60000, 10000)) +
  theme(axis.text.x = element_text(angle = 45))

### Homogeneity of Variance and AVOVA/Tukey Post-hoc

ANOVA_Fit <- aov(Fitness ~ as.factor(Generation), Mixed_2021_2022)
summary(ANOVA_Fit)
TukeyHSD(ANOVA_Fit)
leveneTest(Fitness ~ as.factor(Generation), Mixed_2021_2022)

### QQ-plot for Fitness & Shapiro Normality

T <- Mixed_2021_2022 %>% select(Fitness, Generation, Genotype) %>%  arrange(Generation)
T <- T %>% mutate(Fit = Fitness)
T <- T %>% select(c("Generation", "Genotype", "Fit"))
T$Genotype <- factor(T$Genotype, levels = T$Genotype[order(T$Fit)])

par(mfrow = c(2, 3))
for (i in unique(T$Generation)) {
  p <- T %>% filter(Generation == i) %>% select(Fit)
  qqnorm(p$Fit, main = paste0("Generation ", i))
  qqline(p$Fit)
}

for(i in unique(T$Generation)){
  p <- T %>% filter(Generation == i) %>% select(Fit)
  s <- shapiro.test(p$Fit)
  print(s)
}

###### Kruskal-Wallis Test and Dunn test

kruskal.test(Fitness ~ Generation, Mixed_2021_2022)
dunn.test(Mixed_2021_2022$Fitness, g = Mixed_2021_2022$Generation)
