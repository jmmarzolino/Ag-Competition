library(tidyverse)
library(readr)
library(dplyr)
library(ggpubr)
library(ggplot2)

library(tidyr)
library(car)
library(gridExtra)
library(dunn.test)

### Load Data

Rep_Single <- read_delim("~/Documents/GitHub/Ag-Competition/Rep_Single")

##### Flowering Time Distributions and Statistical Analysis

### Histograms for Average Flowering Time Over Generations Single

ggplot(Rep_Single, aes(x = FT)) +
  geom_histogram(binwidth = .9) +
  facet_grid(~Generation) +
  stat_bin(bins = 50) +
  labs(x = "Average Flowering Time (Days)",
       y = "Frequency",
       title = "Evolution of Flowering Time Over Generations")


### Testing for Homogeneity of Variance and ANOVA for FT

ANOVA_FT <- aov(FT ~ as.factor(Generation), Rep_Single)
summary(ANOVA_FT)
TukeyHSD(ANOVA_FT)
leveneTest(FT ~ as.factor(Generation), Rep_Single)

### QQ-plot for Single FT & Shapiro Normality

T <- Rep_Single %>% select(FT, Generation, Genotype) %>%  arrange(Generation)
T <- T %>% mutate(FT  = FT)
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

table(T$Generation)

### Bar Graphs Single For Each Generation (FT)

for (i in unique(T$Generation)){
  tmp <- ggplot(subset(T, Generation == i), aes(Genotype, FT)) + geom_bar(stat = "identity") + 
    labs(y = "Flowering Time (Days After Sowing)") +
    coord_flip() +
    ggtitle(paste0("Generation ", i)) 
  ggsave(tmp, file = paste0("scripts/plotting/Bar_Plot_FT_Generation_", i, ".png"))
}

### Kruskal-Wallis Test

kruskal.test(FT ~ Generation, Rep_Single)
dunn.test(Rep_Single$FT, g = Rep_Single$Generation)

### Overlapping Histograms for Flowering Time

ggplot(Rep_Single, aes(x = FT, fill = Generation, group = Generation)) +
  geom_histogram(alpha = .5, position = "identity") +
  labs(x = "Average Flowering Time (Days)",
       y = "Frequency",
       title = "Average Flowering Time Over Generations")


##### Total Seed Weight Distributions and Statistical Analysis

### Histograms for Average Total Weight Over Generations (Single) - Normal

Rep_Single <- Rep_Single %>% group_by(Generation) %>% mutate(GenAvg = mean(`Brown Bag Weight`)) %>% ungroup()
ggplot(Rep_Single, aes(x = `Brown Bag Weight`)) +
  geom_histogram(binwidth = 1) +
  facet_grid(~Generation) +
  stat_bin(bins = 60) +
  labs(x = "Average Total Weight (g)",
       y = "Frequency",
       title = "Average Total Weight Over Generations")+
  geom_vline(aes(xintercept = GenAvg), color = "red") 

### Testing for Homogenity of Variance and ANOVA for Average Total Weight

ANOVA_Total_Weight <- aov(`Brown Bag Weight` ~ as.factor(Generation), Rep_Single)
summary(ANOVA_Total_Weight)
TukeyHSD(ANOVA_Total_Weight)
leveneTest(`Brown Bag Weight` ~ as.factor(Generation), Rep_Single)

### QQ plot for Total Weight (Single)

T <- Rep_Single %>% select(`Brown Bag Weight`, Generation, Genotype) %>%  arrange(Generation)
T <- T %>% mutate(TW  = `Brown Bag Weight`)
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

### Bar Graphs Single For Each Generation (Total Weight)

for (i in unique(T$Generation)){
  tmp <- ggplot(subset(T, Generation == i), aes(Genotype, TW)) + geom_bar(stat = "identity") + 
    labs(y = "Total Weight (g)") +
    coord_flip() +
    ggtitle(paste0("Generation ", i)) 
  ggsave(tmp, file = paste0("scripts/plotting/Bar_Plot_Total_Weight_Generation_", i, ".png"))
}


##### Fecundity Distributions and Statistical Analysis

### Histograms for Average Fecundity Over Generations (Single) - Normal

Rep_Single <- Rep_Single %>% group_by(Generation) %>% mutate(GenAvg = mean(Fecundity)) %>% ungroup()
ggplot(Rep_Single, aes(x = Fecundity)) +
  geom_histogram(binwidth = .9) +
  facet_grid(~Generation) +
  stat_bin(bins = 70) +
  labs(x = "Average Fecundity",
       y = "Frequency",
       title = "Average Fecundity Over Generations") +
  geom_vline(aes(xintercept = GenAvg), color = "red") 

### Testing for Homogenity of Variance and ANOVA for Average Fecundity

ANOVA_Avg_Fec <- aov(Fecundity ~ as.factor(Generation), Rep_Single)
summary(ANOVA_Avg_Fec)
TukeyHSD(ANOVA_Avg_Fec)
leveneTest(Fecundity ~ as.factor(Generation), Rep_Single)

### QQ plot Average Fecundity (Single)

T <- Rep_Single %>% select(Fecundity, Generation, Genotype) %>%  arrange(Generation)
T <- T %>% mutate(fec  = Fecundity)
T <- T %>% select(c("Generation", "Genotype", "fec"))
T$Genotype <- factor(T$Genotype, levels = T$Genotype[order(T$fec)])

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

### Bar Graphs Single For Each Generation (Fecundity)

for (i in unique(T$Generation)){
  tmp <- ggplot(subset(T, Generation == i), aes(Genotype, fec)) + geom_bar(stat = "identity") + 
    labs(y = "Fecundity") +
    coord_flip() +
    ggtitle(paste0("Generation ", i)) 
  ggsave(tmp, file = paste0("scripts/plotting/Bar_Plot_Fecundity_Generation_", i, ".png"))
}

##### Fitness Distributions & Statistical Analysis

### Histograms for Average Fit Over Generations (Single) - Normal

Rep_Single <- Rep_Single %>% group_by(Generation) %>% mutate(GenAvg = mean(Fitness)) %>% ungroup()
ggplot(Rep_Single, aes(x = Fitness)) +
  geom_histogram(binwidth = .9) +
  facet_grid(~Generation) +
  stat_bin(bins = 70) +
  labs(x = "Average Fitness",
       y = "Frequency",
       title = "Average Fitness Over Generations") +
  geom_vline(aes(xintercept = GenAvg), color = "red")

### Testing for Homogenity of Variance and ANOVA for Average Fitness

ANOVA_Avg_Fit <- aov(Fitness ~ as.factor(Generation), Rep_Single)
summary(ANOVA_Avg_Fit)
TukeyHSD(ANOVA_Avg_Fit)
leveneTest(Fitness ~ as.factor(Generation), Rep_Single)

### QQ plot Average Fitness (Single)

T <- Rep_Single %>% select(Fitness, Generation, Genotype) %>%  arrange(Generation)
T <- T %>% mutate(fit  = Fitness)
T <- T %>% select(c("Generation", "Genotype", "fit"))
T$Genotype <- factor(T$Genotype, levels = T$Genotype[order(T$fit)])

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

### Bar Graphs Single For Each Generation (Fitness)

for (i in unique(T$Generation)){
  tmp <- ggplot(subset(T, Generation == i), aes(Genotype, fit)) + geom_bar(stat = "identity") + 
    labs(y = "Fitness") +
    coord_flip() +
    ggtitle(paste0("Generation ", i)) 
  ggsave(tmp, file = paste0("scripts/plotting/Bar_Plot_Fitness_Generation_", i, ".png"))
}
