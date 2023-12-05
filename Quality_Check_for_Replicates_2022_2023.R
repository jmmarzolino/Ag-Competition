library(tidyverse)
library(readr)
library(dplyr)
library(ggpubr)
library(ggplot2)
library(googlesheets4)
library(tidyr)
library(car)
library(gridExtra)
library(dunn.test)


### Load Data

Full_Data <- read_delim("~/Documents/GitHub/Ag-Competition/Full_Data")
Side_By_Side_Replicates <- read_delim("~/Documents/GitHub/Ag-Competition/Side_by_Side_Replicates")
Rep_Mixed <- read_delim("~/Documents/GitHub/Ag-Competition/Rep_Mixed")

### Correlation of Total Seed Weight among replicates

ggplot(Side_By_Side_Replicates, aes(x=`Brown Bag Weight`, y=`Brown Bag Weight__2`, add = "reg.line")) +
  geom_point(alpha = .3) +
  geom_abline(slope =1, intercept= 0, color = "red") +
  stat_cor(label.y = 235, label.x = 20) +
  geom_smooth(method = "lm") +
  labs(x = "Total Seed Weight Replicate 1 (grams)",
       y= "Total Seed Weight Replicate 2 (grams)",
       title = "Correlation of Replicates for Total Seed Weight") +
  facet_wrap(~Condition)

ggsave("Correlation_of_Rep1_and_Rep2_Total_SW.png")

res <- Side_By_Side_Replicates %>% summarise(residuals = resid(lm(`Brown Bag Weight__2` ~ `Brown Bag Weight`)))
Side_By_Side_Replicates$Residuals <- res$residuals

### Correlation of Fecundity by Replicate

ggplot(Side_By_Side_Replicates, aes(x=Fecundity, y= Fecundity__2)) +
  geom_point(alpha = .5) +
  stat_cor(label.y = 4200, label.x = 20) +
  geom_abline(slope =1, intercept= 0, color = "red") +
  geom_smooth(method = "lm") +
  labs(x = "Fecundity (replicate 1)",
       y= "Fecundity (replicate 2)",
       title = "Correlation of Replicates for Fecundity") +
  facet_wrap(~Condition)

ggsave("Correlation_of_Rep1_and_Rep2_Fecundity.png")

res <- Side_By_Side_Replicates %>% summarise(residuals = resid(lm(Fecundity__2 ~ Fecundity)))
Side_By_Side_Replicates$Residuals <- res$residuals

### Correlation of Flowering Time

ggplot(Side_By_Side_Replicates, aes(FT_DAYS, FT_DAYS__2)) +
  geom_jitter(alpha = .5) +
  geom_smooth(method = 'lm') +
  stat_cor(label.x = 100, label.y = 135) +
  geom_abline(slope = 1, intercept = 0, col = "red") +
  labs(x = "Flowering Time Replicate 1 (Days)",
       y = "Flowering Time Replicate 2 (Days)",
       title = "Correlation of Replicates for Flowering Time") +
  facet_wrap(~Condition)

### Correlation of Fitness

ggplot(Side_By_Side_Replicates, aes(Fitness, Fitness__2)) +
  geom_jitter(alpha = .5) +
  geom_smooth(method = 'lm') +
  stat_cor(label.x = 100, label.y = 40000) +
  geom_abline(slope = 1, intercept = 0, col = "red") +
  labs(x = "Fitness Replicate 1",
       y = "Fitness Replicate 2",
       title = "Correlation of Replicates for Fitness") +
  facet_wrap(~Condition)

res <- Side_By_Side_Replicates %>% summarise(residuals = resid(lm(Fitness__2 ~ Fitness)))
Side_By_Side_Replicates$Residuals <- res$residuals

ggplot(SW100, aes(Fitness, Residuals)) +
  geom_point(alpha = .5) +
  geom_hline(yintercept = 0, col = "red") +
  labs(title = "Residuals for Replicates of Fitness (Mixed)") 

### Correlation of 100 Seed Weight between Replicates (Single)

r1_100SW <- subset(Full_Data, replicate == "rep 1" & Condition == "single")
r2_100SW <- subset(Full_Data, replicate == "rep 2" & Condition == "single")
colnames(r2_100SW)[1:13] <- paste(colnames(r2_100SW)[c(1:13)], '_2', sep = '_')
SW100 <- inner_join(r1_100SW, r2_100SW, by = c("Genotypes" = "Genotypes__2","Condition" = "Condition__2"))

ggplot(SW100, aes(x=`Brown Bag Weight`, y=`Brown Bag Weight__2`)) +
  geom_jitter(alpha = .5) +
  stat_cor(label.y = 235, label.x = 20) +
  geom_abline(slope =1, intercept= 0, color = "darkgreen") +
  geom_smooth(method = "lm") +
  labs(x = "100 Seed Weight (rep 1)",
       y= "100 Seed Weight (rep 2)",
       title = "Correlation of Single Replicate for 100 Seed Weight")

ggsave("Correlation_100_SW_Single.png")

res <- SW100 %>% summarise(residuals = resid(lm(`Brown Bag Weight__2` ~ `Brown Bag Weight`)))
SW100$Residuals <- res$residuals

### Correlation of 100 Seed Weight by Replicates (Mixed)

r1_100SW_mixed<- subset(Full_Data, replicate == "rep 1" & Condition == "mixed")
r2_100SW_mixed <- subset(Full_Data, replicate == "rep 2" & Condition == "mixed")
colnames(r2_100SW_mixed)[1:13] <- paste(colnames(r2_100SW_mixed)[c(1:13)], '_2', sep = '_')
SW100 <- inner_join(r1_100SW_mixed, r2_100SW_mixed, by = c("Genotypes" = "Genotypes__2","Condition" = "Condition__2"))

y <- ggplot(SW100, aes(x=`Brown Bag Weight`, y=`Brown Bag Weight__2`)) +
  geom_jitter(alpha = .5) +
  stat_cor(label.y = 235, label.x = 20) +
  geom_abline(slope =1, intercept= 0, color = "darkgreen") +
  geom_smooth(method = "lm") +
  labs(x = "100 Seed Weight (rep 1)",
       y= "100 Seed Weight (rep 2)",
       title = "Mixed 100 Seed Weight Correlation by Replicate")

ggsave("Correlation_100_SW_Mixed.png")

res <- SW100 %>% summarise(residuals = resid(lm(`Brown Bag Weight__2` ~ `Brown Bag Weight`)))
SW100$Residuals <- res$residuals

l <- ggplot(SW100, aes(`100 seed weight`, Residuals)) +
  geom_point(alpha = .5) +
  geom_hline(yintercept = 0, col = "red") +
  labs(title = "Residuals for 100 Seed Weight (Mixed)")

grid.arrange(y, l, ncol = 2, nrow = 1) +
ggsave("Correlation_100_SW_By_Replicates_Mixed.png")

### Q-Q Plot For Total Weight ~ Generation

T <- Full_Data %>% select(`Brown Bag Weight`, Generation, Genotypes) %>%  arrange(Generation)
T <- T %>% mutate(Total_Weight = `Brown Bag Weight`)
T <- T %>% select(c("Generation", "Genotypes", "Total_Weight"))


par(mfrow = c(2, 3))
for (i in unique(T$Generation)) {
  p <- T %>% filter(Generation == i) %>% select(Total_Weight)
  qqnorm(p$Total_Weight, main = paste0("Total Weight for Generation ", i))
  qqline(p$Total_Weight)
}

for(i in unique(T$Generation)){
  p <- T %>% filter(Generation == i) %>% select(Total_Weight)
  s <- shapiro.test(p$Total_Weight)
  print(s)
}

### Histograms For Total Weight by Generation

ggplot(Full_Data, aes(x= `Brown Bag Weight`)) +
  geom_histogram(binwidth = 1) +
  stat_bin(bins = 50) +
  labs(x= "Total Weight (grams)",
       y= "Frequency",
       title = "Total Weight over Generations") +
  facet_grid(~Generation)

ANOVA_Total_weight_over_Generations <- aov(`Brown Bag Weight` ~ as.factor(Generation), Full_Data)
summary(ANOVA_Total_weight_over_Generations)
TukeyHSD(ANOVA_Total_weight_over_Generations)
leveneTest(`Brown Bag Weight` ~ as.factor(Generation), Full_Data)

ggsave("Histograms_Total_Weight_Over_Generations.png")

### Overlapping Histograms For Average Total Weight by Generation

Side_By_Side_Replicates$Generation <- as.factor(Side_By_Side_Replicates$Generation)
ggplot(Side_By_Side_Replicates, aes(x = Avg_Total_Weight, group = Generation, fill = Generation)) +
  geom_histogram(alpha = .5, binwidth = .5) +
  stat_bin(bins = 50) +
  scale_fill_brewer(palette = "Blues")

ggsave("Overlapping_Histograms_Total_Weight_Over_Generations.png")

### Histograms For 100 Seed Weight by Generation

ggplot(Side_By_Side_Replicates, aes(x = `Avg_100_SW`)) +
  geom_histogram(binwidth = .1) +
  stat_bin(bins = 50) +
  labs(x = "100 Seed Weight (grams)",
       y = "frequency",
       title = "100 Seed Weight by Generations") +
  facet_grid(~Generation)

### Testing for Homogenity of Variance and ANOVA for 100 SW

ANOVA_100_SW <- aov(`100 seed weight` ~ as.factor(Generation), Side_By_Side_Replicates)
summary(ANOVA_100_SW)
TukeyHSD(ANOVA_100_SW)
leveneTest(`100 seed weight` ~ as.factor(Generation), Side_By_Side_Replicates)

ggsave("Histograms_Average_SW_100_Over_Generations.png")

### Overlapping Histogram for 100 Seed Weight by Generation

Full_Data$Generation <- as.factor(Full_Data$Generation)
ggplot(Full_Data, aes(x = `100 seed weight`, group = Generation, fill = Generation)) +
  geom_histogram(alpha = .5, position = 'identity') +
  scale_fill_brewer(palette = "Greens") +
  labs(x = "100 Seed Weight (grams)",
       y = "frequency",
       title = "100 Seed Weight by Generations")

ggsave("Overplapping_Histograms_SW_100_Over_Generations.png")




### Histograms for Average Flowering Time Over Generations Single

ggplot(Rep_Single, aes(x = Avg_FT)) +
  geom_histogram(binwidth = .9) +
  facet_grid(~Generation) +
  stat_bin(bins = 50) +
  labs(x = "Average Flowering Time (Days)",
       y = "Frequency",
       title = "Average Flowering Time Over Generations (Single)")

### Testing for Homogenity of Variance and ANOVA for FT

ANOVA_FT <- aov(Avg_FT ~ as.factor(Generation), Rep_Single)
summary(ANOVA_FT)
TukeyHSD(ANOVA_FT)
leveneTest(Avg_FT ~ as.factor(Generation), Side_By_Side_Replicates)

### QQ-plot for Single FT
T <- Rep_Single %>% select(Avg_FT, Generation, Genotypes) %>%  arrange(Generation)
T <- T %>% mutate(FT  = Avg_FT)
T <- T %>% select(c("Generation", "Genotypes", "FT"))

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

kruskal.test(Avg_FT ~ Generation, Side_By_Side_Replicates)
dunn.test(Side_By_Side_Replicates$Avg_FT, g = Side_By_Side_Replicates$Generation)

### Overlapping Histograms for Flowering Time

ggplot(Rep_Single, aes(x = Avg_FT, fill = Generation, group = Generation)) +
  geom_histogram(alpha = .5, position = 'identity') +
  labs(x = "Average Flowering Time (Days)",
       y = "Frequency",
       title = "Average Flowering Time Over Generations")

### Histograms for Average Total Weight Over Generations Single

ggplot(Rep_Single, aes(x = Avg_Total_Weight)) +
  geom_histogram(binwidth = .9) +
  facet_grid(~Generation) +
  stat_bin(bins = 60) +
  labs(x = "Average Total Weight (g)",
       y = "Frequency",
       title = "Average Total Weight Over Generations")


### Testing for Homogenity of Variance and ANOVA for Average Total Weight

ANOVA_Total_Weight <- aov(Avg_Total_Weight ~ as.factor(Generation), Rep_Single)
summary(ANOVA_Total_Weight)
TukeyHSD(ANOVA_Total_Weight)
leveneTest(Avg_Total_Weight ~ as.factor(Generation), Rep_Single)

### QQ plot for Total Weight Single

T <- Rep_Single %>% select(Avg_Total_Weight, Generation, Genotypes) %>%  arrange(Generation)
T <- T %>% mutate(TW  = Avg_Total_Weight)
T <- T %>% select(c("Generation", "Genotypes", "TW")) 


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


for (i in unique(T$Generation)){
 tmp <-  ggplot(subset(T, Generation == i), aes(Genotypes, TW)) + geom_bar(stat = 'identity') + 
   theme(axis.text.x = element_text(angle = 60)) + labs(y = "Total Weight (g)") +
   ggtitle(paste0("Generation ", i))
 ggsave(tmp, file = paste0("Bar_Plot_Total_Weight_Generation_", i, ".png"))
}


### Histograms for Average Fecundity Over Generations Single

ggplot(Rep_Single, aes(x = Avg_Fecundity)) +
  geom_histogram(binwidth = .9) +
  facet_grid(~Generation) +
  stat_bin(bins = 70) +
  labs(x = "Average Fecundity",
       y = "Frequency",
       title = "Average Fecundity Over Generations")


### Testing for Homogenity of Variance and ANOVA for Average Fecundity

ANOVA_Avg_Fec <- aov(Avg_Fecundity ~ as.factor(Generation), Rep_Single)
summary(ANOVA_Avg_Fec)
TukeyHSD(ANOVA_Avg_Fec)
leveneTest(Avg_Fecundity ~ as.factor(Generation), Rep_Single)

### QQ plot Average Fec Single

T <- Rep_Single %>% select(Avg_Fecundity, Generation, Genotypes) %>%  arrange(Generation)
T <- T %>% mutate(fec  = Avg_Fecundity)
T <- T %>% select(c("Generation", "Genotypes", "fec"))

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

### Histograms for Average Fit Over Generations Single

ggplot(Rep_Single, aes(x = Avg_Fit)) +
  geom_histogram(binwidth = .9) +
  facet_grid(~Generation) +
  stat_bin(bins = 70) +
  labs(x = "Average Fitness",
       y = "Frequency",
       title = "Average Fitness Over Generations")


### Testing for Homogenity of Variance and ANOVA for Average Fitness

ANOVA_Avg_Fit <- aov(Avg_Fit ~ as.factor(Generation), Rep_Single)
summary(ANOVA_Avg_Fit)
TukeyHSD(ANOVA_Avg_Fit)
leveneTest(Avg_Fit ~ as.factor(Generation), Rep_Single)

### QQ plot Average Fec Single

T <- Rep_Single %>% select(Avg_Fit, Generation, Genotypes) %>%  arrange(Generation)
T <- T %>% mutate(fit  = Avg_Fit)
T <- T %>% select(c("Generation", "Genotypes", "fit"))

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


ggplot(Rep_Single, aes(Genotypes, Avg_Total_Weight)) +
  geom_bar(stat = 'identity') + 
  facet_grid(~Generation)

### Bar Graphs Single For Each Generation (Total Weight)

for (i in unique)

ggplot(Rep_Single)
