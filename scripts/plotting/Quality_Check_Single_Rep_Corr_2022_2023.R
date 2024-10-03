library(tidyverse)
library(readr)
library(dplyr)
library(ggpubr)
library(ggplot2)

library(tidyr)
library(car)
library(gridExtra)
library(dunn.test)

Replicate_corr_tbl_Single <- Replicate_corr_tbl %>% filter(Condition == "single")

### Correlation of Replicates (Total Seed Weight) 

g2 <- ggplot(Replicate_corr_tbl_Single, aes(x=`Brown Bag Weight`, y=`Brown Bag Weight_2`, add = "reg.line")) +
  geom_point(alpha = .5) +
  geom_abline(slope =1, intercept= 0, color = "red") +
  stat_cor(label.y = 235, label.x = 20) +
  geom_smooth(method = "lm") +
  labs(x = "Rep 1",
       y= "Rep 2",
       title = "Total Weight 2022-2023") 
ggsave("scripts/plotting/Correlation_of_Single_Rep_2022_2023_TW.png")

res <- Replicate_corr_tbl %>% summarise(residuals = resid(lm(`Brown Bag Weight_2` ~ `Brown Bag Weight`)))
Replicate_corr_tbl$Residuals <- res$residuals

### Correlation of Replicates (Flowering Time)

g4 <- ggplot(Replicate_corr_tbl_Single, aes(FT, FT_2)) +
  geom_jitter(alpha = .5) +
  geom_smooth(method = 'lm') +
  stat_cor(label.x = 100, label.y = 135) +
  geom_abline(slope = 1, intercept = 0, col = "red") +
  labs(x = "Rep 1",
       y = "Rep 2",
       title = "Flowering Time 2022-2023")
ggsave("scripts/plotting/Correlation_of_Single_Rep_2022_2023_FT.png")

### Correlation of Replicates (Fecundity)

g6 <- ggplot(Replicate_corr_tbl_Single, aes(x=Fecundity, y= Fecundity_2)) +
  geom_point(alpha = .5) +
  stat_cor(label.y = 4200, label.x = 20) +
  geom_abline(slope =1, intercept= 0, color = "red") +
  geom_smooth(method = "lm") +
  labs(x = "Rep 1",
       y= "Rep 2",
       title = "Fecundity 2022-2023") 
ggsave("scripts/plotting/Correlation_of_Single_Rep_2022_2023_Fecundity.png")

res <- Replicate_corr_tbl %>% summarise(residuals = resid(lm(Fecundity_2 ~ Fecundity)))
Replicate_corr_tbl$Residuals <- res$residuals

### Correlation of Replicates (Fitness)

g8 <- ggplot(Replicate_corr_tbl_Single, aes(Fitness, Fitness_2)) +
  geom_jitter(alpha = .5) +
  geom_smooth(method = 'lm') +
  stat_cor(label.x = 100, label.y = 40000) +
  geom_abline(slope = 1, intercept = 0, col = "red") +
  labs(x = "Rep 1",
       y = "Rep 2",
       title = "Fitness 2022-2023")
ggsave("scripts/plotting/Correlation_of_Single_Rep_2022_2023_Fitness.png")

res <- Replicate_corr_tbl %>% summarise(residuals = resid(lm(Fitness_2 ~ Fitness)))
Replicate_corr_tbl$Residuals <- res$residuals


### Correlation of 100 Seed Weight between Replicates (Single)

r1_100SW <- subset(Full_Data, Replicate == "rep 1" & Condition == "single")
r2_100SW <- subset(Full_Data, Replicate == "rep 2" & Condition == "single")
colnames(r2_100SW)[1:13] <- paste(colnames(r2_100SW)[c(1:13)], '_2', sep = '_')
SW100 <- inner_join(r1_100SW, r2_100SW, by = c("Genotype" = "Genotype__2","Condition" = "Condition__2"))

ggplot(SW100, aes(x=`Brown Bag Weight`, y=`Brown Bag Weight__2`)) +
  geom_jitter(alpha = .5) +
  stat_cor(label.y = 235, label.x = 20) +
  geom_abline(slope =1, intercept= 0, color = "darkgreen") +
  geom_smooth(method = "lm") +
  labs(x = "100 Seed Weight (rep 1)",
       y= "100 Seed Weight (rep 2)",
       title = "Correlation of Single Replicate for 100 Seed Weight")

ggsave("scripts/plotting/Correlation_100_SW_Single.png")

res <- SW100 %>% summarise(residuals = resid(lm(`Brown Bag Weight__2` ~ `Brown Bag Weight`)))
SW100$Residuals <- res$residuals

### Correlation of 100 Seed Weight by Replicates (Mixed)

r1_100SW_mixed<- subset(Full_Data, Replicate == "rep 1" & Condition == "mixed")
r2_100SW_mixed <- subset(Full_Data, Replicate == "rep 2" & Condition == "mixed")
colnames(r2_100SW_mixed)[1:13] <- paste(colnames(r2_100SW_mixed)[c(1:13)], '_2', sep = '_')
SW100 <- inner_join(r1_100SW_mixed, r2_100SW_mixed, by = c("Genotype" = "Genotype__2","Condition" = "Condition__2"))

y <- ggplot(SW100, aes(x=`Brown Bag Weight`, y=`Brown Bag Weight__2`)) +
  geom_jitter(alpha = .5) +
  stat_cor(label.y = 235, label.x = 20) +
  geom_abline(slope =1, intercept= 0, color = "darkgreen") +
  geom_smooth(method = "lm") +
  labs(x = "100 Seed Weight (rep 1)",
       y= "100 Seed Weight (rep 2)",
       title = "Mixed 100 Seed Weight Correlation by Replicate")

ggsave("scripts/plotting/Correlation_100_SW_Mixed.png")

res <- SW100 %>% summarise(residuals = resid(lm(`Brown Bag Weight__2` ~ `Brown Bag Weight`)))
SW100$Residuals <- res$residuals

l <- ggplot(SW100, aes(`100 seed weight`, Residuals)) +
  geom_point(alpha = .5) +
  geom_hline(yintercept = 0, col = "red") +
  labs(title = "Residuals for 100 Seed Weight (Mixed)")

grid.arrange(y, l, ncol = 2, nrow = 1) +
ggsave("scripts/plotting/Correlation_100_SW_By_Replicates_Mixed.png")


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

### Q-Q Plot For Total Weight ~ Generation

T <- Full_Data %>% select(`Brown Bag Weight`, Generation, Genotype) %>%  arrange(Generation)
T <- T %>% mutate(Total_Weight = `Brown Bag Weight`)
T <- T %>% select(c("Generation", "Genotype", "Total_Weight"))

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

### Overlapping Histograms For Total Weight by Generation

Average_Haplo_rep$Generation <- as.factor(Average_Haplo_rep$Generation)
ggplot(Average_Haplo_rep, aes(x = `Brown Bag Weight`, group = Generation, fill = Generation)) +
  geom_histogram(alpha = .5, binwidth = .5) +
  stat_bin(bins = 50) +
  scale_fill_brewer(palette = "Blues")

ggsave("scripts/plotting/Overlapping_Histograms_Total_Weight_Over_Generations.png")

### Histograms For 100 Seed Weight by Generation

ggplot(Full_Data, aes(x = `100 seed weight`)) +
  geom_histogram(binwidth = .1) +
  stat_bin(bins = 50) +
  labs(x = "100 Seed Weight (grams)",
       y = "frequency",
       title = "Generational Change in Average 100 Seed Weight") +
  facet_grid(~Generation)

### Testing for Homogenity of Variance and ANOVA for 100 SW

ANOVA_100_SW <- aov(`100 seed weight` ~ as.factor(Generation), Full_Data)
summary(ANOVA_100_SW)
TukeyHSD(ANOVA_100_SW)
leveneTest(`100 seed weight` ~ as.factor(Generation), Full_Data)

ggsave("scripts/plotting/Histograms_Average_SW_100_Over_Generations.png")

### Overlapping Histogram for 100 Seed Weight by Generation

Full_Data$Generation <- as.factor(Full_Data$Generation)
ggplot(Full_Data, aes(x = `100 seed weight`, group = Generation, fill = Generation)) +
  geom_histogram(alpha = .5, position = 'identity') +
  scale_fill_brewer(palette = "Greens") +
  labs(x = "100 Seed Weight (grams)",
       y = "frequency",
       title = "100 Seed Weight by Generations")

ggsave("scripts/plotting/Overplapping_Histograms_SW_100_Over_Generations.png")
