library(tidyverse)
library(readr)
library(dplyr)
library(ggpubr)
library(ggplot2)
library(googlesheets4)
library(tidyr)

### Load Data

Full_Data <- read_delim("~/Documents/GitHub/Ag-Competition/Full_Data")
Side_By_Side_Replicates <- read_delim("~/Documents/GitHub/Ag-Competition/Side_by_Side_Replicates")

### Correlation of Total Seed Weight among replicates

ggplot(Side_By_Side_Replicates, aes(x=`Brown Bag Weight`, y=`Brown Bag Weight__2`, add = "reg.line")) +
  geom_point(alpha = .3) +
  geom_abline(slope =1, intercept= 0, color = "darkgreen") +
  stat_cor(label.y = 235, label.x = 20) +
  geom_smooth(method = "lm") +
  labs(x = "Rep 1 Total Weight (grams)",
       y= "Rep 2 Total Weight (grams)",
       title = "Total Seed Weight Correlation Between Rep 1 and Rep 2") +
  theme(
    axis.title.x = element_text(margin=margin(t=10)),
    panel.background = element_rect(fill=NA),
    legend.title = element_text(size=10),
    legend.key = element_blank()) +
  annotate("text",
           x = 75, 
           y = 250,
           label = paste("Number of OBS: ", nrow(Side_By_Side_Replicates)))

ggsave("Correlation_of_Rep1_and_Rep2_Total_SW.png")

res <- Side_By_Side_Replicates %>% summarise(residuals = resid(lm(`Brown Bag Weight__2` ~ `Brown Bag Weight`)))
Side_By_Side_Replicates$Residuals <- res$residuals

### Correlation of Fecundity by Replicate

ggplot(Side_By_Side_Replicates, aes(x=Fecundity, y= Fecundity__2)) +
  geom_point(alpha = .5) +
  stat_cor(label.y = 3400, label.x = 20) +
  geom_abline(slope =1, intercept= 0, color = "darkgreen") +
  geom_smooth(method = "lm") +
  labs(x = "Fecundity (rep 1)",
       y= "Fecundity (rep 2)",
       title = "Rep 1 vs. Rep 2 Fecundity Correlation") 

ggsave("Correlation_of_Rep1_and_Rep2_Fecundity.png")

res <- Side_By_Side_Replicates %>% summarise(residuals = resid(lm(Fecundity__2 ~ Fecundity)))
Side_By_Side_Replicates$Residuals <- res$residuals

### Correlation of 100 Seed Weight by Replicates (Single)

r1_100SW <- subset(df, replicate == "rep 1" & Condition == "single")
r2_100SW <- subset(df, replicate == "rep 2" & Condition == "single")
colnames(r2_100SW)[1:13] <- paste(colnames(r2_100SW)[c(1:13)], '_2', sep = '_')
SW100 <- inner_join(r1_100SW, r2_100SW, by = c("Genotypes" = "Genotypes__2","Condition" = "Condition__2"))

ggplot(SW100, aes(x=`Brown Bag Weight`, y=`Brown Bag Weight__2`)) +
  geom_jitter(alpha = .5) +
  stat_cor(label.y = 235, label.x = 20) +
  geom_abline(slope =1, intercept= 0, color = "darkgreen") +
  geom_smooth(method = "lm") +
  labs(x = "100 Seed Weight (rep 1)",
       y= "100 Seed Weight (rep 2)",
       title = "Single 100 Seed Weight Correlation by Replicate")

ggsave("Correlation_100_SW_Single.png")

res <- SW100 %>% summarise(residuals = resid(lm(`Brown Bag Weight__2` ~ `Brown Bag Weight`)))
SW100$Residuals <- res$residuals

### Correlation of 100 Seed Weight by Replicates (Mixed)

r1_100SW_mixed<- subset(Full_Data, replicate == "rep 1" & Condition == "mixed")
r2_100SW_mixed <- subset(Full_Data, replicate == "rep 2" & Condition == "mixed")
colnames(r2_100SW_mixed)[1:13] <- paste(colnames(r2_100SW_mixed)[c(1:13)], '_2', sep = '_')
SW100 <- inner_join(r1_100SW_mixed, r2_100SW_mixed, by = c("Genotypes" = "Genotypes__2","Condition" = "Condition__2"))

ggplot(SW100, aes(x=`Brown Bag Weight`, y=`Brown Bag Weight__2`)) +
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

### Histograms For Total Weight by Generation
ggplot(Full_Data, aes(x= `Brown Bag Weight`)) +
  geom_histogram() +
  labs(x= "Total Weight (grams)",
       y= "Frequency",
       title = "Total Weight over Generations") +
  facet_wrap(~Generation)

ANOVA_Total_weight_over_Generations <- aov(`Brown Bag Weight` ~ as.factor(Generation), Full_Data)
summary(ANOVA_Total_weight_over_Generations)
TukeyHSD(ANOVA_Total_weight_over_Generations)

ggsave("Histograms_Total_Weight_Over_Generations.png")

### Overlapping Histograms For Total Weight by Generation
Full_Data$Generation <- as.factor(Full_Data$Generation)
ggplot(Full_Data, aes(x = `Brown Bag Weight`, group = Generation, fill = Generation, group = Generation)) +
  geom_histogram(alpha = .5, position = 'identity', binwidth = 10) +
  scale_fill_brewer(palette = "Blues")

ggsave("Overlapping_Histograms_Total_Weight_Over_Generations.png")

### Histograms For 100 Seed Weight by Generation
ggplot(Full_Data, aes(x = `100 seed weight`)) +
  geom_histogram() +
  labs(x = "100 Seed Weight (grams)",
       y = "frequency",
       title = "100 Seed Weight by Generations") +
  facet_wrap(~Generation)

ANOVA_100_SW <- aov(`100 seed weight` ~ as.factor(Generation), Full_Data)
summary(ANOVA_100_SW)
TukeyHSD(ANOVA_100_SW)

ggsave("Histograms_SW_100_Over_Generations.png")

### Overlapping Histogram for 100 Seed Weight by Generation
Full_Data$Generation <- as.factor(Full_Data$Generation)
ggplot(Full_Data, aes(x = `100 seed weight`, group = Generation, fill = Generation)) +
  geom_histogram(alpha = .5, position = 'identity') +
  scale_fill_brewer(palette = "Greens") +
  labs(x = "100 Seed Weight (grams)",
       y = "frequency",
       title = "100 Seed Weight by Generations")

ggsave("Overplapping_Histograms_SW_100_Over_Generations.png")

