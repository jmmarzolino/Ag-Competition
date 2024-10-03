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
library(gridExtra)
library(car)
library(ggrepel)
library(ggExtra)

setwd("/bigdata/koeniglab/jmarz001/Ag-Competition/data/")
PHENO_FULL <- read_delim("FT_FITNESS.tsv")

# MIXED
# Creating a function to easily graph all phenos
# x = dataframe
# y = phenotype graphed in quotes (included in title)
# z = Season (included in title)

Graphing_Corr <- function(x, y, z){
  ggplot(x, aes(`1`, `2`), add = "reg.line") +
    geom_jitter(alpha = .5) +
    geom_abline(slope = 1, intercept = 0, color = 'red') +
    geom_smooth(method = 'lm') +
    labs(x = "Rep 1",
         y = "Rep 2",
         title = paste("Correlation of", sep = " ", y, "Replicates", z))
}

# Creating a dataframe that filters by mixed condition in the 2021-2022 season

PHENO_MIXED_2022 <- PHENO_FULL %>% filter(Condition == "mixed" & Exp_year == 2022)

# Mixed Correlation 2021-2022 TW

cmp <- PHENO_MIXED_2022 %>% group_by(Genotype, Generation, Replicate) %>% summarise(TOTAL_MASS) %>% spread(key = Replicate, value = TOTAL_MASS) %>% ungroup()
p1 <- Graphing_Corr(cmp, "Total Weight (g)", "2021-2022") +
  stat_cor(label.y = 170) +
  xlim(0, 250) +
  ylim(0, 250)

outlier_upper <- quantile(PHENO_MIXED_2022$TOTAL_MASS, .75, na.rm = T) + (1.5 * IQR(PHENO_MIXED_2022$TOTAL_MASS, na.rm = T))
outlier_lower <- quantile(PHENO_MIXED_2022$TOTAL_MASS, .25, na.rm = T) - (1.5 * IQR(PHENO_MIXED_2022$TOTAL_MASS, na.rm = T))
outlier_data <- PHENO_MIXED_2022 %>% filter(TOTAL_MASS > outlier_upper | TOTAL_MASS < outlier_lower)

# Residual Plotting
cmp <- na.omit(cmp)
residual_df <- cmp %>% reframe(Residuals = resid(lm(`2` ~ `1`)))
cmp$Residuals <- residual_df$Residuals
Standard_dev <- sd(cmp$Residuals)
cmp <- cmp %>% mutate(Residuals = Residuals/Standard_dev)
tmp_outlier <- cmp %>% filter(Residuals > 2 | Residuals < -2)

a <- ggplot(cmp, aes(x = `1`, y = Residuals)) +
  geom_jitter() +
  geom_hline(yintercept = 0, color = 'red') +
  geom_hline(yintercept = 2, color = "blue") +
  geom_hline(yintercept = -2, color = 'blue') +
  scale_y_continuous(breaks = seq(-4, 4, 1)) +
  geom_text_repel(label = ifelse(cmp$Residuals > 2 | cmp$Residuals < -2,
                                 cmp$Genotype,
                                 ""), size = 3, hjust = 1, max.overlaps = 30) +
  labs(x = 'Total Weight (grams)',
       title = "Residual Plot Mixed Total Weight 2022")


# Mixed Correlation 2021-2022 Fecundity

cmp <- PHENO_MIXED_2022 %>% group_by(Genotype, Generation, Replicate) %>% summarise(FECUNDITY) %>% spread(key = Replicate, value = FECUNDITY) %>% ungroup()
p2 <- Graphing_Corr(cmp, "Relative Fecundity", "2021-2022") +
  stat_cor(label.y = 2.5, label.x = 6.5) +
  ylim(-1.1, 3)+
  xlim(-1.1, 3) +
  geom_text_repel(label = ifelse(cmp$`1` > 3 | cmp$`2` > 3,
                                 cmp$Genotype,
                                 ""), size = 3, hjust =1, max.overlaps = 30)

# Residual Plotting
cmp <- na.omit(cmp)
residual_df <- cmp %>% reframe(Residuals = resid(lm(`2` ~ `1`)))
cmp$Residuals <- residual_df$Residuals
Standard_dev <- sd(cmp$Residuals)
cmp <- cmp %>% mutate(Residuals = Residuals/Standard_dev)
tmp_outlier <- cmp %>% filter(Residuals > 2 | Residuals < -2)

b <- ggplot(cmp, aes(x = `1`, y = Residuals)) +
  geom_jitter() +
  geom_hline(yintercept = 0, color = 'red') +
  geom_hline(yintercept = 2, color = "blue") +
  geom_hline(yintercept = -2, color = 'blue') +
  scale_y_continuous(breaks = seq(-4, 4, 1)) +
  geom_text_repel(label = ifelse(cmp$Residuals > 2 | cmp$Residuals < -2,
                                 cmp$Genotype,
                                 ""), size = 3, hjust = 1, max.overlaps = 30) +
  labs(x = 'Centered Fecundity',
       title = "Residual Plot Mixed Centered Fecundity 2022")

# Mixed Correlation 2021-2022 Fitness

cmp <- PHENO_MIXED_2022 %>% group_by(Genotype, Generation, Replicate) %>% summarise(ABS_FITNESS) %>% spread(key = Replicate, value = ABS_FITNESS) %>% ungroup()
p3 <- Graphing_Corr(cmp, "Absolute Fitness", "2021-2022") +
  stat_cor(label.y = 3) +
  xlim(-2.5, 4.2) +
  ylim(-2.5, 4.2)

# Residual Plotting
cmp <- na.omit(cmp)
residual_df <- cmp %>% reframe(Residuals = resid(lm(`2` ~ `1`)))
cmp$Residuals <- residual_df$Residuals
Standard_dev <- sd(cmp$Residuals)
cmp <- cmp %>% mutate(Residuals = Residuals/Standard_dev)
tmp_outlier <- cmp %>% filter(Residuals > 2 | Residuals < -2)

c <- ggplot(cmp, aes(x = `1`, y = Residuals)) +
  geom_jitter() +
  geom_hline(yintercept = 0, color = 'red') +
  geom_hline(yintercept = 2, color = "blue") +
  geom_hline(yintercept = -2, color = 'blue') +
  scale_y_continuous(breaks = seq(-4, 4, 1)) +
  geom_text_repel(label = ifelse(cmp$Residuals > 2 | cmp$Residuals < -2,
                                 cmp$Genotype,
                                 ""), size = 3, hjust = 1, max.overlaps = 30) +
  labs(x = 'Centered Fitness',
       title = "Residual Plot Mixed Centered Fitness 2022")

# Mixed Correlation 2021-2022 Flowering Time

cmp <- PHENO_MIXED_2022 %>% group_by(Genotype, Generation, Replicate) %>% summarise(FT) %>% spread(key = Replicate, value = FT) %>% ungroup()
p4 <- Graphing_Corr(cmp, "Flowering Time", "2021-2022") +
  stat_cor(label.y = 105) +
  xlim(75, 135) +
  ylim(75, 135)

outlier_upper <- quantile(PHENO_MIXED_2022$FT, .75, na.rm = T) + (1.5 * IQR(PHENO_MIXED_2022$FT, na.rm = T))
outlier_lower <- quantile(PHENO_MIXED_2022$FT, .25, na.rm = T) - (1.5 * IQR(PHENO_MIXED_2022$FT, na.rm = T))
outlier_data <- PHENO_MIXED_2022 %>% filter(FT > outlier_upper | FT < outlier_lower)

# Residual Plotting
cmp <- na.omit(cmp)
residual_df <- cmp %>% reframe(Residuals = resid(lm(`2` ~ `1`)))
cmp$Residuals <- residual_df$Residuals
Standard_dev <- sd(cmp$Residuals)
cmp <- cmp %>% mutate(Residuals = Residuals/Standard_dev)
tmp_outlier <- cmp %>% filter(Residuals > 2 | Residuals < -2)

d <- ggplot(cmp, aes(x = `1`, y = Residuals)) +
  geom_jitter() +
  geom_hline(yintercept = 0, color = 'red') +
  geom_hline(yintercept = 2, color = "blue") +
  geom_hline(yintercept = -2, color = 'blue') +
  scale_y_continuous(breaks = seq(-4, 4, 1)) +
  geom_text_repel(label = ifelse(cmp$Residuals > 2 | cmp$Residuals < -2,
                                 cmp$Genotype,
                                 ""), size = 3, hjust = 1, max.overlaps = 30) +
  labs(x = 'Flowering Time',
       title = "Residual Plot Mixed Flowering Time 2022")


# Mixed Correlation 2021-2022 100 SW

cmp <- PHENO_MIXED_2022 %>% group_by(Genotype, Generation, Replicate) %>% summarise(SEED_WEIGHT_100) %>% spread(key = Replicate, value = SEED_WEIGHT_100) %>% ungroup()
p5 <- Graphing_Corr(cmp, "100 Seed Weight (g)", "2021-2022") +
  stat_cor(label.y = 4) +
  xlim(2.4, 6.1) +
  ylim(2.4, 6.1)


outlier_upper <- quantile(PHENO_MIXED_2022$SEED_WEIGHT_100, .75, na.rm = T) + (1.5 * IQR(PHENO_MIXED_2022$SEED_WEIGHT_100, na.rm = T))
outlier_lower <- quantile(PHENO_MIXED_2022$SEED_WEIGHT_100, .25, na.rm = T) - (1.5 * IQR(PHENO_MIXED_2022$SEED_WEIGHT_100, na.rm = T))
outlier_data <- PHENO_MIXED_2022 %>% filter(SEED_WEIGHT_100 > outlier_upper | SEED_WEIGHT_100 < outlier_lower)

# Residual Plotting
cmp <- na.omit(cmp)
residual_df <- cmp %>% reframe(Residuals = resid(lm(`2` ~ `1`)))
cmp$Residuals <- residual_df$Residuals
Standard_dev <- sd(cmp$Residuals)
cmp <- cmp %>% mutate(Residuals = Residuals/Standard_dev)
tmp_outlier <- cmp %>% filter(Residuals > 2 | Residuals < -2)

e <- ggplot(cmp, aes(x = `1`, y = Residuals)) +
  geom_jitter() +
  geom_hline(yintercept = 0, color = 'red') +
  geom_hline(yintercept = 2, color = "blue") +
  geom_hline(yintercept = -2, color = 'blue') +
  scale_y_continuous(breaks = seq(-4, 4, 1)) +
  geom_text_repel(label = ifelse(cmp$Residuals > 2 | cmp$Residuals < -2,
                                 cmp$Genotype,
                                 ""), size = 3, hjust = 1, max.overlaps = 35) +
  labs(x = '100 Seed Weight (grams)',
       title = "Residual Plot Mixed 100 Seed Weight (grams) 2022")

# Creating a dataframe that filters by mixed condition in the 2022-2023 season

PHENO_MIXED_2023 <- PHENO_FULL %>% filter(Condition == "mixed" & Exp_year == 2023)

# Mixed Correlation 2022-2023 TW

cmp <- PHENO_MIXED_2023 %>% group_by(Genotype, Generation, Replicate) %>% summarise(TOTAL_MASS) %>% spread(key = Replicate, value = TOTAL_MASS) %>% ungroup()
p6 <- Graphing_Corr(cmp, "Total Weight (g)", "2022-2023") +
  stat_cor(label.y = 200) +
  xlim(0, 250) +
  ylim(0, 250)

PHENO_tmp <- PHENO2023 %>% filter(Condition == "mixed")

outlier_upper <- quantile(PHENO_tmp$TOTAL_MASS, .75, na.rm = T) + (1.5 * IQR(PHENO_tmp$TOTAL_MASS, na.rm = T))
outlier_lower <- quantile(PHENO_tmp$TOTAL_MASS, .25, na.rm = T) - (1.5 * IQR(PHENO_tmp$TOTAL_MASS, na.rm = T))
outlier_data <- PHENO_tmp %>% filter(TOTAL_MASS > outlier_upper | TOTAL_MASS < outlier_lower)

# Residual Plotting
cmp <- na.omit(cmp)
residual_df <- cmp %>% reframe(Residuals = resid(lm(`2` ~ `1`)))
cmp$Residuals <- residual_df$Residuals
Standard_dev <- sd(cmp$Residuals)
cmp <- cmp %>% mutate(Residuals = Residuals/Standard_dev)
tmp_outlier <- cmp %>% filter(Residuals > 2 | Residuals < -2)

f <- ggplot(cmp, aes(x = `1`, y = Residuals)) +
  geom_jitter() +
  geom_hline(yintercept = 0, color = 'red') +
  geom_hline(yintercept = 2, color = "blue") +
  geom_hline(yintercept = -2, color = 'blue') +
  scale_y_continuous(breaks = seq(-4, 4, 1)) +
  geom_text_repel(label = ifelse(cmp$Residuals > 2 | cmp$Residuals < -2,
                                 cmp$Genotype,
                                 ""), size = 3, hjust = 1, max.overlaps = 35) +
  labs(x = 'Total Weight (grams)',
       title = "Residual Plot Mixed Total Weight (grams) 2023")


# Mixed Correlation 2022-2023 Fecundity

cmp <- PHENO_MIXED_2023 %>% group_by(Genotype, Generation, Replicate) %>% summarise(FECUNDITY) %>% spread(key = Replicate, value = FECUNDITY) %>% ungroup()
p7 <- Graphing_Corr(cmp, "Relative Fecundity", "2022-2023") +
  stat_cor(label.y = 13) +
  xlim(-1.1, 3) +
  ylim(-1.1, 3) +
  geom_text_repel(label = ifelse(cmp$`1` > 3 | cmp$`2` > 3,
                                 cmp$Genotype,
                                 ""), size = 3, hjust =1, max.overlaps = 30)



# Residual Plotting
cmp <- na.omit(cmp)
residual_df <- cmp %>% reframe(Residuals = resid(lm(`2` ~ `1`)))
cmp$Residuals <- residual_df$Residuals
Standard_dev <- sd(cmp$Residuals)
cmp <- cmp %>% mutate(Residuals = Residuals/Standard_dev)
tmp_outlier <- cmp %>% filter(Residuals > 2 | Residuals < -2)

g <- ggplot(cmp, aes(x = `1`, y = Residuals)) +
  geom_jitter() +
  geom_hline(yintercept = 0, color = 'red') +
  geom_hline(yintercept = 2, color = "blue") +
  geom_hline(yintercept = -2, color = 'blue') +
  scale_y_continuous(breaks = seq(-4, 4, 1)) +
  geom_text_repel(label = ifelse(cmp$Residuals > 2 | cmp$Residuals < -2,
                                 cmp$Genotype,
                                 ""), size = 3, hjust = 1, max.overlaps = 35) +
  labs(x = 'Centered Fecundity',
       title = "Residual Plot Mixed Centered Fecundity 2023")

# Mixed Correlation 2022-2023 Fitness

cmp <- PHENO_MIXED_2023 %>% group_by(Genotype, Generation, Replicate) %>% summarise(ABS_FITNESS) %>% spread(key = Replicate, value = ABS_FITNESS) %>% ungroup()
p8 <- Graphing_Corr(cmp, "Absolute Fitness", "2022-2023") +
  stat_cor(label.y = 3) +
  xlim(-2.5, 4.2) +
  ylim(-2.5, 4.2)


# Residual Plotting
cmp <- na.omit(cmp)
residual_df <- cmp %>% reframe(Residuals = resid(lm(`2` ~ `1`)))
cmp$Residuals <- residual_df$Residuals
Standard_dev <- sd(cmp$Residuals)
cmp <- cmp %>% mutate(Residuals = Residuals/Standard_dev)
tmp_outlier <- cmp %>% filter(Residuals > 2 | Residuals < -2)

h <- ggplot(cmp, aes(x = `1`, y = Residuals)) +
  geom_jitter() +
  geom_hline(yintercept = 0, color = 'red') +
  geom_hline(yintercept = 2, color = "blue") +
  geom_hline(yintercept = -2, color = 'blue') +
  scale_y_continuous(breaks = seq(-4, 4, 1)) +
  geom_text_repel(label = ifelse(cmp$Residuals > 2 | cmp$Residuals < -2,
                                 cmp$Genotype,
                                 ""), size = 3, hjust = 1, max.overlaps = 35) +
  labs(x = 'Centered Absolute Fitness',
       title = "Residual Plot Mixed Centered Absolute Fitness 2023")

# Mixed Correlation 2022-2023 Flowering Time

cmp <- PHENO_MIXED_2023 %>% group_by(Genotype, Generation, Replicate) %>% summarise(FT) %>% spread(key = Replicate, value = FT) %>% ungroup()
p9 <- Graphing_Corr(cmp, "Flowering Time", "2022-2023") +
  stat_cor(label.y = 125) +
  xlim(75, 135) +
  ylim(75, 135)


PHENO_tmp <- PHENO2023 %>% filter(Condition == "mixed")
outlier_upper <- quantile(PHENO_tmp$FT, .75, na.rm = T) + (1.5 * IQR(PHENO_tmp$FT, na.rm = T))
outlier_lower <- quantile(PHENO_tmp$FT, .25, na.rm = T) - (1.5 * IQR(PHENO_tmp$FT, na.rm = T))
outlier_data <- PHENO_tmp %>% filter(FT > outlier_upper | FT < outlier_lower)


# Residual Plotting
cmp <- na.omit(cmp)
residual_df <- cmp %>% reframe(Residuals = resid(lm(`2` ~ `1`)))
cmp$Residuals <- residual_df$Residuals
Standard_dev <- sd(cmp$Residuals)
cmp <- cmp %>% mutate(Residuals = Residuals/Standard_dev)
tmp_outlier <- cmp %>% filter(Residuals > 2 | Residuals < -2)

i <- ggplot(cmp, aes(x = `1`, y = Residuals)) +
  geom_jitter() +
  geom_hline(yintercept = 0, color = 'red') +
  geom_hline(yintercept = 2, color = "blue") +
  geom_hline(yintercept = -2, color = 'blue') +
  scale_y_continuous(breaks = seq(-4, 4, 1)) +
  geom_text_repel(label = ifelse(cmp$Residuals > 2 | cmp$Residuals < -2,
                                 cmp$Genotype,
                                 ""), size = 3, hjust = 1, max.overlaps = 35) +
  labs(x = 'Flowering Time',
       title = "Residual Plot Mixed Flowering Time 2023")

# Mixed Correlation 2022-2023 100 SW

cmp <- PHENO_MIXED_2023 %>% group_by(Genotype, Generation, Replicate) %>% summarise(SEED_WEIGHT_100) %>% spread(key = Replicate, value = SEED_WEIGHT_100) %>% ungroup()
p10 <- Graphing_Corr(cmp, "100 Seed Weight (g)", "2022-2023") +
  stat_cor(label.y = 6) +
  xlim(2.4, 6.1) +
  ylim(2.4, 6.1)

PHENO_tmp <- PHENO2023 %>% filter(Condition == "mixed")
outlier_upper <- quantile(PHENO_tmp$SEED_WEIGHT_100, .75, na.rm = T) + (1.5 * IQR(PHENO_tmp$SEED_WEIGHT_100, na.rm = T))
outlier_lower <- quantile(PHENO_tmp$SEED_WEIGHT_100, .25, na.rm = T) - (1.5 * IQR(PHENO_tmp$SEED_WEIGHT_100, na.rm = T))
outlier_data <- PHENO_tmp %>% filter(SEED_WEIGHT_100 > outlier_upper | SEED_WEIGHT_100 < outlier_lower)

# Residual Plotting
cmp <- na.omit(cmp)
residual_df <- cmp %>% reframe(Residuals = resid(lm(`2` ~ `1`)))
cmp$Residuals <- residual_df$Residuals
Standard_dev <- sd(cmp$Residuals)
cmp <- cmp %>% mutate(Residuals = Residuals/Standard_dev)
tmp_outlier <- cmp %>% filter(Residuals > 2 | Residuals < -2)

j <- ggplot(cmp, aes(x = `1`, y = Residuals)) +
  geom_jitter() +
  geom_hline(yintercept = 0, color = 'red') +
  geom_hline(yintercept = 2, color = "blue") +
  geom_hline(yintercept = -2, color = 'blue') +
  scale_y_continuous(breaks = seq(-4, 4, 1)) +
  geom_text_repel(label = ifelse(cmp$Residuals > 2 | cmp$Residuals < -2,
                                 cmp$Genotype,
                                 ""), size = 3, hjust = 1, max.overlaps = 35) +
  labs(x = '100 Seed Weight (grams)',
       title = "Residual Plot Mixed 100 Seed Weight (grams) 2023")

y <- arrangeGrob(a,b,c,d,e,f,g,h,i,j, nrow = 2, ncol = 5)
ggsave("/bigdata/koeniglab/jmarz001/Ag-Competition/results/Residuals_Mixed_Both_Years.png", y, width = 26, height =14)

z <- arrangeGrob(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, nrow = 2, ncol = 5, top = "Correlation of Mixed Replicates")
ggsave("/bigdata/koeniglab/jmarz001/Ag-Competition/results/Correlation_plots_Mixed.png", z, width = 32, height = 16)




##############################################################


# SINGLE

# Creating a function to easily graph all phenos
# x = dataframe
# y = phenotype graphed in quotes (included in title)
# z = Season (included in title)

Graphing_Corr <- function(x, y, z){
  ggplot(x, aes(`1`, `2`), add = "reg.line") +
    geom_jitter(alpha = .5) +
    geom_abline(slope = 1, intercept = 0, color = 'red') +
    geom_smooth(method = 'lm') +
    labs(x = "Rep 1",
         y = "Rep 2",
         title = paste("Correlation of", sep = " ", y, "Replicates", z))
}

# Load Data and filter for 2021-2022 Year and Single Condition

PHENO_SINGLE_2022 <- PHENO_FULL %>% filter(Condition == "single" & Exp_year == 2022)

# Single Correlation 2021-2022 TW

cmp <- PHENO_SINGLE_2022 %>% group_by(Genotype, Generation, Replicate) %>% summarise(TOTAL_MASS) %>% spread(key = Replicate, value = TOTAL_MASS) %>% ungroup()
p1 <- Graphing_Corr(cmp, "Total Weight (g)", "2021-2022") +
  stat_cor(label.y = 170) +
  ylim(0, 230) +
  xlim(0, 230)

outlier_upper <- quantile(PHENO_SINGLE_2022$TOTAL_MASS, .75, na.rm = T) + (1.5 * IQR(PHENO_SINGLE_2022$TOTAL_MASS, na.rm = T))
outlier_lower <- quantile(PHENO_SINGLE_2022$TOTAL_MASS, .25, na.rm = T) - (1.5 * IQR(PHENO_SINGLE_2022$TOTAL_MASS, na.rm = T))
outlier_data <- PHENO_SINGLE_2022 %>% filter(TOTAL_MASS > outlier_upper | TOTAL_MASS < outlier_lower)

ggplot(PHENO_SINGLE_2022, aes(x = TOTAL_MASS)) +
  geom_histogram(bins = 60) +
  geom_vline(xintercept = outlier_upper, color = "red")

# Residual Plotting
cmp <- na.omit(cmp)
residual_df <- cmp %>% reframe(Residuals = resid(lm(`2` ~ `1`)))
cmp$Residuals <- residual_df$Residuals
Standard_dev <- sd(cmp$Residuals)
cmp <- cmp %>% mutate(Residuals = Residuals/Standard_dev)
tmp_outlier <- cmp %>% filter(Residuals > 2 | Residuals < -2)

a <- ggplot(cmp, aes(x = `1`, y = Residuals)) +
  geom_jitter() +
  geom_hline(yintercept = 0, color = 'red') +
  geom_hline(yintercept = 2, color = "blue") +
  geom_hline(yintercept = -2, color = 'blue') +
  scale_y_continuous(breaks = seq(-4, 4, 1)) +
  geom_text_repel(label = ifelse(cmp$Residuals > 2 | cmp$Residuals < -2,
                                 cmp$Genotype,
                                 ""), size = 3, hjust = 1, max.overlaps = 30) +
  labs(x = 'Total Weight (grams)',
       title = "Residual Plot Single Total Weight 2022")

# Single Correlation 2021-2022 Fecundity

cmp <- PHENO_SINGLE_2022 %>% group_by(Genotype, Generation, Replicate) %>% summarise(FECUNDITY) %>% spread(key = Replicate, value = FECUNDITY) %>% ungroup()
p2 <- Graphing_Corr(cmp, "Relative Fecundity", "2021-2022") +
  stat_cor(label.y = 2, label.x = 1) +
  xlim(-1.1, 3) +
  ylim(-1.1, 3) +
  geom_text_repel(label = ifelse(cmp$`1` > 3 | cmp$`2` > 3,
                                 cmp$Genotype,
                                 ""), size = 3, hjust =1, max.overlaps = 30)

# Residual Plotting
cmp <- na.omit(cmp)
residual_df <- cmp %>% reframe(Residuals = resid(lm(`2` ~ `1`)))
cmp$Residuals <- residual_df$Residuals
Standard_dev <- sd(cmp$Residuals)
cmp <- cmp %>% mutate(Residuals = Residuals/Standard_dev)
tmp_outlier <- cmp %>% filter(Residuals > 2 | Residuals < -2)

b <- ggplot(cmp, aes(x = `1`, y = Residuals)) +
  geom_jitter() +
  geom_hline(yintercept = 0, color = 'red') +
  geom_hline(yintercept = 2, color = "blue") +
  geom_hline(yintercept = -2, color = 'blue') +
  scale_y_continuous(breaks = seq(-4, 4, 1)) +
  geom_text_repel(label = ifelse(cmp$Residuals > 2 | cmp$Residuals < -2,
                                 cmp$Genotype,
                                 ""), size = 3, hjust = .8, max.overlaps = 20) +
  labs(x = "Centered Fecundity",
       title = "Residuals Single Centered Fecundity 2022")

# Single Correlation 2021-2022 Fitness

cmp <- PHENO_SINGLE_2022 %>% group_by(Genotype, Generation, Replicate) %>% summarise(ABS_FITNESS) %>% spread(key = Replicate, value = ABS_FITNESS) %>% ungroup()
p3 <- Graphing_Corr(cmp, "Absolute Fitness", "2021-2022") +
  stat_cor(label.y = 3) +
  ylim(-1.8, 5) +
  xlim(-1.8, 5)

# Residual Plotting
cmp <- na.omit(cmp)
residual_df <- cmp %>% reframe(Residuals = resid(lm(`2` ~ `1`)))
cmp$Residuals <- residual_df$Residuals
Standard_dev <- sd(cmp$Residuals)
cmp <- cmp %>% mutate(Residuals = Residuals/Standard_dev)
tmp_outlier <- cmp %>% filter(Residuals > 2 | Residuals < -2)

c <- ggplot(cmp, aes(x = `1`, y = Residuals)) +
  geom_jitter() +
  geom_hline(yintercept = 0, color = 'red') +
  geom_hline(yintercept = 2, color = "blue") +
  geom_hline(yintercept = -2, color = 'blue') +
  scale_y_continuous(breaks = seq(-4, 4, 1)) +
  geom_text_repel(label = ifelse(cmp$Residuals > 2 | cmp$Residuals < -2,
                                 cmp$Genotype,
                                 ""), size = 3, hjust = .8, max.overlaps = 20) +
  labs(x = "Centered Absolute Fitness",
       title = "Residuals Single Centered Absolute Fitness 2022")

# Single Correlation 2021-2022 Flowering Time

cmp <- PHENO_SINGLE_2022 %>% group_by(Genotype, Generation, Replicate) %>% summarise(FT) %>% spread(key = Replicate, value = FT) %>% ungroup()
p4 <- Graphing_Corr(cmp, "Flowering Time", "2021-2022") +
  stat_cor(label.y = 105) +
  xlim(65, 145) +
  ylim(65, 145)

outlier_upper <- quantile(PHENO_SINGLE_2022$FT, .75, na.rm = T) + (1.5 * (quantile(PHENO_SINGLE_2022$FT, .75, na.rm = T) - quantile(PHENO_SINGLE_2022$FT, .25, na.rm =T)))
outlier_lower <- quantile(PHENO_SINGLE_2022$FT, .25, na.rm = T) - (1.5 * (quantile(PHENO_SINGLE_2022$FT, .75, na.rm = T) - quantile(PHENO_SINGLE_2022$FT, .25, na.rm =T)))
outlier_data <- PHENO_SINGLE_2022 %>% filter(FT > outlier_upper | FT < outlier_lower)

ggplot(PHENO_SINGLE_2022, aes(x = FT)) +
  geom_histogram(bins = 50) +
  geom_vline(xintercept = outlier_upper, color = "red")


# Residual Plotting
cmp <- na.omit(cmp)
residual_df <- cmp %>% reframe(Residuals = resid(lm(`2` ~ `1`)))
cmp$Residuals <- residual_df$Residuals
Standard_dev <- sd(cmp$Residuals)
cmp <- cmp %>% mutate(Residuals = Residuals/Standard_dev)
tmp_outlier <- cmp %>% filter(Residuals > 2 | Residuals < -2)

d <- ggplot(cmp, aes(x = `1`, y = Residuals)) +
  geom_jitter() +
  geom_hline(yintercept = 0, color = 'red') +
  geom_hline(yintercept = 2, color = "blue") +
  geom_hline(yintercept = -2, color = 'blue') +
  scale_y_continuous(breaks = seq(-4, 4, 1)) +
  geom_text_repel(label = ifelse(cmp$Residuals > 2 | cmp$Residuals < -2,
                                 cmp$Genotype,
                                 ""), size = 3, hjust = .5, max.overlaps = 20) +
  labs(x = "Flowering Time",
       title = "Residual Plot Single Flowering Time 2022")

# Single Correlation 2021-2022 100 SW

cmp <- PHENO_SINGLE_2022 %>% group_by(Genotype, Generation, Replicate) %>% summarise(SEED_WEIGHT_100) %>% spread(key = Replicate, value = SEED_WEIGHT_100) %>% ungroup()
p5 <- Graphing_Corr(cmp, "100 Seed Weight (g)", "2021-2022") +
  stat_cor(label.y = 5.5) +
  ylim(2.5, 7) +
  xlim(2.5,7)

outlier_upper <- quantile(PHENO_SINGLE_2022$SEED_WEIGHT_100, .75, na.rm = T) + (1.5 * (quantile(PHENO_SINGLE_2022$SEED_WEIGHT_100, .75, na.rm = T) - quantile(PHENO_SINGLE_2022$SEED_WEIGHT_100, .25, na.rm =T)))
outlier_lower <- quantile(PHENO_SINGLE_2022$SEED_WEIGHT_100, .25, na.rm = T) - (1.5 * (quantile(PHENO_SINGLE_2022$SEED_WEIGHT_100, .75, na.rm = T) - quantile(PHENO_SINGLE_2022$SEED_WEIGHT_100, .25, na.rm =T)))
outlier_data <- PHENO_SINGLE_2022 %>% filter(SEED_WEIGHT_100 > outlier_upper | SEED_WEIGHT_100 < outlier_lower)

ggplot(PHENO_SINGLE_2022, aes(x = SEED_WEIGHT_100)) +
  geom_histogram(bins = 400) +
  geom_vline(xintercept = outlier_upper, color = "red")

# Residual Plotting
cmp <- na.omit(cmp)
residual_df <- cmp %>% reframe(Residuals = resid(lm(`2` ~ `1`)))
cmp$Residuals <- residual_df$Residuals
Standard_dev <- sd(cmp$Residuals)
cmp <- cmp %>% mutate(Residuals = Residuals/Standard_dev)
tmp_outlier <- cmp %>% filter(Residuals > 2 | Residuals < -2)

e <- ggplot(cmp, aes(x = `1`, y = Residuals)) +
  geom_jitter() +
  geom_hline(yintercept = 0, color = 'red') +
  geom_hline(yintercept = 2, color = "blue") +
  geom_hline(yintercept = -2, color = 'blue') +
  scale_y_continuous(breaks = seq(-4, 4, 1)) +
  geom_text_repel(label = ifelse(cmp$Residuals > 2 | cmp$Residuals < -2,
                                 cmp$Genotype,
                                 ""), size = 3, hjust = .5, max.overlaps = 20) +
  labs(x = "100 Seed Weight",
       title = "Residuals Single 100 Seed Weigh 2022")


# Creating a Dataframe that filters by 2022-2023 season and single condition

PHENO_SINGLE_2023 <- PHENO_FULL %>% filter(Condition == "single" & Exp_year == 2023)

# Single Correlation 2022-2023 TW

cmp <- PHENO_SINGLE_2023 %>% group_by(Genotype, Generation, Replicate) %>% reframe(TOTAL_MASS) %>% spread(key = Replicate, value = TOTAL_MASS) %>% ungroup()
p6 <- Graphing_Corr(cmp, "Total Weight (g)", "2022-2023") +
  stat_cor(label.y = 170) +
  xlim(0, 230) +
  ylim(0, 230)

# Identifying Outliers
PHENO_tmp <- PHENO2023 %>% filter(Condition == "single")

outlier_upper <- quantile(PHENO_tmp$TOTAL_MASS, .75, na.rm = T) + (1.5 * (quantile(PHENO_tmp$TOTAL_MASS, .75, na.rm = T) - quantile(PHENO_tmp$TOTAL_MASS, .25, na.rm =T)))
outlier_lower <- quantile(PHENO_tmp$TOTAL_MASS, .25, na.rm = T) - (1.5 * (quantile(PHENO_tmp$TOTAL_MASS, .75, na.rm = T) - quantile(PHENO_tmp$TOTAL_MASS, .25, na.rm =T)))
outlier_data <- PHENO_tmp %>% filter(TOTAL_MASS > outlier_upper | TOTAL_MASS < outlier_lower)


ggplot(PHENO_SINGLE_2023, aes(x = TOTAL_MASS)) +
  geom_histogram(bins = 100) +
  geom_vline(xintercept = outlier_upper, color = "red") +
  geom_vline(xintercept = outlier_lower, color = 'blue')

# Residual Plotting
cmp <- na.omit(cmp)
residual_df <- cmp %>% reframe(Residuals = resid(lm(`2` ~ `1`)))
cmp$Residuals <- residual_df$Residuals
Standard_dev <- sd(cmp$Residuals)
cmp <- cmp %>% mutate(Residuals = Residuals/Standard_dev)
tmp_outlier <- cmp %>% filter(Residuals > 2 | Residuals < -2)

f <- ggplot(cmp, aes(x = `1`, y = Residuals)) +
  geom_jitter() +
  geom_hline(yintercept = 0, color = 'red') +
  geom_hline(yintercept = 2, color = "blue") +
  geom_hline(yintercept = -2, color = 'blue') +
  scale_y_continuous(breaks = seq(-4, 4, 1)) +
  geom_text_repel(label = ifelse(cmp$Residuals > 2 | cmp$Residuals < -2,
                                 cmp$Genotype,
                                 ""), size = 3, hjust = .8, max.overlaps = 20) +
  labs(x = "Total Weight (grams)",
       title = "Residuals Total Weight 2023")

# Single Correlation 2022-2023 Fecundity

cmp <- PHENO_SINGLE_2023 %>% group_by(Genotype, Generation, Replicate) %>% summarise(FECUNDITY) %>% spread(key = Replicate, value = FECUNDITY) %>% ungroup()
p7 <- Graphing_Corr(cmp, "Relative Fecundity", "2022-2023") +
  stat_cor(label.y = 2.5) +
  xlim(-1.1, 3) +
  ylim(-1.1, 3) +
  geom_text_repel(label = ifelse(cmp$`1` > 3 | cmp$`2` > 3,
                                 cmp$Genotype,
                                 ""), size = 3, hjust =1, max.overlaps = 30)

# Residual Plotting
cmp <- na.omit(cmp)
residual_df <- cmp %>% reframe(Residuals = resid(lm(`2` ~ `1`)))
cmp$Residuals <- residual_df$Residuals
Standard_dev <- sd(cmp$Residuals)
cmp <- cmp %>% mutate(Residuals = Residuals/Standard_dev)
tmp_outlier <- cmp %>% filter(Residuals > 2 | Residuals < -2)

g <- ggplot(cmp, aes(x = `1`, y = Residuals)) +
  geom_jitter() +
  geom_hline(yintercept = 0, color = 'red') +
  geom_hline(yintercept = 2, color = "blue") +
  geom_hline(yintercept = -2, color = 'blue') +
  scale_y_continuous(breaks = seq(-4, 4, 1)) +
  geom_text_repel(label = ifelse(cmp$Residuals > 2 | cmp$Residuals < -2,
                                 cmp$Genotype,
                                 ""), size = 3, hjust = .8, max.overlaps = 20) +
  labs(x = "Centered Fecundity",
       title = "Residuals Single Fecundity 2023")

# Single Correlation 2022-2023 Fitness

cmp <- PHENO_SINGLE_2023 %>% group_by(Genotype, Generation, Replicate) %>% summarise(ABS_FITNESS) %>% spread(key = Replicate, value = ABS_FITNESS) %>% ungroup()
p8 <- Graphing_Corr(cmp, "Absolute Fitness", "2022-2023") +
  stat_cor(label.y = 3) +
  ylim(-1.8, 5) +
  xlim(-1.8, 5)

# Residual Plotting
cmp <- na.omit(cmp)
residual_df <- cmp %>% reframe(Residuals = resid(lm(`2` ~ `1`)))
cmp$Residuals <- residual_df$Residuals
Standard_dev <- sd(cmp$Residuals)
cmp <- cmp %>% mutate(Residuals = Residuals/Standard_dev)
tmp_outlier <- cmp %>% filter(Residuals > 2 | Residuals < -2)

h <- ggplot(cmp, aes(x = `1`, y = Residuals)) +
  geom_jitter() +
  geom_hline(yintercept = 0, color = 'red') +
  geom_hline(yintercept = 2, color = "blue") +
  geom_hline(yintercept = -2, color = 'blue') +
  scale_y_continuous(breaks = seq(-4, 4, 1)) +
  geom_text_repel(label = ifelse(cmp$Residuals > 2 | cmp$Residuals < -2,
                                 cmp$Genotype,
                                 ""), size = 3, hjust = .8, max.overlaps = 20) +
  labs(x = "Centered Absolute Fitness",
       title = "Residuals Single Centered Absolute Fitness 2023")

# Single Correlation 2022-2023 Flowering Time

cmp <- PHENO_SINGLE_2023 %>% group_by(Genotype, Generation, Replicate) %>% summarise(FT) %>% spread(key = Replicate, value = FT) %>% ungroup()
p9 <- Graphing_Corr(cmp, "Flowering Time", "2022-2023") +
  stat_cor(label.y = 125) +
  ylim(65, 145) +
  xlim(65, 145)

# Residual Plotting
cmp <- na.omit(cmp)
residual_df <- cmp %>% reframe(Residuals = resid(lm(`2` ~ `1`)))
cmp$Residuals <- residual_df$Residuals
Standard_dev <- sd(cmp$Residuals)
cmp <- cmp %>% mutate(Residuals = Residuals/Standard_dev)
tmp_outlier <- cmp %>% filter(Residuals > 2 | Residuals < -2)

i <- ggplot(cmp, aes(x = `1`, y = Residuals)) +
  geom_jitter() +
  geom_hline(yintercept = 0, color = 'red') +
  geom_hline(yintercept = 2, color = "blue") +
  geom_hline(yintercept = -2, color = 'blue') +
  scale_y_continuous(breaks = seq(-4, 4, 1)) +
  geom_text_repel(label = ifelse(cmp$Residuals > 2 | cmp$Residuals < -2,
                                 cmp$Genotype,
                                 ""), size = 3, hjust = .8, max.overlaps = 20) +
  labs(x = "Flowering Time",
       title = "Residuals Single Flowering Time 2023")


PHENO_tmp <- PHENO2023 %>% filter(Condition == "single")
outlier_upper <- quantile(PHENO_tmp$FT, .75, na.rm = T) + (1.5 * (quantile(PHENO_tmp$FT, .75, na.rm = T) - quantile(PHENO_tmp$FT, .25, na.rm =T)))
outlier_lower <- quantile(PHENO_tmp$FT, .25, na.rm = T) - (1.5 * (quantile(PHENO_tmp$FT, .75, na.rm = T) - quantile(PHENO_tmp$FT, .25, na.rm =T)))
outlier_data <- PHENO_tmp %>% filter(FT > outlier_upper | FT < outlier_lower)


# Single Correlation 2022-2023 100 SW

cmp <- PHENO_SINGLE_2023 %>% group_by(Genotype, Generation, Replicate) %>% summarise(SEED_WEIGHT_100) %>% spread(key = Replicate, value = SEED_WEIGHT_100) %>% ungroup()
p10 <- Graphing_Corr(cmp, "100 Seed Weight (g)", "2022-2023") +
  stat_cor(label.y = 5.5) +
  ylim(2.5, 7) +
  xlim(2.5,7)

PHENO_tmp <- PHENO2023 %>% filter(Condition == "single")
outlier_upper <- quantile(PHENO_tmp$SEED_WEIGHT_100, .75, na.rm = T) + (1.5 * (quantile(PHENO_tmp$SEED_WEIGHT_100, .75, na.rm = T) - quantile(PHENO_tmp$SEED_WEIGHT_100, .25, na.rm =T)))
outlier_lower <- quantile(PHENO_tmp$SEED_WEIGHT_100, .25, na.rm = T) - (1.5 * (quantile(PHENO_tmp$SEED_WEIGHT_100, .75, na.rm = T) - quantile(PHENO_tmp$SEED_WEIGHT_100, .25, na.rm =T)))
outlier_data <- PHENO_tmp %>% filter(SEED_WEIGHT_100 > outlier_upper | SEED_WEIGHT_100 < outlier_lower)

# Residual Plotting
cmp <- na.omit(cmp)
residual_df <- cmp %>% reframe(Residuals = resid(lm(`2` ~ `1`)))
cmp$Residuals <- residual_df$Residuals
Standard_dev <- sd(cmp$Residuals)
cmp <- cmp %>% mutate(Residuals = Residuals/Standard_dev)
tmp_outlier <- cmp %>% filter(Residuals > 2 | Residuals < -2)

j <- ggplot(cmp, aes(x = `1`, y = Residuals)) +
  geom_jitter() +
  geom_hline(yintercept = 0, color = 'red') +
  geom_hline(yintercept = 2, color = "blue") +
  geom_hline(yintercept = -2, color = 'blue') +
  scale_y_continuous(breaks = seq(-4, 4, 1)) +
  geom_text_repel(label = ifelse(cmp$Residuals > 2 | cmp$Residuals < -2,
                                 cmp$Genotype,
                                 ""), size = 3, hjust = .8, max.overlaps = 20) +
  labs(x = "100 Seed Weight (grams)",
       title = "Residuals Single 100 Seed Weight 2023")

g <- arrangeGrob(a,b,c,d,e,f,g,h,i,j, nrow = 2, ncol = 5)
ggsave("/bigdata/koeniglab/jmarz001/Ag-Competition/results/Residuals_Single_Both_Years.png", g, width = 26, height = 16)

z <- arrangeGrob(p1, p2, p3, p4, p5, p6, p7,p8, p9, p10, nrow = 2, ncol = 5, top = "Correlation of Single Replicates")
ggsave("/bigdata/koeniglab/jmarz001/Ag-Competition/results/Correlation_Plot_Single.png", z, width = 32, height = 16)

############################################

# BETWEEN YEARS

#Creating a function to easily graph all phenos
# x = dataframe
# y = phenotype graphed in quotes (included in title)

new_graph <- function(x, y){
  ggplot(x, aes(`1`, `2`), add = "reg.line") +
    geom_jitter(alpha = .5) +
    geom_abline(slope = 1, intercept = 0, color = 'red') +
    geom_smooth(method = 'lm') +
    labs(x = "Rep 1",
         y = "Rep 2",
         title = paste("Correlation of", sep = " ", y, "Single Replicates"))}

# Averaging the Replicates between the years

tmp <- PHENO_FULL %>% group_by(Genotype, Replicate, Generation, Condition) %>% summarise(across(.cols = where(is.numeric), .fns = mean, na.rm = T)) %>% ungroup

### SINGLE

smp <- tmp %>% filter(Condition == 'single')

# TOTAL WEIGHT

mp <- smp %>% group_by(Genotype, Generation, Replicate) %>% summarise(TOTAL_MASS) %>% spread(key = Replicate, value = TOTAL_MASS) %>% ungroup()

a1 <- new_graph(mp, "Total Weight (grams)") +
  stat_cor(label.y = 170) +
  xlim(0, 210) +
  ylim(0, 210)

# CENTERED FECUNDITY

mp <- smp %>% group_by(Genotype, Generation, Replicate) %>% summarise(FECUNDITY) %>% spread(key = Replicate, value = FECUNDITY) %>% ungroup()

a2 <- new_graph(mp, "Centered Fecundity") +
  stat_cor(label.y = 5, label.x = 1.5) +
  xlim(-1, 8.5) +
  ylim(-1, 8.5)

# ABSOLUTE FITNESS

mp <- smp %>% group_by(Genotype, Generation, Replicate) %>% summarise(ABS_FITNESS) %>% spread(key = Replicate, value = ABS_FITNESS) %>% ungroup()

a3 <- new_graph(mp, "Centered Absolute Fitness") +
  stat_cor() +
  xlim(-2.5, 3) +
  ylim(-2.5, 3)

# FLOWERING TIME

mp <- smp %>% group_by(Genotype, Generation, Replicate) %>% summarise(FT) %>% spread(key = Replicate, value = FT) %>% ungroup()

a4 <- new_graph(mp, "Flowering Time") +
  stat_cor() +
  xlim(92, 125) +
  ylim(92, 125)

# 100 SEED WEIGHT

mp <- smp %>% group_by(Genotype, Generation, Replicate) %>% summarise(SEED_WEIGHT_100) %>% spread(key = Replicate, value = SEED_WEIGHT_100) %>% ungroup()

a5 <- new_graph(mp, "100 Seed Weight") +
  stat_cor() +
  xlim(2.9, 6.6) +
  ylim(2.9, 6.6)

####################################################

new_graph <- function(x, y){
  ggplot(x, aes(`1`, `2`), add = "reg.line") +
    geom_jitter(alpha = .5) +
    geom_abline(slope = 1, intercept = 0, color = 'red') +
    geom_smooth(method = 'lm') +
    labs(x = "Rep 1",
         y = "Rep 2",
         title = paste("Correlation of", sep = " ", y, "Mixed Replicates"))}

### MIXED

smp <- tmp %>% filter(Condition == 'mixed')

# TOTAL WEIGHT
mp <- smp %>% group_by(Genotype, Generation, Replicate) %>% summarise(TOTAL_MASS) %>% spread(key = Replicate, value = TOTAL_MASS) %>% ungroup()

a6 <- new_graph(mp, "Total Weight (grams)") +
  stat_cor(label.y = 170) +
  xlim(0, 210) +
  ylim(0, 210)

# CENTERED FECUNDITY

mp <- smp %>% group_by(Genotype, Generation, Replicate) %>% summarise(FECUNDITY) %>% spread(key = Replicate, value = FECUNDITY) %>% ungroup()

a7 <- new_graph(mp, "Centered Fecundity") +
  stat_cor(label.y = 5, label.x = 1.5) +
  xlim(-1, 8.5) +
  ylim(-1, 8.5)

# ABSOLUTE FITNESS

mp <- smp %>% group_by(Genotype, Generation, Replicate) %>% summarise(ABS_FITNESS) %>% spread(key = Replicate, value = ABS_FITNESS) %>% ungroup()

a8 <- new_graph(mp, "Centered Absolute Fitness") +
  stat_cor() +
  xlim(-2.5, 3) +
  ylim(-2.5, 3)

# FLOWERING TIME

mp <- smp %>% group_by(Genotype, Generation, Replicate) %>% summarise(FT) %>% spread(key = Replicate, value = FT) %>% ungroup()

a9 <- new_graph(mp, "Flowering Time") +
  stat_cor() +
  xlim(92, 125) +
  ylim(92, 125)

# 100 SEED WEIGHT

mp <- smp %>% group_by(Genotype, Generation, Replicate) %>% summarise(SEED_WEIGHT_100) %>% spread(key = Replicate, value = SEED_WEIGHT_100) %>% ungroup()

a10 <- new_graph(mp, "100 Seed Weight") +
  stat_cor() +
  xlim(2.9, 6.6) +
  ylim(2.9, 6.6)

n <- arrangeGrob(a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, top = "Correlation Between Years", nrow = 2, ncol = 5)
ggsave("/bigdata/koeniglab/jmarz001/Ag-Competition/results/Correlation_Between_Years.png", n, width =32, height = 16)

# Distribution of number of plants survived

tok <- PHENO_FULL %>% filter(Plants != 0)
k <- tok$Plants
l <- test$Plants

mean_val <- mean(tok$Plants)
median_val <- median(tok$Plants)

ggplot() +
  geom_histogram(aes(k), alpha = .5) +
  geom_histogram(aes(l), fill = 'red', alpha = .3) +
  scale_y_continuous(breaks = seq(0, 1000, 50)) +
  scale_x_continuous(breaks = seq(0, 12, 1)) +
  geom_vline(aes(xintercept = mean_val, color = 'mean')) +
  geom_vline(aes(xintercept = median_val, color = 'median')) +
  labs(colour = "Mean and Median",
       x = "number of plants",
       title = "Distribution of number of plants")
ggsave('/bigdata/koeniglab/jmarz001/Ag-Competition/data/Distribution_Plants.pdf', width = 12, height = 10)
