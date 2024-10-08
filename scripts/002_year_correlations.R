#!/usr/bin/env Rscript
#SBATCH --ntasks=1
#SBATCH --mem=50G
#SBATCH --time=02:00:00
#SBATCH --job-name="Ag-Competition"
#SBATCH --output=/rhome/jmarz001/bigdata/Ag-Competition/scripts/002_year_correlations.stdout
#SBATCH -p short


library(tidyverse)
library(gridExtra)
library(car)
library(ggrepel)
library(data.table)

setwd("/rhome/jmarz001/bigdata/Ag-Competition/data")
df <- fread("JOINED_PHENOTYPES.tsv")






# average phenotypes across replicates 
## wait that shouldn"t average anything then, get rid of Replicate grouping?
df_mean <- df %>% group_by(Genotype, Condition, Exp_year) %>% summarise(across(where(is.numeric), \(x) mean(x, na.rm=T)))



#Creating a function to easily graph all phenos
# x = dataframe
# y = phenotype graphed in quotes (included in title)

new_graph <- function(x, y){
  ggplot(x, aes(`1`, `2`), add = "reg.line") +
    geom_jitter(alpha = .5) +
    geom_abline(slope = 1, intercept = 0, color = "red") +
    geom_smooth(method = "lm") +
    labs(x = "1",
         y = "2",
         title = paste("Correlation of", sep = " ", y, "Single Replicates"))}

# Averaging the Replicates between the years

tmp <- df %>% group_by(Genotype, Replicate, Generation, Condition) %>% summarise(across(.cols = where(is.numeric), .fns = mean, na.rm = T)) %>% ungroup

### SINGLE

smp <- tmp %>% filter(Condition == "single")

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


new_graph <- function(x, y){
  ggplot(x, aes(`1`, `2`), add = "reg.line") +
    geom_jitter(alpha = .5) +
    geom_abline(slope = 1, intercept = 0, color = "red") +
    geom_smooth(method = "lm") +
    labs(x = "1",
         y = "2",
         title = paste("Correlation of", sep = " ", y, "Mixed Replicates"))}









### MIXED

smp <- tmp %>% filter(Condition == "mixed")

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

