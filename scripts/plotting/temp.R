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

### (2022 rep 1 vs. 2023 rep 1)

#TW

tmp <- PHENO_FULL %>% filter(Condition == "single" & Replicate == 1)
tmp <- tmp %>% group_by(Genotype, Replicate, Exp_year) %>% summarise(TOTAL_MASS) %>% spread(key = Exp_year, value = TOTAL_MASS)

a1 <- ggplot(tmp, aes(`2022`, `2023`), add = 'reg.line') +
  geom_jitter() +
  geom_smooth(method = 'lm') +
  stat_cor() +
  labs(x = "2022 rep 1",
       y = "2023 rep 1",
       title = "Total Weight") +
  xlim(0, 230) +
  ylim(0, 230)

# FEC

tmp <- PHENO_FULL %>% filter(Condition == "single" & Replicate == 1)
tmp <- tmp %>% group_by(Genotype, Replicate, Exp_year) %>% summarise(FECUNDITY) %>% spread(key = Exp_year, value = FECUNDITY)

a2 <- ggplot(tmp, aes(`2022`, `2023`), add = 'reg.line') +
  geom_jitter() +
  geom_smooth(method = 'lm') +
  stat_cor() +
  labs(x = "2022 rep 1",
       y = "2023 rep 1",
       title = "Centered Fecundity") +
  xlim(-1, 12) +
  ylim(-1, 12)

# Flower

tmp <- PHENO_FULL %>% filter(Condition == "single" & Replicate == 1)
tmp <- tmp %>% group_by(Genotype, Replicate, Exp_year) %>% summarise(FT) %>% spread(key = Exp_year, value = FT)

a3 <- ggplot(tmp, aes(`2022`, `2023`), add = 'reg.line') +
  geom_jitter() +
  geom_smooth(method = 'lm') +
  stat_cor() +
  labs(x = "2022 rep 1",
       y = "2023 rep 1",
       title = "Flowering Time") +
  xlim(65, 145) +
  ylim(65, 145)

# FIT

tmp <- PHENO_FULL %>% filter(Condition == "single" & Replicate == 1)
tmp <- tmp %>% group_by(Genotype, Replicate, Exp_year) %>% summarise(ABS_FITNESS) %>% spread(key = Exp_year, value = ABS_FITNESS)

a4 <- ggplot(tmp, aes(`2022`, `2023`), add = 'reg.line') +
  geom_jitter() +
  geom_smooth(method = 'lm') +
  stat_cor() +
  labs(x = "2022 rep 1",
       y = "2023 rep 1",
       title = "Centered Absolute Fitness") +
  ylim(-1.8, 5) +
  xlim(-1.8, 5)

# 100 SW

tmp <- PHENO_FULL %>% filter(Condition == "single" & Replicate == 1)
tmp <- tmp %>% group_by(Genotype, Replicate, Exp_year) %>% summarise(SEED_WEIGHT_100) %>% spread(key = Exp_year, value = SEED_WEIGHT_100)

a5 <- ggplot(tmp, aes(`2022`, `2023`), add = 'reg.line') +
  geom_jitter() +
  geom_smooth(method = 'lm') +
  stat_cor() +
  labs(x = "2022 rep 1",
       y = "2023 rep 1",
       title = "100 Seed Weight") +
  ylim(2.5, 7) +
  xlim(2.5,7)

###################################################
### (2022 rep 2 vs 2023 rep 2)

# TW

tmp <- PHENO_FULL %>% filter(Condition == "single" & Replicate == 2)
tmp <- tmp %>% group_by(Genotype, Replicate, Exp_year) %>% summarise(TOTAL_MASS) %>% spread(key = Exp_year, value = TOTAL_MASS)

a6 <- ggplot(tmp, aes(`2022`, `2023`), add = 'reg.line') +
  geom_jitter() +
  geom_smooth(method = 'lm') +
  stat_cor() +
  labs(x = '2022 rep 2', y = "2023 rep 2", title = "Total Weight") +
  xlim(0, 230) +
  ylim(0, 230)

# FEC

tmp <- PHENO_FULL %>% filter(Condition == "single" & Replicate == 2)
tmp <- tmp %>% group_by(Genotype, Replicate, Exp_year) %>% summarise(FECUNDITY) %>% spread(key = Exp_year, value = FECUNDITY)

a7 <- ggplot(tmp, aes(`2022`, `2023`), add = 'reg.line') +
  geom_jitter() +
  geom_smooth(method = 'lm') +
  stat_cor() +
  labs(x = '2022 rep 2', y = "2023 rep 2", title = "Centered Fecundity") +
  xlim(-1, 12) +
  ylim(-1, 12)

# FLOWER

tmp <- PHENO_FULL %>% filter(Condition == "single" & Replicate == 2)
tmp <- tmp %>% group_by(Genotype, Replicate, Exp_year) %>% summarise(FT) %>% spread(key = Exp_year, value = FT)

a8 <- ggplot(tmp, aes(`2022`, `2023`), add = 'reg.line') +
  geom_jitter() +
  geom_smooth(method = 'lm') +
  stat_cor() +
  labs(x = '2022 rep 2', y = "2023 rep 2", title = "Flowering Time") +
  xlim(65, 145) +
  ylim(65, 145)

# FITNESS

tmp <- PHENO_FULL %>% filter(Condition == "single" & Replicate == 2)
tmp <- tmp %>% group_by(Genotype, Replicate, Exp_year) %>% summarise(ABS_FITNESS) %>% spread(key = Exp_year, value = ABS_FITNESS)

a9 <- ggplot(tmp, aes(`2022`, `2023`), add = 'reg.line') +
  geom_jitter() +
  geom_smooth(method = 'lm') +
  stat_cor() +
  labs(x = '2022 rep 2', y = "2023 rep 2", title = "Centered Absolute Fitness") +
  ylim(-1.8, 5) +
  xlim(-1.8, 5)

# 100 SW

tmp <- PHENO_FULL %>% filter(Condition == "single" & Replicate == 2)
tmp <- tmp %>% group_by(Genotype, Replicate, Exp_year) %>% summarise(SEED_WEIGHT_100) %>% spread(key = Exp_year, value = SEED_WEIGHT_100)

a10 <- ggplot(tmp, aes(`2022`, `2023`), add = 'reg.line') +
  geom_jitter() +
  geom_smooth(method = 'lm') +
  stat_cor() +
  labs(x = '2022 rep 2', y = "2023 rep 2", title = "100 Seed Weight") +
  ylim(2.5, 7) +
  xlim(2.5,7)

############################################################################
# Rep 1 2022 vs Rep 2 2023

# TW
tmp <- PHENO_FULL %>% filter(Condition == 'single') %>% group_by(Genotype, Replicate, Exp_year) %>% summarise(TOTAL_MASS)
tmp <- tmp %>% filter(Replicate == 1 & Exp_year == 2022 | Replicate == 2 & Exp_year ==2023)
O <- subset(tmp, Exp_year == 2022)
O <- O %>% filter(Genotype != '2_168')
O1 <- subset(tmp, Exp_year == 2023)
tmp <- cbind(O, O1)

a11 <- ggplot(tmp, aes(TOTAL_MASS...4, TOTAL_MASS...8), add = "reg.line") +
  geom_jitter() +
  geom_smooth(method = "lm") +
  stat_cor() +
  labs(x = "Rep 1 2022", y = "Rep 2 2023", title = "Total Weight") +
  xlim(0, 230) +
  ylim(0, 230)

# FEC

tmp <- PHENO_FULL %>% filter(Condition == 'single') %>% group_by(Genotype, Replicate, Exp_year) %>% summarise(FECUNDITY)
tmp <- tmp %>% filter(Replicate == 1 & Exp_year == 2022 | Replicate == 2 & Exp_year ==2023)
O <- subset(tmp, Exp_year == 2022)
O <- O %>% filter(Genotype != '2_168')
O1 <- subset(tmp, Exp_year == 2023)
tmp <- cbind(O, O1)

a12 <- ggplot(tmp, aes(FECUNDITY...4, FECUNDITY...8), add = "reg.line") +
  geom_jitter() +
  geom_smooth(method = "lm") +
  stat_cor() +
  labs(x = "Rep 1 2022", y = "Rep 2 2023", title = "Centered Fecundity") +
  xlim(-1, 12) +
  ylim(-1, 12)

# FLOWER

tmp <- PHENO_FULL %>% filter(Condition == 'single') %>% group_by(Genotype, Replicate, Exp_year) %>% summarise(FT)
tmp <- tmp %>% filter(Replicate == 1 & Exp_year == 2022 | Replicate == 2 & Exp_year ==2023)
O <- subset(tmp, Exp_year == 2022)
O <- O %>% filter(Genotype != '2_168')
O1 <- subset(tmp, Exp_year == 2023)
tmp <- cbind(O, O1)

a13 <- ggplot(tmp, aes(FT...4, FT...8), add = "reg.line") +
  geom_jitter() +
  geom_smooth(method = "lm") +
  stat_cor() +
  labs(x = "Rep 1 2022", y = "Rep 2 2023", title = "Flowering Time") +
  xlim(65, 145) +
  ylim(65, 145)

# FIT

tmp <- PHENO_FULL %>% filter(Condition == 'single') %>% group_by(Genotype, Replicate, Exp_year) %>% summarise(ABS_FITNESS)
tmp <- tmp %>% filter(Replicate == 1 & Exp_year == 2022 | Replicate == 2 & Exp_year ==2023)
O <- subset(tmp, Exp_year == 2022)
O <- O %>% filter(Genotype != '2_168')
O1 <- subset(tmp, Exp_year == 2023)
tmp <- cbind(O, O1)

a14 <- ggplot(tmp, aes(ABS_FITNESS...4, ABS_FITNESS...8), add = "reg.line") +
  geom_jitter() +
  geom_smooth(method = "lm") +
  stat_cor() +
  labs(x = "Rep 1 2022", y = "Rep 2 2023", title = "Centered Absolute Fitness") +
  ylim(-1.8, 5) +
  xlim(-1.8, 5)

# 100 SW

tmp <- PHENO_FULL %>% filter(Condition == 'single') %>% group_by(Genotype, Replicate, Exp_year) %>% summarise(SEED_WEIGHT_100)
tmp <- tmp %>% filter(Replicate == 1 & Exp_year == 2022 | Replicate == 2 & Exp_year ==2023)
O <- subset(tmp, Exp_year == 2022)
O <- O %>% filter(Genotype != '2_168')
O1 <- subset(tmp, Exp_year == 2023)
tmp <- cbind(O, O1)

a15 <- ggplot(tmp, aes(`100_seed_weight...4`, `100_seed_weight...8`), add = "reg.line") +
  geom_jitter() +
  geom_smooth(method = "lm") +
  stat_cor() +
  labs(x = "Rep 1 2022", y = "Rep 2 2023", title = '100 Seed Weight') +
  ylim(2.5, 7) +
  xlim(2.5,7)

################################################################################

# Rep 2 2022 vs. Rep 1 2023

# TW
tmp <- PHENO_FULL %>% filter(Condition == 'single') %>% group_by(Genotype, Replicate, Exp_year) %>% summarise(TOTAL_MASS)
tmp <- tmp %>% filter(Replicate == 2 & Exp_year == 2022 | Replicate == 1 & Exp_year ==2023)
O <- subset(tmp, Exp_year == 2022)
O <- O %>% filter(Genotype != '2_168')
O1 <- subset(tmp, Exp_year == 2023)
tmp <- cbind(O, O1)

a16 <- ggplot(tmp, aes(TOTAL_MASS...4, TOTAL_MASS...8), add = "reg.line") +
  geom_jitter() +
  geom_smooth(method = "lm") +
  stat_cor() +
  labs(x = "Rep 2 2022", y = "Rep 1 2023", title = "Total Weight") +
  xlim(0, 230) +
  ylim(0, 230)

# FEC

tmp <- PHENO_FULL %>% filter(Condition == 'single') %>% group_by(Genotype, Replicate, Exp_year) %>% summarise(FECUNDITY)
tmp <- tmp %>% filter(Replicate == 2 & Exp_year == 2022 | Replicate == 1 & Exp_year ==2023)
O <- subset(tmp, Exp_year == 2022)
O <- O %>% filter(Genotype != '2_168')
O1 <- subset(tmp, Exp_year == 2023)
tmp <- cbind(O, O1)

a17 <- ggplot(tmp, aes(FECUNDITY...4, FECUNDITY...8), add = "reg.line") +
  geom_jitter() +
  geom_smooth(method = "lm") +
  stat_cor() +
  labs(x = "Rep 2 2022", y = "Rep 1 2023", title = "Centered Fecundity") +
  xlim(-1, 12) +
  ylim(-1, 12)

# FLOWER

tmp <- PHENO_FULL %>% filter(Condition == 'single') %>% group_by(Genotype, Replicate, Exp_year) %>% summarise(FT)
tmp <- tmp %>% filter(Replicate == 2 & Exp_year == 2022 | Replicate == 1 & Exp_year ==2023)
O <- subset(tmp, Exp_year == 2022)
O <- O %>% filter(Genotype != '2_168')
O1 <- subset(tmp, Exp_year == 2023)
tmp <- cbind(O, O1)

a18 <- ggplot(tmp, aes(FT...4, FT...8), add = "reg.line") +
  geom_jitter() +
  geom_smooth(method = "lm") +
  stat_cor() +
  labs(x = "Rep 2 2022", y = "Rep 1 2023", title = "Flowering Time") +
  xlim(65, 145) +
  ylim(65, 145)

# FIT

tmp <- PHENO_FULL %>% filter(Condition == 'single') %>% group_by(Genotype, Replicate, Exp_year) %>% summarise(ABS_FITNESS)
tmp <- tmp %>% filter(Replicate == 2 & Exp_year == 2022 | Replicate == 1 & Exp_year ==2023)
O <- subset(tmp, Exp_year == 2022)
O <- O %>% filter(Genotype != '2_168')
O1 <- subset(tmp, Exp_year == 2023)
tmp <- cbind(O, O1)

a19 <- ggplot(tmp, aes(ABS_FITNESS...4, ABS_FITNESS...8), add = "reg.line") +
  geom_jitter() +
  geom_smooth(method = "lm") +
  stat_cor() +
  labs(x = "Rep 2 2022", y = "Rep 1 2023", title = "Centered Absolute Fitness") +
  ylim(-1.8, 5) +
  xlim(-1.8, 5)

# 100 SW

tmp <- PHENO_FULL %>% filter(Condition == 'single') %>% group_by(Genotype, Replicate, Exp_year) %>% summarise(SEED_WEIGHT_100)
tmp <- tmp %>% filter(Replicate == 2 & Exp_year == 2022 | Replicate == 1 & Exp_year ==2023)
O <- subset(tmp, Exp_year == 2022)
O <- O %>% filter(Genotype != '2_168')
O1 <- subset(tmp, Exp_year == 2023)
tmp <- cbind(O, O1)

a20 <- ggplot(tmp, aes(`100_seed_weight...4`, `100_seed_weight...8`), add = "reg.line") +
  geom_jitter() +
  geom_smooth(method = "lm") +
  stat_cor() +
  labs(x = "Rep 2 2022", y = "Rep 1 2023", title = "100 Seed Weight") +
  ylim(2.5, 7) +
  xlim(2.5,7)

mew <- arrangeGrob(a1, a6, a11, a16, a2, a7, a12, a17, a3, a8, a13, a18, a4, a9, a14, a19, a5, a10, a15, a20, top = "Single Correlations All Combos", nrow = 5, ncol = 4)
ggsave('/bigdata/koeniglab/jmarz001/Ag-Competition/data/Single_Corr_All_Combos.pdf', mew, width = 42, height = 26)


########################################################################################

## MIXED


### (2022 rep 1 vs. 2023 rep 1)

#TW

tmp <- PHENO_FULL %>% filter(Condition == "mixed" & Replicate == 1)
tmp <- tmp %>% group_by(Genotype, Replicate, Exp_year) %>% summarise(TOTAL_MASS) %>% spread(key = Exp_year, value = TOTAL_MASS)

a1 <- ggplot(tmp, aes(`2022`, `2023`), add = 'reg.line') +
  geom_jitter() +
  geom_smooth(method = 'lm') +
  stat_cor() +
  labs(x = "2022 rep 1",
       y = "2023 rep 1",
       title = "Total Weight") +
  xlim(0, 250) +
  ylim(0, 250)

# FEC

tmp <- PHENO_FULL %>% filter(Condition == "mixed" & Replicate == 1)
tmp <- tmp %>% group_by(Genotype, Replicate, Exp_year) %>% summarise(FECUNDITY) %>% spread(key = Exp_year, value = FECUNDITY)

a2 <- ggplot(tmp, aes(`2022`, `2023`), add = 'reg.line') +
  geom_jitter() +
  geom_smooth(method = 'lm') +
  stat_cor() +
  labs(x = "2022 rep 1",
       y = "2023 rep 1",
       title = "Centered Fecundity") +
  ylim(-1.5, 15)+
  xlim(-1.5, 15)

# Flower

tmp <- PHENO_FULL %>% filter(Condition == "mixed" & Replicate == 1)
tmp <- tmp %>% group_by(Genotype, Replicate, Exp_year) %>% summarise(FT) %>% spread(key = Exp_year, value = FT)

a3 <- ggplot(tmp, aes(`2022`, `2023`), add = 'reg.line') +
  geom_jitter() +
  geom_smooth(method = 'lm') +
  stat_cor() +
  labs(x = "2022 rep 1",
       y = "2023 rep 1",
       title = "Flowering Time") +
  xlim(75, 135) +
  ylim(75, 135)

# FIT

tmp <- PHENO_FULL %>% filter(Condition == "mixed" & Replicate == 1)
tmp <- tmp %>% group_by(Genotype, Replicate, Exp_year) %>% summarise(ABS_FITNESS) %>% spread(key = Exp_year, value = ABS_FITNESS)

a4 <- ggplot(tmp, aes(`2022`, `2023`), add = 'reg.line') +
  geom_jitter() +
  geom_smooth(method = 'lm') +
  stat_cor() +
  labs(x = "2022 rep 1",
       y = "2023 rep 1",
       title = "Centered Absolute Fitness") +
  xlim(-2.5, 4.2) +
  ylim(-2.5, 4.2)

# 100 SW

tmp <- PHENO_FULL %>% filter(Condition == "mixed" & Replicate == 1)
tmp <- tmp %>% group_by(Genotype, Replicate, Exp_year) %>% summarise(SEED_WEIGHT_100) %>% spread(key = Exp_year, value = SEED_WEIGHT_100)

a5 <- ggplot(tmp, aes(`2022`, `2023`), add = 'reg.line') +
  geom_jitter() +
  geom_smooth(method = 'lm') +
  stat_cor() +
  labs(x = "2022 rep 1",
       y = "2023 rep 1",
       title = "100 Seed Weight") +
  xlim(2.4, 6.1) +
  ylim(2.4, 6.1)

###################################################
### (2022 rep 2 vs 2023 rep 2)

# TW

tmp <- PHENO_FULL %>% filter(Condition == "mixed" & Replicate == 2)
tmp <- tmp %>% group_by(Genotype, Replicate, Exp_year) %>% summarise(TOTAL_MASS) %>% spread(key = Exp_year, value = TOTAL_MASS)

a6 <- ggplot(tmp, aes(`2022`, `2023`), add = 'reg.line') +
  geom_jitter() +
  geom_smooth(method = 'lm') +
  stat_cor() +
  labs(x = '2022 rep 2', y = "2023 rep 2", title = "Total Weight") +
  xlim(0, 250) +
  ylim(0, 250)

# FEC

tmp <- PHENO_FULL %>% filter(Condition == "mixed" & Replicate == 2)
tmp <- tmp %>% group_by(Genotype, Replicate, Exp_year) %>% summarise(FECUNDITY) %>% spread(key = Exp_year, value = FECUNDITY)

a7 <- ggplot(tmp, aes(`2022`, `2023`), add = 'reg.line') +
  geom_jitter() +
  geom_smooth(method = 'lm') +
  stat_cor() +
  labs(x = '2022 rep 2', y = "2023 rep 2", title = "Centered Fecundity") +
  ylim(-1.5, 15)+
  xlim(-1.5, 15)

# FLOWER

tmp <- PHENO_FULL %>% filter(Condition == "mixed" & Replicate == 2)
tmp <- tmp %>% group_by(Genotype, Replicate, Exp_year) %>% summarise(FT) %>% spread(key = Exp_year, value = FT)

a8 <- ggplot(tmp, aes(`2022`, `2023`), add = 'reg.line') +
  geom_jitter() +
  geom_smooth(method = 'lm') +
  stat_cor() +
  labs(x = '2022 rep 2', y = "2023 rep 2", title = "Flowering Time") +
  xlim(75, 135) +
  ylim(75, 135)

# FITNESS

tmp <- PHENO_FULL %>% filter(Condition == "mixed" & Replicate == 2)
tmp <- tmp %>% group_by(Genotype, Replicate, Exp_year) %>% summarise(ABS_FITNESS) %>% spread(key = Exp_year, value = ABS_FITNESS)

a9 <- ggplot(tmp, aes(`2022`, `2023`), add = 'reg.line') +
  geom_jitter() +
  geom_smooth(method = 'lm') +
  stat_cor() +
  labs(x = '2022 rep 2', y = "2023 rep 2", title = "Centered Absolute Fitness") +
  xlim(-2.5, 4.2) +
  ylim(-2.5, 4.2)

# 100 SW

tmp <- PHENO_FULL %>% filter(Condition == "mixed" & Replicate == 2)
tmp <- tmp %>% group_by(Genotype, Replicate, Exp_year) %>% summarise(SEED_WEIGHT_100) %>% spread(key = Exp_year, value = SEED_WEIGHT_100)

a10 <- ggplot(tmp, aes(`2022`, `2023`), add = 'reg.line') +
  geom_jitter() +
  geom_smooth(method = 'lm') +
  stat_cor() +
  labs(x = '2022 rep 2', y = "2023 rep 2", title = "100 Seed Weight") +
  xlim(2.4, 6.1) +
  ylim(2.4, 6.1)

############################################################################
# Rep 1 2022 vs Rep 2 2023

# TW
tmp <- PHENO_FULL %>% filter(Condition == 'mixed') %>% group_by(Genotype, Replicate, Exp_year) %>% summarise(TOTAL_MASS)
tmp <- tmp %>% filter(Replicate == 1 & Exp_year == 2022 | Replicate == 2 & Exp_year ==2023)
O <- subset(tmp, Exp_year == 2022)
O <- O %>% filter(Genotype != '2_168')
O1 <- subset(tmp, Exp_year == 2023)
tmp <- cbind(O, O1)

a11 <- ggplot(tmp, aes(TOTAL_MASS...4, TOTAL_MASS...8), add = "reg.line") +
  geom_jitter() +
  geom_smooth(method = "lm") +
  stat_cor() +
  labs(x = "Rep 1 2022", y = "Rep 2 2023", title = "Total Weight") +
  xlim(0, 250) +
  ylim(0, 250)

# FEC

tmp <- PHENO_FULL %>% filter(Condition == 'mixed') %>% group_by(Genotype, Replicate, Exp_year) %>% summarise(FECUNDITY)
tmp <- tmp %>% filter(Replicate == 1 & Exp_year == 2022 | Replicate == 2 & Exp_year ==2023)
O <- subset(tmp, Exp_year == 2022)
O <- O %>% filter(Genotype != '2_168')
O1 <- subset(tmp, Exp_year == 2023)
tmp <- cbind(O, O1)

a12 <- ggplot(tmp, aes(FECUNDITY...4, FECUNDITY...8), add = "reg.line") +
  geom_jitter() +
  geom_smooth(method = "lm") +
  stat_cor() +
  labs(x = "Rep 1 2022", y = "Rep 2 2023", title = "Centered Fecundity") +
  ylim(-1.5, 15)+
  xlim(-1.5, 15)

# FLOWER

tmp <- PHENO_FULL %>% filter(Condition == 'mixed') %>% group_by(Genotype, Replicate, Exp_year) %>% summarise(FT)
tmp <- tmp %>% filter(Replicate == 1 & Exp_year == 2022 | Replicate == 2 & Exp_year ==2023)
O <- subset(tmp, Exp_year == 2022)
O <- O %>% filter(Genotype != '2_168')
O1 <- subset(tmp, Exp_year == 2023)
tmp <- cbind(O, O1)

a13 <- ggplot(tmp, aes(FT...4, FT...8), add = "reg.line") +
  geom_jitter() +
  geom_smooth(method = "lm") +
  stat_cor() +
  labs(x = "Rep 1 2022", y = "Rep 2 2023", title = "Flowering Time") +
  xlim(75, 135) +
  ylim(75, 135)

# FIT

tmp <- PHENO_FULL %>% filter(Condition == 'mixed') %>% group_by(Genotype, Replicate, Exp_year) %>% summarise(ABS_FITNESS)
tmp <- tmp %>% filter(Replicate == 1 & Exp_year == 2022 | Replicate == 2 & Exp_year ==2023)
O <- subset(tmp, Exp_year == 2022)
O <- O %>% filter(Genotype != '2_168')
O1 <- subset(tmp, Exp_year == 2023)
tmp <- cbind(O, O1)

a14 <- ggplot(tmp, aes(ABS_FITNESS...4, ABS_FITNESS...8), add = "reg.line") +
  geom_jitter() +
  geom_smooth(method = "lm") +
  stat_cor() +
  labs(x = "Rep 1 2022", y = "Rep 2 2023", title = "Centered Absolute Fitness") +
  xlim(-2.5, 4.2) +
  ylim(-2.5, 4.2)

# 100 SW

tmp <- PHENO_FULL %>% filter(Condition == 'mixed') %>% group_by(Genotype, Replicate, Exp_year) %>% summarise(SEED_WEIGHT_100)
tmp <- tmp %>% filter(Replicate == 1 & Exp_year == 2022 | Replicate == 2 & Exp_year ==2023)
O <- subset(tmp, Exp_year == 2022)
O <- O %>% filter(Genotype != '2_168')
O1 <- subset(tmp, Exp_year == 2023)
tmp <- cbind(O, O1)

a15 <- ggplot(tmp, aes(`100_seed_weight...4`, `100_seed_weight...8`), add = "reg.line") +
  geom_jitter() +
  geom_smooth(method = "lm") +
  stat_cor() +
  labs(x = "Rep 1 2022", y = "Rep 2 2023", title = '100 Seed Weight') +
  xlim(2.4, 6.1) +
  ylim(2.4, 6.1)

################################################################################

# Rep 2 2022 vs. Rep 1 2023

# TW
tmp <- PHENO_FULL %>% filter(Condition == 'mixed') %>% group_by(Genotype, Replicate, Exp_year) %>% summarise(TOTAL_MASS)
tmp <- tmp %>% filter(Replicate == 2 & Exp_year == 2022 | Replicate == 1 & Exp_year ==2023)
O <- subset(tmp, Exp_year == 2022)
O <- O %>% filter(Genotype != '2_168')
O1 <- subset(tmp, Exp_year == 2023)
tmp <- cbind(O, O1)

a16 <- ggplot(tmp, aes(TOTAL_MASS...4, TOTAL_MASS...8), add = "reg.line") +
  geom_jitter() +
  geom_smooth(method = "lm") +
  stat_cor() +
  labs(x = "Rep 2 2022", y = "Rep 1 2023", title = "Total Weight") +
  xlim(0, 250) +
  ylim(0, 250)

# FEC

tmp <- PHENO_FULL %>% filter(Condition == 'mixed') %>% group_by(Genotype, Replicate, Exp_year) %>% summarise(FECUNDITY)
tmp <- tmp %>% filter(Replicate == 2 & Exp_year == 2022 | Replicate == 1 & Exp_year ==2023)
O <- subset(tmp, Exp_year == 2022)
O <- O %>% filter(Genotype != '2_168')
O1 <- subset(tmp, Exp_year == 2023)
tmp <- cbind(O, O1)

a17 <- ggplot(tmp, aes(FECUNDITY...4, FECUNDITY...8), add = "reg.line") +
  geom_jitter() +
  geom_smooth(method = "lm") +
  stat_cor() +
  labs(x = "Rep 2 2022", y = "Rep 1 2023", title = "Centered Fecundity") +
  ylim(-1.5, 15)+
  xlim(-1.5, 15)

# FLOWER

tmp <- PHENO_FULL %>% filter(Condition == 'mixed') %>% group_by(Genotype, Replicate, Exp_year) %>% summarise(FT)
tmp <- tmp %>% filter(Replicate == 2 & Exp_year == 2022 | Replicate == 1 & Exp_year ==2023)
O <- subset(tmp, Exp_year == 2022)
O <- O %>% filter(Genotype != '2_168')
O1 <- subset(tmp, Exp_year == 2023)
tmp <- cbind(O, O1)

a18 <- ggplot(tmp, aes(FT...4, FT...8), add = "reg.line") +
  geom_jitter() +
  geom_smooth(method = "lm") +
  stat_cor() +
  labs(x = "Rep 2 2022", y = "Rep 1 2023", title = "Flowering Time") +
  xlim(75, 135) +
  ylim(75, 135)

# FIT

tmp <- PHENO_FULL %>% filter(Condition == 'mixed') %>% group_by(Genotype, Replicate, Exp_year) %>% summarise(ABS_FITNESS)
tmp <- tmp %>% filter(Replicate == 2 & Exp_year == 2022 | Replicate == 1 & Exp_year ==2023)
O <- subset(tmp, Exp_year == 2022)
O <- O %>% filter(Genotype != '2_168')
O1 <- subset(tmp, Exp_year == 2023)
tmp <- cbind(O, O1)

a19 <- ggplot(tmp, aes(ABS_FITNESS...4, ABS_FITNESS...8), add = "reg.line") +
  geom_jitter() +
  geom_smooth(method = "lm") +
  stat_cor() +
  labs(x = "Rep 2 2022", y = "Rep 1 2023", title = "Centered Absolute Fitness") +
  xlim(-2.5, 4.2) +
  ylim(-2.5, 4.2)

# 100 SW

tmp <- PHENO_FULL %>% filter(Condition == 'mixed') %>% group_by(Genotype, Replicate, Exp_year) %>% summarise(SEED_WEIGHT_100)
tmp <- tmp %>% filter(Replicate == 2 & Exp_year == 2022 | Replicate == 1 & Exp_year ==2023)
O <- subset(tmp, Exp_year == 2022)
O <- O %>% filter(Genotype != '2_168')
O1 <- subset(tmp, Exp_year == 2023)
tmp <- cbind(O, O1)

a20 <- ggplot(tmp, aes(`100_seed_weight...4`, `100_seed_weight...8`), add = "reg.line") +
  geom_jitter() +
  geom_smooth(method = "lm") +
  stat_cor() +
  labs(x = "Rep 2 2022", y = "Rep 1 2023", title = "100 Seed Weight") +
  xlim(2.4, 6.1) +
  ylim(2.4, 6.1)

mew <- arrangeGrob(a1, a6, a11, a16, a2, a7, a12, a17, a3, a8, a13, a18, a4, a9, a14, a19, a5, a10, a15, a20, top = "Mixed Correlations All Combos", nrow = 5, ncol = 4)
ggsave('/bigdata/koeniglab/jmarz001/Ag-Competition/data/Mixed_Corr_All_Combos.pdf', mew, width = 42, height = 26)
