#!/usr/bin/env Rscript

#SBATCH --ntasks=1
#SBATCH --mem=30G
#SBATCH --time=02:00:00
#SBATCH --output=/rhome/jmarz001/bigdata/Ag-Competition/Comparing_Haplotypes.stdout
#SBATCH -p koeniglab

library(tidyverse)
library(readr)
library(dplyr)
library(ggpubr)
library(ggplot2)
library(googlesheets4)
library(tidyr)

# 2021-2022

Single_2021_2022 <- PHENO_FULL_AVERAGE %>% filter(Condition == 'single' & Exp_year == 2022)

### 5a_Comparing_TW_Between_Haplotypes.R

Single_2021_2022 <- Single_2021_2022 %>% filter(Haplotype != "NA")
Single_2021_2022$Haplotype <- as.factor(Single_2021_2022$Haplotype)
Haplo_graph <- Single_2021_2022 %>% group_by(Generation, Haplotype) %>% summarise(Avg_TW = mean(total_seed_mass_g, na.rm = T)) %>% ungroup()
Haplo_graph <- Haplo_graph %>% group_by(Generation) %>% mutate(GenAvg = mean(Avg_TW, na.rm = T)) %>% ungroup()
Haplo_graph$`Avg Yield > Pop & Gen Avg Yield` <- ifelse(Haplo_graph$Avg_TW > Haplo_graph$GenAvg & Haplo_graph$Avg_TW > mean(Haplo_graph$Avg_TW, na.rm = T), TRUE, FALSE)
f_AY <- Haplo_graph %>% filter(`Avg Yield > Pop & Gen Avg Yield` == TRUE) %>% filter(!duplicated(Haplotype)) %>% select(Haplotype)

ggplot(Haplo_graph,aes(x= Haplotype, y = Avg_TW, fill = `Avg Yield > Pop & Gen Avg Yield`)) +
  geom_bar(stat = 'identity') +
  geom_hline(yintercept = mean(Haplo_graph$Avg_TW, na.rm = T), color = "green") +
  geom_hline(aes(yintercept = GenAvg), color = 'blue', linetype = 2) +
  labs(y = "Average Total Seed Weight",
       title = "Average Total Seed Weight between Haplotypes 2021-2022") +
  scale_y_continuous(breaks = seq(0, 200, 25)) +
  facet_wrap(~Generation, scales = "free_x") +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_fill_manual(values = c("black", "red"))
ggsave("scripts/plotting/05a_Average_TW_Haplotypes_2021_2022.png", width = 16, height = 12)

### 5ai_Comparing_FT_Between_Haplotypes.R

Haplo_graph <- Single_2021_2022 %>% group_by(Generation, Haplotype) %>% summarise(Avg_FT = mean(FT_DAYS, na.rm = T)) %>% ungroup()
Haplo_graph <- Haplo_graph %>% group_by(Generation) %>% mutate(GenAvg = mean(Avg_FT, na.rm = T)) %>% ungroup()
Haplo_graph$`Avg FT > Pop & Gen Avg FT` <- ifelse(Haplo_graph$Avg_FT > Haplo_graph$GenAvg & Haplo_graph$Avg_FT > mean(Haplo_graph$Avg_FT, na.rm = T), TRUE, FALSE)
f_FT <- Haplo_graph %>% filter(`Avg FT > Pop & Gen Avg FT` == TRUE) %>% filter(!duplicated(Haplotype)) %>% select(Haplotype)

ggplot(Haplo_graph,aes(x= Haplotype, y = Avg_FT, fill = `Avg FT > Pop & Gen Avg FT`)) +
  geom_bar(stat = 'identity') +
  geom_hline(yintercept = mean(Haplo_graph$Avg_FT, na.rm = T), color = "green") +
  geom_hline(aes(yintercept = GenAvg), color = 'blue', linetype = 2) +
  labs(y = "Average Flowering Time (Days)",
       title = "Average Flowering Time between Haplotypes 2021-2022") +
  facet_wrap(~Generation, scales = "free_x") +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_fill_manual(values = c("black", "red"))
ggsave("scripts/plotting/05ai_Average_FT_Haplotypes_2021_2022.png", width = 16, height = 12)

### 5aii_Comparing_Fecundity_Between_Haplotypes.R

Haplo_graph <- Single_2021_2022 %>% group_by(Generation, Haplotype) %>% summarise(Avg_Fec = mean(FECUNDITY, na.rm = T)) %>% ungroup()
Haplo_graph <- Haplo_graph %>% filter(Avg_Fec != "NA")
Haplo_graph <- Haplo_graph %>% group_by(Generation) %>% mutate(GenAvg = mean(Avg_Fec, na.rm = T)) %>% ungroup()
Haplo_graph$`Avg Fec > Pop & Gen Avg Fec` <- ifelse(Haplo_graph$Avg_Fec > Haplo_graph$GenAvg & Haplo_graph$Avg_Fec > mean(Haplo_graph$Avg_Fec, na.rm = T), TRUE, FALSE)
f_Fe <- Haplo_graph %>% filter(`Avg Fec > Pop & Gen Avg Fec` == TRUE) %>% filter(!duplicated(Haplotype)) %>% select(Haplotype)

ggplot(Haplo_graph,aes(x= Haplotype, y = Avg_Fec, fill = `Avg Fec > Pop & Gen Avg Fec`)) +
  geom_bar(stat = 'identity') +
  geom_hline(yintercept = mean(Haplo_graph$Avg_Fec, na.rm = T), color = "green") +
  geom_hline(aes(yintercept = GenAvg), color = 'blue', linetype = 2) +
  labs(y = "Average Centered Fecundity",
       title = "Average Fecundity between Haplotypes 2021-2022") +
  facet_wrap(~Generation, scales = "free_x") +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_fill_manual(values = c("black", "red"))
ggsave("scripts/plotting/05aii_Average_Fec_Haplotypes_2021_2022.png", width = 16, height = 12)

### 5aiii_Comparing_Fitness_Between_Haplotypes.R

Haplo_graph <- Single_2021_2022 %>% group_by(Generation, Haplotype) %>% summarise(Avg_Fit = mean(ABS_FITNESS, na.rm = T)) %>% ungroup()
Haplo_graph <- Haplo_graph %>% filter(Avg_Fit != "NA")
Haplo_graph <- Haplo_graph %>% group_by(Generation) %>% mutate(GenAvg = mean(Avg_Fit, na.rm = T)) %>% ungroup()
Haplo_graph$`Avg Fit > Pop & Gen Avg Fit` <- ifelse(Haplo_graph$Avg_Fit > Haplo_graph$GenAvg & Haplo_graph$Avg_Fit > mean(Haplo_graph$Avg_Fit, na.rm = T), TRUE, FALSE)
f_Fi <- Haplo_graph %>% filter(`Avg Fit > Pop & Gen Avg Fit` == TRUE) %>% filter(!duplicated(Haplotype)) %>% select(Haplotype)
ui <- Reduce(intersect, list(f_AY, f_Fe, f_Fi))

ggplot(Haplo_graph,aes(x= Haplotype, y = Avg_Fit, fill = `Avg Fit > Pop & Gen Avg Fit`)) +
  geom_bar(stat = 'identity') +
  geom_hline(yintercept = mean(Haplo_graph$Avg_Fit, na.rm = TRUE), color = "green") +
  geom_hline(aes(yintercept = GenAvg), color = 'blue', linetype = 2) +
  labs(y = "Average Fitness",
       title = "Average Absolute Fitness between Haplotypes 2021-2022") +
  facet_wrap(~Generation, scales = "free_x") +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_fill_manual(values = c("black", "red"))
ggsave("scripts/plotting/05aiii_Average_Fit_Haplotypes_2021_2022.png", width = 16, height = 12)

### 5aiiii_Comparing_100SW_Between_Haplotypes.R

Haplo_graph <- Single_2021_2022 %>% group_by(Generation, Haplotype) %>% summarise(Avg_100SW = mean(`100_seed_weight`, na.rm = T)) %>% ungroup()
Haplo_graph <- Haplo_graph %>% filter(Avg_100SW != "NA")
Haplo_graph <- Haplo_graph %>% group_by(Generation) %>% mutate(GenAvg = mean(Avg_100SW, na.rm = T)) %>% ungroup()
Haplo_graph$`Avg Fit > Pop & Gen Avg Fit` <- ifelse(Haplo_graph$Avg_100SW > Haplo_graph$GenAvg & Haplo_graph$Avg_100SW > mean(Haplo_graph$Avg_100SW, na.rm = T), TRUE, FALSE)
f_Fd <- Haplo_graph %>% filter(`Avg Fit > Pop & Gen Avg Fit` == TRUE) %>% filter(!duplicated(Haplotype)) %>% select(Haplotype)
ui <- Reduce(intersect, list(f_AY, f_Fe, f_Fi, f_Fd))

ggplot(Haplo_graph,aes(x= Haplotype, y = Avg_100SW, fill = `Avg Fit > Pop & Gen Avg Fit`)) +
  geom_bar(stat = 'identity') +
  geom_hline(yintercept = mean(Haplo_graph$Avg_100SW, na.rm = TRUE), color = "green") +
  geom_hline(aes(yintercept = GenAvg), color = 'blue', linetype = 2) +
  labs(y = "Average 100 Seed Weight (grams)",
       title = "Average 100 Seed Weigh between Haplotypes 2021-2022") +
  facet_wrap(~Generation, scales = "free_x") +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_fill_manual(values = c("black", "red"))
ggsave("scripts/plotting/05aiii_Average_Fit_Haplotypes_2021_2022.png", width = 16, height = 12)


###########################################

# 2022-2023

Single_2022_2023 <- PHENO_FULL_AVERAGE %>% filter(Condition == 'single' & Exp_year == 2023)

### 5a_Comparing_TW_Between_Haplotypes.R

Single_2022_2023 <- Single_2022_2023 %>% filter(Haplotype != "NA")
Single_2022_2023$Haplotype <- as.factor(Single_2022_2023$Haplotype)
Haplo_graph <- Single_2022_2023 %>% group_by(Generation, Haplotype) %>% summarise(Avg_TW = mean(total_seed_mass_g, na.rm = T)) %>% ungroup()
Haplo_graph <- Haplo_graph %>% group_by(Generation) %>% mutate(GenAvg = mean(Avg_TW, na.rm = T)) %>% ungroup()
Haplo_graph$`Avg Yield > Pop & Gen Avg Yield` <- ifelse(Haplo_graph$Avg_TW > Haplo_graph$GenAvg & Haplo_graph$Avg_TW > mean(Haplo_graph$Avg_TW, na.rm = T), TRUE, FALSE)
f_AY <- Haplo_graph %>% filter(`Avg Yield > Pop & Gen Avg Yield` == TRUE) %>% filter(!duplicated(Haplotype)) %>% select(Haplotype)

ggplot(Haplo_graph,aes(x= Haplotype, y = Avg_TW, fill = `Avg Yield > Pop & Gen Avg Yield`)) +
  geom_bar(stat = 'identity') +
  geom_hline(yintercept = mean(Haplo_graph$Avg_TW, na.rm = T), color = "green") +
  geom_hline(aes(yintercept = GenAvg), color = 'blue', linetype = 2) +
  labs(y = "Average Total Seed Weight",
       title = "Average Total Seed Weight between Haplotypes 2022-2023") +
  scale_y_continuous(breaks = seq(0, 200, 25)) +
  facet_wrap(~Generation, scales = "free_x") +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_fill_manual(values = c("black", "red"))
ggsave("scripts/plotting/05a_Average_TW_Haplotypes_2021_2022.png", width = 16, height = 12)

### 5ai_Comparing_FT_Between_Haplotypes.R

Haplo_graph <- Single_2022_2023 %>% group_by(Generation, Haplotype) %>% summarise(Avg_FT = mean(FT_DAYS, na.rm = T)) %>% ungroup()
Haplo_graph <- Haplo_graph %>% group_by(Generation) %>% mutate(GenAvg = mean(Avg_FT, na.rm = T)) %>% ungroup()
Haplo_graph$`Avg FT > Pop & Gen Avg FT` <- ifelse(Haplo_graph$Avg_FT > Haplo_graph$GenAvg & Haplo_graph$Avg_FT > mean(Haplo_graph$Avg_FT, na.rm = T), TRUE, FALSE)
f_FT <- Haplo_graph %>% filter(`Avg FT > Pop & Gen Avg FT` == TRUE) %>% filter(!duplicated(Haplotype)) %>% select(Haplotype)

ggplot(Haplo_graph,aes(x= Haplotype, y = Avg_FT, fill = `Avg FT > Pop & Gen Avg FT`)) +
  geom_bar(stat = 'identity') +
  geom_hline(yintercept = mean(Haplo_graph$Avg_FT, na.rm = T), color = "green") +
  geom_hline(aes(yintercept = GenAvg), color = 'blue', linetype = 2) +
  labs(y = "Average Flowering Time (Days)",
       title = "Average Flowering Time between Haplotypes 2022-2023") +
  facet_wrap(~Generation, scales = "free_x") +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_fill_manual(values = c("black", "red"))
ggsave("scripts/plotting/05ai_Average_FT_Haplotypes_2021_2022.png", width = 16, height = 12)

### 5aii_Comparing_Fecundity_Between_Haplotypes.R

Haplo_graph <- Single_2022_2023 %>% group_by(Generation, Haplotype) %>% summarise(Avg_Fec = mean(FECUNDITY, na.rm = T)) %>% ungroup()
Haplo_graph <- Haplo_graph %>% filter(Avg_Fec != "NA")
Haplo_graph <- Haplo_graph %>% group_by(Generation) %>% mutate(GenAvg = mean(Avg_Fec, na.rm = T)) %>% ungroup()
Haplo_graph$`Avg Fec > Pop & Gen Avg Fec` <- ifelse(Haplo_graph$Avg_Fec > Haplo_graph$GenAvg & Haplo_graph$Avg_Fec > mean(Haplo_graph$Avg_Fec, na.rm = T), TRUE, FALSE)
f_Fe <- Haplo_graph %>% filter(`Avg Fec > Pop & Gen Avg Fec` == TRUE) %>% filter(!duplicated(Haplotype)) %>% select(Haplotype)

ggplot(Haplo_graph,aes(x= Haplotype, y = Avg_Fec, fill = `Avg Fec > Pop & Gen Avg Fec`)) +
  geom_bar(stat = 'identity') +
  geom_hline(yintercept = mean(Haplo_graph$Avg_Fec, na.rm = T), color = "green") +
  geom_hline(aes(yintercept = GenAvg), color = 'blue', linetype = 2) +
  labs(y = "Average Centered Fecundity",
       title = "Average Fecundity between Haplotypes 2022-2023") +
  facet_wrap(~Generation, scales = "free_x") +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_fill_manual(values = c("black", "red"))
ggsave("scripts/plotting/05aii_Average_Fec_Haplotypes_2021_2022.png", width = 16, height = 12)

### 5aiii_Comparing_Fitness_Between_Haplotypes.R

Haplo_graph <- Single_2022_2023 %>% group_by(Generation, Haplotype) %>% summarise(Avg_Fit = mean(ABS_FITNESS, na.rm = T)) %>% ungroup()
Haplo_graph <- Haplo_graph %>% filter(Avg_Fit != "NA")
Haplo_graph <- Haplo_graph %>% group_by(Generation) %>% mutate(GenAvg = mean(Avg_Fit, na.rm = T)) %>% ungroup()
Haplo_graph$`Avg Fit > Pop & Gen Avg Fit` <- ifelse(Haplo_graph$Avg_Fit > Haplo_graph$GenAvg & Haplo_graph$Avg_Fit > mean(Haplo_graph$Avg_Fit, na.rm = T), TRUE, FALSE)
f_Fi <- Haplo_graph %>% filter(`Avg Fit > Pop & Gen Avg Fit` == TRUE) %>% filter(!duplicated(Haplotype)) %>% select(Haplotype)
ui <- Reduce(intersect, list(f_AY, f_Fe, f_Fi))

ggplot(Haplo_graph,aes(x= Haplotype, y = Avg_Fit, fill = `Avg Fit > Pop & Gen Avg Fit`)) +
  geom_bar(stat = 'identity') +
  geom_hline(yintercept = mean(Haplo_graph$Avg_Fit, na.rm = TRUE), color = "green") +
  geom_hline(aes(yintercept = GenAvg), color = 'blue', linetype = 2) +
  labs(y = "Average Fitness",
       title = "Average Absolute Fitness between Haplotypes 2022-2023") +
  facet_wrap(~Generation, scales = "free_x") +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_fill_manual(values = c("black", "red"))
ggsave("scripts/plotting/05aiii_Average_Fit_Haplotypes_2021_2022.png", width = 16, height = 12)

### 5aiiii_Comparing_100SW_Between_Haplotypes.R

Haplo_graph <- Single_2022_2023 %>% group_by(Generation, Haplotype) %>% summarise(Avg_100SW = mean(`100_seed_weight`, na.rm = T)) %>% ungroup()
Haplo_graph <- Haplo_graph %>% filter(Avg_100SW != "NA")
Haplo_graph <- Haplo_graph %>% group_by(Generation) %>% mutate(GenAvg = mean(Avg_100SW, na.rm = T)) %>% ungroup()
Haplo_graph$`Avg Fit > Pop & Gen Avg Fit` <- ifelse(Haplo_graph$Avg_100SW > Haplo_graph$GenAvg & Haplo_graph$Avg_100SW > mean(Haplo_graph$Avg_100SW, na.rm = T), TRUE, FALSE)
f_Fd <- Haplo_graph %>% filter(`Avg Fit > Pop & Gen Avg Fit` == TRUE) %>% filter(!duplicated(Haplotype)) %>% select(Haplotype)
ui <- Reduce(intersect, list(f_AY, f_Fe, f_Fi, f_Fd))

ggplot(Haplo_graph,aes(x= Haplotype, y = Avg_100SW, fill = `Avg Fit > Pop & Gen Avg Fit`)) +
  geom_bar(stat = 'identity') +
  geom_hline(yintercept = mean(Haplo_graph$Avg_100SW, na.rm = TRUE), color = "green") +
  geom_hline(aes(yintercept = GenAvg), color = 'blue', linetype = 2) +
  labs(y = "Average 100 Seed Weight (grams)",
       title = "Average 100 Seed Weigh between Haplotypes 2022-2023") +
  facet_wrap(~Generation, scales = "free_x") +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_fill_manual(values = c("black", "red"))
ggsave("scripts/plotting/05aiii_Average_Fit_Haplotypes_2021_2022.png", width = 16, height = 12)
