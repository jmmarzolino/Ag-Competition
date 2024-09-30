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
library(car)
library(gridExtra)

### 2021-2022

tmp2021 <- PHENO_FULL_AVERAGE %>% filter(Exp_year == 2022)
tmp2021$Generation <- as.numeric(tmp2021$Generation)
tmp2021$Genotypes <- as.factor(tmp2021$Genotypes)

### 3a_Bar_Graph_Avg_TW_Both_Conditions.R

ggplot(tmp2021, aes(Genotypes, total_seed_mass_g, color = Condition, fill = Condition)) +
  geom_bar(stat = 'identity', position = position_dodge(), alpha = .5, width = .5) +
  labs(y = "Total Seed Weight (g)",
       title = "Comparing Average Total Weight Between Conditions 2021-2022") +
  facet_wrap(~Generation, scales = "free_x") +
  theme(axis.text.x = element_text(size = 5)) +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_y_continuous(breaks = seq(0, 200, 10))
ggsave("scripts/plotting/03a_Bar_Graph_TW_Between_Conditions_2021_2022.png", width = 30, height = 10)

### 3ai_Bar_Graph_Avg_FT_Both_Conditions.R

ggplot(tmp2021, aes(Genotypes, FT_DAYS, color = Condition, fill = Condition)) +
  geom_bar(stat = 'identity', position = position_dodge(), alpha = .3, width = .5) +
  labs(y = "Flowering Time (Days)",
       title = "Comparing Average Flowering Time Between Conditions 2021-2022") +
  facet_wrap(~Generation, scales = "free_x") +
  theme(axis.text.x = element_text(size = 5)) +
  theme(axis.text.x = element_text(angle = 90))
ggsave("scripts/plotting/03ai_Bar_Graph_FT_Between_Conditions_2021_2022.png", width = 18, height = 10)

### 3aii_Bar_Graph_Avg_Fec_Both_Conditions.R

ggplot(tmp2021, aes(Genotypes, FECUNDITY, color = Condition, fill = Condition)) +
  geom_bar(stat = 'identity', position = position_dodge(), alpha = .3, width = .5) +
  labs(y = "Fecundity",
       title = "Comparing Average Fecundity Between Conditions 2021-2022") +
  facet_wrap(~Generation, scales = "free_x") +
  theme(axis.text.x = element_text(size = 5)) +
  theme(axis.text.x = element_text(angle = 90))
ggsave("scripts/plotting/03aii_Bar_Graph_Fec_Between_Conditions_2021_2022.png", width = 18, height = 10)

### 3aiii_Bar_Graph_Avg_Fit_Both_Conditions.R

ggplot(tmp2021, aes(Genotypes, ABS_FITNESS, color = Condition, fill = Condition)) +
  geom_bar(stat = 'identity', position = position_dodge(), alpha = .3, width = .5) +
  labs(y = "Fitness",
       title = "Comparing Average Absolute Fitness Between Conditions 2021-2022") +
  facet_wrap(~Generation, scales = "free_x") +
  theme(axis.text.x = element_text(size = 5)) +
  theme(axis.text.x = element_text(angle = 90))
ggsave("scripts/plotting/03aiii_Bar_Graph_Fit_Between_Conditions_2021_2022.png",  width = 18, height = 10)

### 3aiiii_Bar_Graph_Avg_100SW_Both_Conditions.R

ggplot(tmp2021, aes(Genotypes, `100_seed_weight`, color = Condition, fill = Condition)) +
  geom_bar(stat = 'identity', position = position_dodge(), alpha = .3, width =.5) +
  labs(y = "Average 100 Seed Weight",
       title = "Comparing Average 100 Seed Weight Between Conditions 2021-2022") +
  facet_wrap(~Generation, scales = "free_x") +
  theme(axis.text.x = element_text(size = 5)) +
  theme(axis.text.x = element_text(angle = 90))
ggsave("scripts/plotting/03aiii_Bar_Graph_100SW_Between_Conditions_2021_2022.png", width = 18, height = 10)


### 3b_Average_TW_Over_Generations.R

fa <- ggplot(tmp2021, aes(x = Generation, y = total_seed_mass_g, color = Condition, add = "reg.line")) +
  geom_jitter(alpha = .6) +
  geom_smooth(method = lm) +
  stat_regline_equation() +
  labs(x = "Generation",
       y = "Average Total Weight (grams)",
       title = "2021-2021")
ggsave("scripts/plotting/03c_Scatterplot_Avg_TW_by_Condition_2021_2022.png")

### 3bi_Average_FT_Over_Generations.R

fb <- ggplot(tmp2021, aes(x = Generation, y = FT_DAYS, color = Condition, add = "reg.line")) +
  geom_jitter(alpha = .6) +
  geom_smooth(method = lm) +
  stat_regline_equation() +
  labs(x = "Generation",
       y = "Average Flowering Time (Days After Sowing)",
       title = "2021-2021")
ggsave("scripts/plotting/03ci_Scatterplot_Avg_FT_by_Condition_2021_2022.png")

### 3bii_Average_Fecundity_Over_Generations.R

fc <- ggplot(tmp2021, aes(x = Generation, y = FECUNDITY, color = Condition, add = "reg.line")) +
  geom_jitter(alpha = .6) +
  geom_smooth(method = lm) +
  geom_hline(yintercept = 0, color = 'red') +
  stat_regline_equation() +
  labs(x = "Generation",
       y = "Average Centered Fecundity",
       title = "2021-2021")
ggsave("scripts/plotting/03cii_Scatterplot_Avg_Fec_by_Condition_2021_2022.png")

### 3biii_Average_Fitness_Over_Generations.R

fd <- ggplot(tmp2021, aes(x = Generation, y = ABS_FITNESS, color = Condition, add = "reg.line")) +
  geom_jitter(alpha = .6) +
  geom_smooth(method = lm) +
  geom_hline(yintercept = 0, color = 'red') +
  stat_regline_equation() +
  labs(x = "Generation",
       y = "Average Centered Absolute Fitness",
       title = "2021-2021")
ggsave("scripts/plotting/03ciii_Scatterplot_Avg_Fit_by_Condition_2021_2022.png")

### 3biii_Average_100SW_OVer_Generations.R

fe <- ggplot(tmp2021, aes(Generation, `100_seed_weight`, color = Condition, add = 'reg.line')) +
  geom_jitter(alpha = .6) +
  geom_smooth(method = 'lm') +
  stat_regline_equation() +
  labs(y = "Average 100 Seed Weight (grams)",
       title = "2021-2022")

### 3ci_T_test_Mixed_vs_Single_TW

t.test(total_seed_mass_g ~ Condition, tmp2021)

### 3cii_T_test_Mixed_vs_Single_Fitness

t.test(ABS_FITNESS ~ Condition, tmp2021)

### 3ciii_T_test_Mixed_vs_Single_Fecundity

t.test(FECUNDITY ~ Condition, tmp2021)

### 3ciiii_T_test_Mixed_vs_Single_FT

t.test(FT_DAYS ~ Condition, tmp2021)

### 3ciiiii_T_test_Mixed_vs_Single_100SW

t.test(`100_seed_weight` ~ Condition, tmp2021)

### 3d_Seed_TW_Per_Genotype

ggplot(tmp2021, aes(x = reorder(Genotypes, + total_seed_mass_g), total_seed_mass_g, fill = Condition, group = Generation)) +
  geom_bar(stat = 'identity', position = position_dodge()) +
  labs(y = "Average Total Seed Weight",
       title = "Average Total Weight of Genotypes Across Generations") +
  facet_wrap(~Generation, scales = "free_x") +
  theme(axis.text.x = element_text(angle = 90))

### 3e_Intermediate_FT_reproductive_success.R

# FT vs. Fit
ggplot(tmp2021, aes(FT_DAYS, ABS_FITNESS, add = "reg.line")) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ x + I(x^2), size = 1) +
  facet_wrap(~Generation, scales = "free_x")

# FT vs. Fec
ggplot(tmp2021, aes(FT_DAYS, FECUNDITY)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ x + I(x^2), size = 1) +
  facet_wrap(~Generation, scales = "free_x")

# FT vs. TW
ggplot(tmp2021, aes(FT_DAYS, total_seed_mass_g)) +
  geom_point()+
  geom_smooth(method = "lm", formula = y ~ x + I(x^2), size = 1) +
  facet_wrap(~Generation, scales = "free_x") +
  labs(y = "Average Total Weight (grams")

# FT vs. 100SW

ggplot(tmp2021, aes(FT_DAYS, `100_seed_weight`)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ x + I(x^2), size = 1) +
  facet_wrap(~Generation, scales = "free_x") +
  labs(y = "Average 100 SW")


########################################################################

# 2022-2023


tmp2023 <- PHENO_FULL_AVERAGE %>% filter(Exp_year == 2023)
tmp2023$Generation <- as.numeric(tmp2023$Generation)
tmp2023$Genotypes <- as.factor(tmp2023$Genotypes)

### 3f_Bar_Graph_Avg_TW_Both_Conditions.R

ggplot(tmp2023, aes(Genotypes, total_seed_mass_g, color = Condition, fill = Condition)) +
  geom_bar(stat = 'identity', position = position_dodge(), alpha = .5, width = .5) +
  labs(y = "Total Seed Weight (g)",
       title = "Comparing Average Total Weight Between Conditions 2022-2023") +
  facet_wrap(~Generation, scales = "free_x") +
  theme(axis.text.x = element_text(size = 5)) +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_y_continuous(breaks = seq(0, 200, 10))
ggsave("scripts/plotting/03a_Bar_Graph_TW_Between_Conditions_2022_2023.png", width = 30, height = 10)

### 3fi_Bar_Graph_Avg_FT_Both_Conditions.R

ggplot(tmp2023, aes(Genotypes, FT_DAYS, color = Condition, fill = Condition)) +
  geom_bar(stat = 'identity', position = position_dodge(), alpha = .3, width = .5) +
  labs(y = "Flowering Time (Days)",
       title = "Comparing Average Flowering Time Between Conditions 2022-2023") +
  facet_wrap(~Generation, scales = "free_x") +
  theme(axis.text.x = element_text(size = 5)) +
  theme(axis.text.x = element_text(angle = 90))
ggsave("scripts/plotting/03ai_Bar_Graph_FT_Between_Conditions_2022_2023.png", width = 18, height = 10)

### 3fii_Bar_Graph_Avg_Fec_Both_Conditions.R

ggplot(tmp2023, aes(Genotypes, FECUNDITY, color = Condition, fill = Condition)) +
  geom_bar(stat = 'identity', position = position_dodge(), alpha = .3, width = .5) +
  labs(y = "Fecundity",
       title = "Comparing Average Fecundity Between Conditions 2022-2023") +
  facet_wrap(~Generation, scales = "free_x") +
  theme(axis.text.x = element_text(size = 5)) +
  theme(axis.text.x = element_text(angle = 90))
ggsave("scripts/plotting/03aii_Bar_Graph_Fec_Between_Conditions_2022_2023.png", width = 18, height = 10)

### 3fiii_Bar_Graph_Avg_Fit_Both_Conditions.R

ggplot(tmp2023, aes(Genotypes, ABS_FITNESS, color = Condition, fill = Condition)) +
  geom_bar(stat = 'identity', position = position_dodge(), alpha = .3, width = .5) +
  labs(y = "Fitness",
       title = "Comparing Average Absolute Fitness Between Conditions 2022-2023") +
  facet_wrap(~Generation, scales = "free_x") +
  theme(axis.text.x = element_text(size = 5)) +
  theme(axis.text.x = element_text(angle = 90))
ggsave("scripts/plotting/03aiii_Bar_Graph_Fit_Between_Conditions_2022_2023.png",  width = 18, height = 10)

### 3fiiii_Bar_Graph_Avg_100SW_Both_Conditions.R

ggplot(tmp2023, aes(Genotypes, `100_seed_weight`, color = Condition, fill = Condition)) +
  geom_bar(stat = 'identity', position = position_dodge(), alpha = .3, width =.5) +
  labs(y = "Average 100 Seed Weight",
       title = "Comparing Average 100 Seed Weight Between Conditions 2022-2023") +
  facet_wrap(~Generation, scales = "free_x") +
  theme(axis.text.x = element_text(size = 5)) +
  theme(axis.text.x = element_text(angle = 90))
ggsave("scripts/plotting/03aiii_Bar_Graph_100SW_Between_Conditions_2022_2023.png", width = 18, height = 10)

### 3g_Average_TW_Over_Generations.R

ff <- ggplot(tmp2023, aes(x = Generation, y = total_seed_mass_g, color = Condition, add = "reg.line")) +
  geom_jitter(alpha = .6) +
  geom_smooth(method = lm) +
  stat_regline_equation() +
  labs(x = "Generation",
       y = "Average Total Weight (grams)",
       title = "2022-2023")
ggsave("scripts/plotting/03c_Scatterplot_Avg_TW_by_Condition_2022_2023.png")

### 3gi_Average_FT_Over_Generations.R

fg <- ggplot(tmp2023, aes(x = Generation, y = FT_DAYS, color = Condition, add = "reg.line")) +
  geom_jitter(alpha = .6) +
  geom_smooth(method = lm) +
  stat_regline_equation() +
  labs(x = "Generation",
       y = "Average Flowering Time (Days After Sowing)",
       title = "2022-2023")
ggsave("scripts/plotting/03ci_Scatterplot_Avg_FT_by_Condition_2022_2023.png")

### 3gii_Average_Fecundity_Over_Generations.R

fh <- ggplot(tmp2023, aes(x = Generation, y = FECUNDITY, color = Condition, add = "reg.line")) +
  geom_jitter() +
  geom_smooth(method = lm) +
  geom_hline(yintercept = 0, color = 'red') +
  stat_regline_equation() +
  labs(x = "Generation",
       y = "Average Centered Fecundity",
       title = "2022-2023")
ggsave("scripts/plotting/03cii_Scatterplot_Avg_Fec_by_Condition_2022_2023.png")

### 3giii_Average_Fitness_Over_Generations.R

fi <- ggplot(tmp2023, aes(x = Generation, y = ABS_FITNESS, color = Condition, add = "reg.line")) +
  geom_jitter(alpha =.6) +
  geom_smooth(method = lm) +
  geom_hline(yintercept = 0, color = 'red') +
  stat_regline_equation() +
  labs(x = "Generation",
       y = "Average Centered Absolute Fitness",
       title = "2022-2023")
ggsave("scripts/plotting/03ciii_Scatterplot_Avg_Fit_by_Condition_2022_2023.png")

### 3giiii_Average_100SW_Over_Generations.R

fj <- ggplot(tmp2023, aes(Generation, `100_seed_weight`, color = Condition, add = "reg.line")) +
  geom_jitter(alpha =.6)+
  geom_smooth(method = "lm") +
  stat_regline_equation()+
  labs(y = "Average 100 Seed Weight (grams)",
       title = "2022-2023")

### 3giiii-_Combined_Mixed_v_Single_Scatterplots.R

y <- grid.arrange(fa, fb, fc, fd, fe, ff, fg, fh, fi, fj, top = "Evolution of our Four Measured Phenotypes Between Conditions", nrow = 2, ncol = 5)
ggsave("scripts/plotting/03giiii_Combined_Mixed_v_Single_Scatterplots.png", y, width = 30, height = 18)

### 3h_T_test_Mixed_vs_Single_TW

t.test(total_seed_mass_g ~ Condition, tmp2023)

### 3hi_T_test_Mixed_vs_Single_Fitness

t.test(ABS_FITNESS ~ Condition, tmp2023)

### 3hii_T_test_Mixed_vs_Single_Fecundity

t.test(FECUNDITY ~ Condition, tmp2023)

### 3hiii_T_test_Mixed_vs_Single_FT

t.test(FT_DAYS ~ Condition, tmp2023)

### 3hiiii_T_test_Mixed_vs_Single_100SW

t.test(`100_seed_weight` ~ Condition, tmp2023)

### 3k_Seed_TW_Per_Genotype

ggplot(Averaged_Full_2021_2022, aes(x = reorder(Genotypes, + total_seed_mass_g), total_seed_mass_g, fill = Condition, group = Generation)) +
  geom_bar(stat = 'identity', position = position_dodge()) +
  labs(y = "Average Total Seed Weight",
       title = "Average Total Weight of Genotypes Across Generations") +
  facet_wrap(~Generation, scales = "free_x") +
  theme(axis.text.x = element_text(angle = 90))

### 3l_Intermediate_FT_reproductive_success.R

# 2022-2023

# FT vs. Fit
ggplot(tmp2023, aes(FT_DAYS, ABS_FITNESS, add = "reg.line")) +
  geom_point() +
  geom_smooth(method = "lm")

# FT vs. Fec
ggplot(tmp2023, aes(FT_DAYS, FECUNDITY)) +
  geom_point() +
  geom_smooth(method = "lm")

# FT vs. Total Weight
ggplot(tmp2023, aes(FT_DAYS, total_seed_mass_g)) +
  geom_point()+
  geom_smooth(method = "lm") +
  ylim(0,200)

# FT vs. 100 SW

ggplot(tmp2023, aes(FT_DAYS, `100_seed_weight`)) +
  geom_point()+
  geom_smooth(method = "lm")

# 2021-2022

# FT vs. Fit

ggplot(tmp2021, aes(FT_DAYS, ABS_FITNESS, add = "reg.line")) +
  geom_point() +
  geom_smooth(method = "lm")

# FT vs. Fec
ggplot(tmp2021, aes(FT_DAYS, FECUNDITY)) +
  geom_point() +
  geom_smooth(method = "lm")


# FT vs. Total Weight
ggplot(tmp2021, aes(FT_DAYS, total_seed_mass_g)) +
  geom_point()+
  geom_smooth(method = "lm") +
  ylim(0, 200)

# FT vs. 100 SW

ggplot(tmp2021, aes(FT_DAYS, `100_seed_weight`)) +
  geom_point()+
  geom_smooth(method = "lm")
