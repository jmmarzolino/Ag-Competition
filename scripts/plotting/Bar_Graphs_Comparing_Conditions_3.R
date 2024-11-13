#!/usr/bin/env Rscript
#SBATCH --mem=30G
#SBATCH --time=02:00:00
#SBATCH --output=/rhome/jmarz001/bigdata/Ag-Competition/scripts/Bar_Graphs_Comparing_Conditions.stdout
#SBATCH -p koeniglab


library(tidyverse)
library(ggpubr)
library(data.table)

### 
setwd("")
df <- fread("FITNESS.tsv")
df$Generation <- as.numeric(df$Generation)
df$Genotype <- as.factor(df$Genotype)

### 3a_Bar_Graph_Avg_TW_Both_Conditions.R

ggplot(df, aes(Genotype, TOTAL_MASS, color = Condition, fill = Condition)) +
  geom_bar(stat = "identity", position = position_dodge(), alpha = .5, width = .5) +
  labs(y = "Total Seed Weight (g)",
       title = "Comparing Average Total Weight Between Conditions ") +
  facet_wrap(~Generation, scales = "free_x") +
  theme(axis.text.x = element_text(size = 5)) +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_y_continuous(breaks = seq(0, 200, 10))
ggsave("scripts/plotting/03a_Bar_Graph_TW_Between_Conditions.png", width = 30, height = 10)


### 3ai_Bar_Graph_Avg_FT_Both_Conditions.R

ggplot(df, aes(Genotype, FT, color = Condition, fill = Condition)) +
  geom_bar(stat = "identity", position = position_dodge(), alpha = .3, width = .5) +
  labs(y = "Flowering Time (Days)",
       title = "Comparing Average Flowering Time Between Conditions ") +
  facet_wrap(~Generation, scales = "free_x") +
  theme(axis.text.x = element_text(size = 5)) +
  theme(axis.text.x = element_text(angle = 90))
ggsave("scripts/plotting/03ai_Bar_Graph_FT_Between_Conditions.png", width = 18, height = 10)

### 3aii_Bar_Graph_Avg_Fec_Both_Conditions.R

ggplot(df, aes(Genotype, FECUNDITY, color = Condition, fill = Condition)) +
  geom_bar(stat = "identity", position = position_dodge(), alpha = .3, width = .5) +
  labs(y = "Fecundity",
       title = "Comparing Average Fecundity Between Conditions ") +
  facet_wrap(~Generation, scales = "free_x") +
  theme(axis.text.x = element_text(size = 5)) +
  theme(axis.text.x = element_text(angle = 90))
ggsave("scripts/plotting/03aii_Bar_Graph_Fec_Between_Conditions.png", width = 18, height = 10)

### 3aiii_Bar_Graph_Avg_Fit_Both_Conditions.R

ggplot(df, aes(Genotype, ABS_FITNESS, color = Condition, fill = Condition)) +
  geom_bar(stat = "identity", position = position_dodge(), alpha = .3, width = .5) +
  labs(y = "Fitness",
       title = "Comparing Average Absolute Fitness Between Conditions ") +
  facet_wrap(~Generation, scales = "free_x") +
  theme(axis.text.x = element_text(size = 5)) +
  theme(axis.text.x = element_text(angle = 90))
ggsave("scripts/plotting/03aiii_Bar_Graph_Fit_Between_Conditions.png",  width = 18, height = 10)

### 3aiiii_Bar_Graph_Avg_100SW_Both_Conditions.R

ggplot(df, aes(Genotype, SEED_WEIGHT_100, color = Condition, fill = Condition)) +
  geom_bar(stat = "identity", position = position_dodge(), alpha = .3, width =.5) +
  labs(y = "Average 100 Seed Weight",
       title = "Comparing Average 100 Seed Weight Between Conditions ") +
  facet_wrap(~Generation, scales = "free_x") +
  theme(axis.text.x = element_text(size = 5)) +
  theme(axis.text.x = element_text(angle = 90))
ggsave("scripts/plotting/03aiii_Bar_Graph_100SW_Between_Conditions.png", width = 18, height = 10)


### 3b_Average_TW_Over_Generations.R

fa <- ggplot(df, aes(x = Generation, y = TOTAL_MASS, color = Condition, add = "reg.line")) +
  geom_jitter(alpha = .6) +
  geom_smooth(method = lm) +
  stat_regline_equation() +
  labs(x = "Generation",
       y = "Average Total Weight (grams)",
       title = "2021-2021")
ggsave("scripts/plotting/03c_Scatterplot_Avg_TW_by_Condition.png")

### 3bi_Average_FT_Over_Generations.R

fb <- ggplot(df, aes(x = Generation, y = FT, color = Condition, add = "reg.line")) +
  geom_jitter(alpha = .6) +
  geom_smooth(method = lm) +
  stat_regline_equation() +
  labs(x = "Generation",
       y = "Average Flowering Time (Days After Sowing)",
       title = "2021-2021")
ggsave("scripts/plotting/03ci_Scatterplot_Avg_FT_by_Condition.png")

### 3bii_Average_Fecundity_Over_Generations.R

fc <- ggplot(df, aes(x = Generation, y = FECUNDITY, color = Condition, add = "reg.line")) +
  geom_jitter(alpha = .6) +
  geom_smooth(method = lm) +
  geom_hline(yintercept = 0, color = "red") +
  stat_regline_equation() +
  labs(x = "Generation",
       y = "Average Centered Fecundity",
       title = "2021-2021")
ggsave("scripts/plotting/03cii_Scatterplot_Avg_Fec_by_Condition.png")

### 3biii_Average_Fitness_Over_Generations.R

fd <- ggplot(df, aes(x = Generation, y = ABS_FITNESS, color = Condition, add = "reg.line")) +
  geom_jitter(alpha = .6) +
  geom_smooth(method = lm) +
  geom_hline(yintercept = 0, color = "red") +
  stat_regline_equation() +
  labs(x = "Generation",
       y = "Average Centered Absolute Fitness",
       title = "2021-2021")
ggsave("scripts/plotting/03ciii_Scatterplot_Avg_Fit_by_Condition.png")

### 3biii_Average_100SW_OVer_Generations.R

fe <- ggplot(df, aes(Generation, SEED_WEIGHT_100, color = Condition, add = "reg.line")) +
  geom_jitter(alpha = .6) +
  geom_smooth(method = "lm") +
  stat_regline_equation() +
  labs(y = "Average 100 Seed Weight (grams)",
       title = "")







### 3giiii-_Combined_Mixed_v_Single_Scatterplots.R

y <- grid.arrange(fa, fb, fc, fd, fe, ff, fg, fh, fi, fj, top = "Evolution of our Four Measured Phenotypes Between Conditions", nrow = 2, ncol = 5)
ggsave("scripts/plotting/03giiii_Combined_Mixed_v_Single_Scatterplots.png", y, width = 30, height = 18)









### 3ci_T_test_Mixed_vs_Single_TW

t.test(TOTAL_MASS ~ Condition, df)

### 3cii_T_test_Mixed_vs_Single_Fitness

t.test(ABS_FITNESS ~ Condition, df)

### 3ciii_T_test_Mixed_vs_Single_Fecundity

t.test(FECUNDITY ~ Condition, df)

### 3ciiii_T_test_Mixed_vs_Single_FT

t.test(FT ~ Condition, df)

### 3ciiiii_T_test_Mixed_vs_Single_100SW

t.test(SEED_WEIGHT_100 ~ Condition, df)

### 3d_Seed_TW_Per_Genotype

ggplot(df, aes(x = reorder(Genotype, + TOTAL_MASS), TOTAL_MASS, fill = Condition, group = Generation)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  labs(y = "Average Total Seed Weight",
       title = "Average Total Weight of Genotype Across Generations") +
  facet_wrap(~Generation, scales = "free_x") +
  theme(axis.text.x = element_text(angle = 90))

### 3e_Intermediate_FT_reproductive_success.R

# FT vs. Fit
ggplot(df, aes(FT, ABS_FITNESS, add = "reg.line")) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ x + I(x^2), size = 1) +
  facet_wrap(~Generation, scales = "free_x")

# FT vs. Fec
ggplot(df, aes(FT, FECUNDITY)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ x + I(x^2), size = 1) +
  facet_wrap(~Generation, scales = "free_x")

# FT vs. TW
ggplot(df, aes(FT, TOTAL_MASS)) +
  geom_point()+
  geom_smooth(method = "lm", formula = y ~ x + I(x^2), size = 1) +
  facet_wrap(~Generation, scales = "free_x") +
  labs(y = "Average Total Weight (grams")

# FT vs. 100SW

ggplot(df, aes(FT, SEED_WEIGHT_100)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ x + I(x^2), size = 1) +
  facet_wrap(~Generation, scales = "free_x") +
  labs(y = "Average 100 SW")

# FT vs. Fit
ggplot(tmp, aes(FT, ABS_FITNESS, add = "reg.line")) +
  geom_point() +
  geom_smooth(method = "lm")

# FT vs. Fec
ggplot(tmp, aes(FT, FECUNDITY)) +
  geom_point() +
  geom_smooth(method = "lm")

# FT vs. Total Weight
ggplot(tmp, aes(FT, TOTAL_MASS)) +
  geom_point()+
  geom_smooth(method = "lm") +
  ylim(0,200)

# FT vs. 100 SW

ggplot(tmp, aes(FT, SEED_WEIGHT_100)) +
  geom_point()+
  geom_smooth(method = "lm")

# 

