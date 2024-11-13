#!/usr/bin/env Rscript

#SBATCH --ntasks=1
#SBATCH --mem=30G
#SBATCH --time=02:00:00
#SBATCH --output=/rhome/jmarz001/bigdata/Ag-Competition/scripts/CCII_Ag_Comp_Graphs_3.stdout
#SBATCH -p koeniglab


library(tidyverse)
library(ggpubr)
library(ggplot2)

setwd("/bigdata/koeniglab/jmarz001/Ag-Competition/data")
fitness_df <- read_delim("FITNESS.tsv")



### 3a_bar_graph_avg_tw

ggplot(fitness_df, aes(Genotype, TOTAL_MASS, fill = Condition)) +
  geom_bar(stat = "identity", position = position_dodge(), width = .5) +
  labs(y = "Average Total Seed Weight (grams)",
       title = "Comparing Average Total Seed Weight Between Individual Genotype") +
  facet_wrap(~Generation, scales = "free_x") +
  theme(axis.text.x = element_text(size = 5)) +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_y_continuous(breaks = seq(0, 200, 10)) +
  theme_bw()
ggsave("/bigdata/koeniglab/jmarz001/Ag-Competition/results/03a_bar_graphs_avg_tw.png", width = 18, height = 10)

### 3a_bar_graph_avg_ft

ggplot(fitness_df, aes(Genotype, FT, fill = Condition)) +
  geom_bar(stat = "identity", position = position_dodge(), width = .5) +
  labs(y = "Flowering Time (Days)",
       title = "Comparing Flowering Time Between Individual Genotype") +
  facet_wrap(~Generation, scales = "free_x") +
  theme(axis.text.x = element_text(size = 5)) +
  theme(axis.text.x = element_text(angle = 90)) +
  theme_bw()
ggsave("/bigdata/koeniglab/jmarz001/Ag-Competition/results/03a_bar_graphs_avg_ft.png", width = 18, height = 10)

### 3a_bar_graphs_avg_fec

ggplot(fitness_df, aes(Genotype, FECUNDITY, fill = Condition)) +
  geom_bar(stat = "identity", position = position_dodge(), width = .5) +
  labs(y = "Fecundity",
       title = "Comparing Fecundity Between Individual Genotype") +
  facet_wrap(~Generation, scales = "free_x") +
  theme(axis.text.x = element_text(size = 5)) +
  theme(axis.text.x = element_text(angle = 90)) +
  theme_bw()
ggsave("/bigdata/koeniglab/jmarz001/Ag-Competition/results/3a_bar_graphs_avg_fec.png"width = 18, height = 10)
### 3a_bar_graph_avg_fit

ggplot(fitness_df, aes(Genotype, FITNESS, fill = Condition)) +
  geom_bar(stat = "identity", position = position_dodge(), width = .5) +
  labs(y = "Fitness",
       title = "Comparing Fitness Between Individual Genotype") +
  facet_wrap(~Generation, scales = "free_x") +
  theme(axis.text.x = element_text(size = 5)) +
  theme(axis.text.x = element_text(angle = 90)) +
  theme_bw()
ggsave("/bigdata/koeniglab/jmarz001/Ag-Competition/results/03_bar_graphs_avg_fit.png", width = 18, height = 10)

### 3a_bar_graph_avg_100sw





## Expected yields
### 3b_Exp_Yield_Per_Plant_by_Condition

ggplot(fitness_df, aes(Generation, Exp_TW_Per_Plant, add = "reg.line")) +
  geom_jitter() +
  geom_smooth(method = lm) +
  stat_regline_equation() +
  labs(y = "Expected Yield Per Plant (g)",
       title = "Generational Change in Expected Yield Per Plant")
ggsave("/bigdata/koeniglab/jmarz001/Ag-Competition/results/03b_Exp_Per_Plant_Yield_by_Condition.png")

### 3bi_Exp_Fec_Per_Plant_by_Condition

ggplot(Average_Haplo_rep, aes(Generation, Exp_Fec_Per_Plant, add = "reg.line")) +
  geom_jitter() +
  geom_smooth(method = lm) +
  stat_regline_equation() +
  labs(y = "Expected Fecundity Per Plant",
       title = "Generational Change in Expected Fecundity Per Plant")
ggsave("/bigdata/koeniglab/jmarz001/Ag-Competition/results/03bi_Exp_Per_Plant_Fec_by_Condition.png")

### 3bii_Exp_Fit_by_Condition

ggplot(Average_Haplo_rep, aes(Generation, Exp_Fit_Per_Plant, add = "reg.line")) +
  geom_jitter() +
  geom_smooth(method = lm) +
  stat_regline_equation() +
  labs(y = "Expected Fitness Per Plant",
       title = "Generational Change in Expected Fitness Per Plant")
ggsave("/bigdata/koeniglab/jmarz001/Ag-Competition/results/03bii_Exp_Per_Plant_Fit_by_Condition.png")

### 3c_Average_TW_Over_Generations

fa <- ggplot(Average_Haplo_rep, aes(x = Generation, y = Centered_TW, add = "reg.line")) +
  geom_jitter() +
  geom_smooth(method = lm) +
  geom_hline(yintercept = 0, color = "red") +
  stat_regline_equation() +
  labs(x = "Generation",
       y = "Average Total Seed Weight (grams)")
ggsave("/bigdata/koeniglab/jmarz001/Ag-Competition/results/03c_Scatterplot_Avg_TW_by_Condition.png")

### 3ci_Average_FT_Over_Generations

fb <- ggplot(Average_Haplo_rep, aes(x = Generation, y = Centered_FT, add = "reg.line")) +
  geom_jitter() +
  geom_smooth(method = lm) +
  geom_hline(yintercept = 0, color = "red") +
  stat_regline_equation() +
  labs(x = "Generation",
       y = "Average Flowering Time (Days After Sowing)") +
  scale_y_continuous(breaks = seq(0, 140, 10))
ggsave("scripts/plotting/03ci_Scatterplot_Avg_FT_by_Condition.png")

### 3cii_Average_Fecundity_Over_Generations

fc <- ggplot(Average_Haplo_rep, aes(x = Generation, y = Centered_Fec, add = "reg.line")) +
  geom_jitter() +
  geom_smooth(method = lm) +
  geom_hline(yintercept = 0, color = "red") +
  stat_regline_equation() +
  labs(x = "Generation",
       y = "Average Fecundity")
ggsave("scripts/plotting/03cii_Scatterplot_Avg_Fec_by_Condition.png")

### 3ciii_Average_Fitness_Over_Generations

fd <- ggplot(Average_Haplo_rep, aes(x = Generation, y = Centered_Fit, add = "reg.line")) +
  geom_jitter() +
  geom_smooth(method = lm) +
  geom_hline(yintercept = 0, color = "red") +
  stat_regline_equation() +
  labs(x = "Generation",
       y = "Average Fitness")
ggsave("/bigdata/koeniglab/jmarz001/Ag-Competition/results/03ciii_Scatterplot_Avg_Fit_by_Condition.png")

### 3ciiii_Combined_Mixed_v_Single_Scatterplots

g <- grid.arrange(fa, fb, fc, fd, top = "Evolution of our Four Measured Phenotypes Between Condition 2022-2023 Season")
ggsave("/bigdata/koeniglab/jmarz001/Ag-Competition/results/03ciiii_Combined_Mixed_v_Single_Scatterplots.png", g, width = 14, height = 10)

### 3cc_T_test_Mixed_vs_Single_Yield

t.test(TOTAL_WEIGHT ~ Condition, Average_Haplo_rep)

### 3cci_T_test_Mixed_vs_Single_Fitness

t.test(Fitness ~ Condition, Average_Haplo_rep)

### 3ccii_T_test_Mixed_vs_Single_Fecundity

t.test(Fecundity ~ Condition, Average_Haplo_rep)

### 3d_Seed_TW_Per_Genotype

ggplot(Average_Haplo_rep, aes(x = reorder(Genotype, +TOTAL_WEIGHT), TOTAL_WEIGHT, fill = Condition, group = Generation)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  labs(y = "Average Total Seed Weight",
       title = "Average Total Weight of Genotype Across Generations") +
  facet_wrap(~Generation, scales = "free_x") +
  theme(axis.text.x = element_text(angle = 90))

### 3e_Intermediate_FT_reproductive_success

# FT vs. Fit
ggplot(Rep_Single, aes(FT, Fitness, add = "reg.line")) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ x + I(x^2), size = 1) +
  facet_wrap(~Generation, scales = "free_x")
ggsave("/bigdata/koeniglab/jmarz001/Ag-Competition/results/03e_Int_FT_vs_Fit.png")

# FT vs. Fec
ggplot(Rep_Single, aes(FT, Fecundity)) +
  geom_point() +
  geom_smooth() +
  facet_wrap(~Generation, scales = "free_x")
ggsave("/bigdata/koeniglab/jmarz001/Ag-Competition/results/03e_Int_FT_vs_fec.png")

# FT vs. TW
ggplot(Rep_Single, aes(FT, TOTAL_WEIGHT)) +
  geom_point()+
  geom_smooth() +
  labs(y = "Total Seed Weight (g)") +
  facet_wrap(~Generation, scales = "free_x")
ggsave("/bigdata/koeniglab/jmarz001/Ag-Competition/results/03e_Int_FT_vs_TW.png")






#### THESE ARE ADDITIONAL GRAPHS FOR MANUSCRIPT (Scatterplot Graphs comparing between the two conditions, then one combined graph)

### Avg_Total_Seed_Weight_Over_Time

fe <- ggplot(Average_Haplo_rep, aes(Generation, TOTAL_WEIGHT, add = "reg.line")) +
  geom_jitter() +
  geom_smooth(method = lm) +
  stat_regline_equation(label.y = 200) +
  labs(x = "Generation",
       y = "Average Total Seed Weight (grams)") +
  facet_wrap(~Condition)

ff <- ggplot(Average_Haplo_rep, aes(Generation, TOTAL_WEIGHT, col = Condition)) +
  geom_jitter() +
  geom_smooth(method = lm) +
  labs(x = "Generation",
       y = "Average Total Seed Weight (grams)") 

Compare_TW <- grid.arrange(fe,ff)
ggsave("scripts/plotting/Extra_3_Average_Yield_Comparisons.png", Compare_TW)

### Avg_FT_Over_Generations

fe <- ggplot(Average_Haplo_rep, aes(Generation, FT, add = "reg.line")) +
  geom_jitter() +
  geom_smooth(method = lm) +
  stat_regline_equation(label.y = 130) +
  labs(x = "Generation",
       y = "Average Flowering Time (Days)") +
  facet_wrap(~Condition)

ff <- ggplot(Average_Haplo_rep, aes(Generation, FT, col = Condition)) +
  geom_jitter() +
  geom_smooth(method = lm) +
  labs(x = "Generation",
       y = "Average Total Seed Weight (grams)") 

Compare_TW <- grid.arrange(fe,ff)
ggsave("scripts/plotting/Extra_3_Average_FT_Comparisons.png", Compare_TW)

### Avg_Fec_Over_Generations

fe <- ggplot(Average_Haplo_rep, aes(Generation, Fecundity, add = "reg.line")) +
  geom_jitter() +
  geom_smooth(method = lm) +
  stat_regline_equation(label.y = 3500) +
  labs(x = "Generation",
       y = "Average Fecundity") +
  facet_wrap(~Condition)

ff <- ggplot(Average_Haplo_rep, aes(Generation, Fecundity, col = Condition)) +
  geom_jitter() +
  geom_smooth(method = lm) +
  labs(x = "Generation",
       y = "Average Fecundity") 

Compare_TW <- grid.arrange(fe,ff)
ggsave("scripts/plotting/Extra_3_Average_Fecundity_Comparisons.png", Compare_TW)

### Avg_Fit_Over_Generations

fe <- ggplot(Average_Haplo_rep, aes(Generation, Fitness, add = "reg.line")) +
  geom_jitter() +
  geom_smooth(method = lm) +
  stat_regline_equation(label.y = 31000) +
  labs(x = "Generation",
       y = "Average Fitness") +
  facet_wrap(~Condition)

ff <- ggplot(Average_Haplo_rep, aes(Generation, Fitness, col = Condition)) +
  geom_jitter() +
  geom_smooth(method = lm) +
  labs(x = "Generation",
       y = "Average Fitness") 

Compare_TW <- grid.arrange(fe,ff)
ggsave("/bigdata/koeniglab/jmarz001/Ag-Competition/results/Extra_3_Average_Fitness_Comparisons.png", Compare_TW)
