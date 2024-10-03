library(tidyverse)
library(readr)
library(dplyr)
library(ggpubr)
library(ggplot2)
library(tidyr)
library(car)
library(gridExtra)

test <- Average_Haplo_rep %>% select(c("Genotype", "Generation", "Condition", "Brown Bag Weight", "Fitness", "Fecundity", "FT")) %>% group_by(Genotype) %>% 
  pivot_wider(names_from = "Condition", values_from = c("Brown Bag Weight", "Fecundity", "Fitness"))
test <- ifelse(test$)


test$new <- ifelse(test$`Brown Bag Weight_single` > test$`Brown Bag Weight_single`, 1,0)


### 3a_Bar_Graph_Avg_TW_Between_Genotype.R

ggplot(Averaged_Full_2021_2022, aes(Genotype, TOTAL_MASS, color = Condition, fill = Condition)) +
  geom_bar(stat = 'identity', position = position_dodge(), alpha = .5, width = .5) +
  labs(y = "Total Seed Weight (g)",
       title = "Comparing Average Total Seed Weight Between Individual Genotype") +
  facet_wrap(~Generation, scales = "free_x") +
  theme(axis.text.x = element_text(size = 5)) +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_y_continuous(breaks = seq(0, 200, 10))
ggsave("scripts/plotting/03a_Bar_Graph_TW_Between_Conditions_2021_2022.png", width = 18, height = 10)

### 3ai_Bar_Graph_Avg_FT_Between_Genotype.R

ggplot(Averaged_Full_2021_2022, aes(Genotype, FT, color = Condition, fill = Condition)) +
  geom_bar(stat = 'identity', position = position_dodge(), alpha = .3, width = .5) +
  labs(y = "Flowering Time (Days)",
       title = "Comparing Flowering Time Between Individual Genotype") +
  facet_wrap(~Generation, scales = "free_x") +
  theme(axis.text.x = element_text(size = 5)) +
  theme(axis.text.x = element_text(angle = 90))
ggsave("scripts/plotting/03ai_Bar_Graph_FT_Between_Conditions_2021_2022.png", width = 18, height = 10)

### 3aii_Bar_Graph_Avg_Fec_Between_Genotype.R

ggplot(Averaged_Full_2021_2022, aes(Genotype, Fecundity, color = Condition, fill = Condition)) +
  geom_bar(stat = 'identity', position = position_dodge(), alpha = .3, width = .5) +
  labs(y = "Fecundity",
       title = "Comparing Fecundity Between Individual Genotype") +
  facet_wrap(~Generation, scales = "free_x") +
  theme(axis.text.x = element_text(size = 5)) +
  theme(axis.text.x = element_text(angle = 90))
ggsave("scripts/plotting/03aii_Bar_Graph_Fec_Between_Conditions_2021_2022.png", width = 18, height = 10)

### 3aiii_Bar_Graph_Avg_Fit_Between_Genotype.R

ggplot(Averaged_Full_2021_2022, aes(Genotype, Fitness, color = Condition, fill = Condition)) +
  geom_bar(stat = 'identity', position = position_dodge(), alpha = .3, width = .5) +
  labs(y = "Fitness",
       title = "Comparing Fitness Between Individual Genotype") +
  facet_wrap(~Generation, scales = "free_x") +
  theme(axis.text.x = element_text(size = 5)) +
  theme(axis.text.x = element_text(angle = 90))
ggsave("scripts/plotting/03aiii_Bar_Graph_Fit_Between_Conditions_2021_2022.png",  width = 18, height = 10)

### 3b_Exp_TW_Per_Plant_by_Condition.R

ggplot(Averaged_Full_2021_2022, aes(Generation, Exp_TW_Per_Plant, color = Condition, add = "reg.line")) +
  geom_jitter() +
  geom_smooth(method = lm) +
  stat_regline_equation() +
  labs(y = "Expected Total Seed Weight Per Plant (g)",
       title = "Generational Change in Expected Total Seed Weight Per Plant")
ggsave("scripts/plotting/03b_Exp_Per_Plant_Yield_by_Condition_2021_2022.png")

### 3bi_Exp_Fec_Per_Plant_by_Condition.R

ggplot(Averaged_Full_2021_2022, aes(Generation, Exp_Fec_Per_Plant, color = Condition, add = "reg.line")) +
  geom_jitter() +
  geom_smooth(method = lm) +
  stat_regline_equation() +
  labs(y = "Expected Fecundity Per Plant",
       title = "Generational Change in Expected Fecundity Per Plant")
ggsave("scripts/plotting/03bi_Exp_Per_Plant_Fec_by_Condition_2021_2022.png")

### 3bii_Exp_Fit_by_Condition.R

ggplot(Averaged_Full_2021_2022, aes(Generation, Exp_Fit_Per_Plant, color = Condition, add = "reg.line")) +
  geom_jitter() +
  geom_smooth(method = lm) +
  stat_regline_equation() +
  labs(y = "Expected Fitness Per Plant",
       title = "Generational Change in Expected Fitness Per Plant")
ggsave("scripts/plotting/03bii_Exp_Per_Plant_Fit_by_Condition_2021_2022.png")

### 3c_Average_TW_Over_Generations.R

fa <- ggplot(Averaged_Full_2021_2022, aes(x = Generation, y = Centered_TW, color = Condition, add = "reg.line")) +
  geom_jitter() +
  geom_smooth(method = lm) +
  geom_hline(yintercept = 0, color = 'red') +
  stat_regline_equation() +
  labs(x = "Generation",
       y = "Average Total Weight (grams)")
ggsave("scripts/plotting/03c_Scatterplot_Avg_TW_by_Condition_2021_2022.png")

### 3ci_Average_FT_Over_Generations.R

fb <- ggplot(Averaged_Full_2021_2022, aes(x = Generation, y = Centered_FT, color = Condition, add = "reg.line")) +
  geom_jitter() +
  geom_smooth(method = lm) +
  geom_hline(yintercept = 0, color = 'red') +
  stat_regline_equation() +
  labs(x = "Generation",
       y = "Average Flowering Time (Days After Sowing)") 
ggsave("scripts/plotting/03ci_Scatterplot_Avg_FT_by_Condition_2021_2022.png")

### 3cii_Average_Fecundity_Over_Generations.R

fc <- ggplot(Averaged_Full_2021_2022, aes(x = Generation, y = Centered_Fec, color = Condition, add = "reg.line")) +
  geom_jitter() +
  geom_smooth(method = lm) +
  geom_hline(yintercept = 0, color = 'red') +
  stat_regline_equation() +
  labs(x = "Generation",
       y = "Average Fecundity")
ggsave("scripts/plotting/03cii_Scatterplot_Avg_Fec_by_Condition_2021_2022.png")

### 3ciii_Average_Fitness_Over_Generations.R

fd <- ggplot(Averaged_Full_2021_2022, aes(x = Generation, y = Centered_Fit, color = Condition, add = "reg.line")) +
  geom_jitter() +
  geom_smooth(method = lm) +
  geom_hline(yintercept = 0, color = 'red') +
  stat_regline_equation() +
  labs(x = "Generation",
       y = "Average Fitness")
ggsave("scripts/plotting/03ciii_Scatterplot_Avg_Fit_by_Condition_2021_2022.png")

### 3ciiii_Combined_Mixed_v_Single_Scatterplots.R

g <- grid.arrange(fa, fb, fc, fd, top = "Evolution of our Four Measured Phenotypes Between Condition 2021-2022 Season")
ggsave("scripts/plotting/03ciiii_Combined_Mixed_v_Single_Scatterplots_2021_2022.png", g, width = 14, height = 10)

### 3cc_T_test_Mixed_vs_Single_Yield

t.test(TOTAL_MASS ~ Condition, Averaged_Full_2021_2022)

### 3cci_T_test_Mixed_vs_Single_Fitness

t.test(Fitness ~ Condition, Averaged_Full_2021_2022)

### 3ccii_T_test_Mixed_vs_Single_Fecundity

t.test(Fecundity ~ Condition, Averaged_Full_2021_2022)

### 3d_Seed_TW_Per_Genotype

ggplot(Averaged_Full_2021_2022, aes(x = reorder(Genotype, + TOTAL_MASS), TOTAL_MASS, fill = Condition, group = Generation)) +
  geom_bar(stat = 'identity', position = position_dodge()) +
  labs(y = "Average Total Seed Weight",
       title = "Average Total Weight of Genotype Across Generations") +
  facet_wrap(~Generation, scales = "free_x") +
  theme(axis.text.x = element_text(angle = 90))

### 3e_Intermediate_FT_reproductive_success.R

# FT vs. Fit
ggplot(Single_2021_2022, aes(FT, Fitness, add = "reg.line")) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ x + I(x^2), size = 1) +
  facet_wrap(~Generation, scales = "free_x")
ggsave("scripts/plotting/03e_Int_FT_vs_Fit_2021_2022.png")

# FT vs. Fec
ggplot(Single_2021_2022, aes(FT, Fecundity)) +
  geom_point() +
  geom_smooth() +
  facet_wrap(~Generation, scales = "free_x")
ggsave("scripts/plotting/03e_Int_FT_vs_fec_2021_2022.png")

# FT vs. Yield
ggplot(Single_2021_2022, aes(FT, TOTAL_MASS)) +
  geom_point()+
  geom_smooth() +
  labs(y = "Yield") +
  facet_wrap(~Generation, scales = "free_x")
ggsave("scripts/plotting/03e_Int_FT_vs_Yield_2021_2022.png")






#### THESE ARE ADDITIONAL GRAPHS FOR MANUSCRIPT (Scatterplot Graphs comparing between the two conditions, then one combined graph)

### Avg_Total_Seed_Weight_Over_Time

fe <- ggplot(Average_Haplo_rep, aes(Generation, `Brown Bag Weight`, add = "reg.line")) +
  geom_jitter(alpha = .3) +
  geom_smooth(method = lm) +
  stat_regline_equation(label.y = 200) +
  labs(x = "Generation",
       y = "Average Total Seed Weight (grams)") +
  facet_wrap(~Condition)

ff <- ggplot(Average_Haplo_rep, aes(Generation, `Brown Bag Weight`, col = Condition)) +
  geom_jitter(alpha = .5) +
  geom_smooth(method = lm) +
  labs(x = "Generation",
       y = "Average Total Seed Weight (grams)") 

Compare_TW <- grid.arrange(fe,ff)
ggsave("scripts/plotting/Extra_3_Average_TW_Comparisons_2021_2022.png", Compare_TW)

### Avg_FT_Over_Generations.R

fe <- ggplot(Average_Haplo_rep, aes(Generation, FT, add = "reg.line")) +
  geom_jitter(alpha = .3) +
  geom_smooth(method = lm) +
  stat_regline_equation(label.y = 130) +
  labs(x = "Generation",
       y = "Average Flowering Time (Days)") +
  facet_wrap(~Condition)

ff <- ggplot(Average_Haplo_rep, aes(Generation, FT, col = Condition)) +
  geom_jitter(alpha = .5) +
  geom_smooth(method = lm) +
  labs(x = "Generation",
       y = "Average Total Seed Weight (grams)") 

Compare_TW <- grid.arrange(fe,ff)
ggsave("scripts/plotting/Extra_3_Average_FT_Comparisons_2021_2022.png", Compare_TW)

### Avg_Fec_Over_Generations.R

fe <- ggplot(Average_Haplo_rep, aes(Generation, Fecundity, add = "reg.line")) +
  geom_jitter(alpha = .3) +
  geom_smooth(method = lm) +
  stat_regline_equation(label.y = 3500) +
  labs(x = "Generation",
       y = "Average Fecundity") +
  facet_wrap(~Condition)

ff <- ggplot(Average_Haplo_rep, aes(Generation, Fecundity, col = Condition)) +
  geom_jitter(alpha = .5) +
  geom_smooth(method = lm) +
  labs(x = "Generation",
       y = "Average Fecundity") 

Compare_TW <- grid.arrange(fe,ff)
ggsave("scripts/plotting/Extra_3_Average_Fecundity_Comparisons_2021_2022.png", Compare_TW)

### Avg_Fit_Over_Generations.R

fe <- ggplot(Average_Haplo_rep, aes(Generation, Fitness, add = "reg.line")) +
  geom_jitter(alpha = .3) +
  geom_smooth(method = lm) +
  stat_regline_equation(label.y = 31000) +
  labs(x = "Generation",
       y = "Average Fitness") +
  facet_wrap(~Condition)

ff <- ggplot(Average_Haplo_rep, aes(Generation, Fitness, col = Condition)) +
  geom_jitter(alpha = .5) +
  geom_smooth(method = lm) +
  labs(x = "Generation",
       y = "Average Fitness") 

Compare_TW <- grid.arrange(fe,ff)
ggsave("scripts/plotting/Extra_3_Average_Fitness_Comparisons.png", Compare_TW)
