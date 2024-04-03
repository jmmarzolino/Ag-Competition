library(tidyverse)
library(readr)
library(dplyr)
library(ggpubr)
library(ggplot2)
library(tidyr)
library(car)

### Load Data

Full_Data <- read_delim("~/Documents/GitHub/Ag-Competition/Full_Data")
Average_Haplo_rep <- read_delim("~/Documents/GitHub/Ag-Competition/Average_Haplo_rep")
Rep_Mixed <- read_delim("~/Documents/GitHub/Ag-Competition/Rep_Mixed")
Rep_Single <- read_delim("~/Documents/GitHub/Ag-Competition/Rep_Single")


test <- Average_Haplo_rep %>% select(c("Genotypes", "Generation", "Condition", "Brown Bag Weight", "Fitness", "Fecundity", "FT_DAYS")) %>% group_by(Genotypes) %>% 
  pivot_wider(names_from = "Condition", values_from = c("Brown Bag Weight", "Fecundity", "Fitness"))
test <- ifelse(test$)


test$new <- ifelse(test$`Brown Bag Weight_single` > test$`Brown Bag Weight_single`, 1,0)


### 3a_Bar_Graph_Avg_Yield_Between_Genotypes.R

ggplot(Average_Haplo_rep, aes(Genotypes, `Brown Bag Weight`, color = Condition, fill = Condition)) +
  geom_bar(stat = 'identity', position = position_dodge(), alpha = .5, width = .5) +
  labs(y = "Yield (grams)",
       title = "Comparing Yield Between Individual Genotypes") +
  facet_wrap(~Generation, scales = "free_x") +
  theme(axis.text.x = element_text(size = 5)) +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_y_continuous(breaks = seq(0, 200, 10))
ggsave("scripts/plotting/03a_Bar_Graph_Yield_Between_Conditions.png", width = 18, height = 10)

### 3ai_Bar_Graph_Avg_FT_Between_Genotypes.R

ggplot(Average_Haplo_rep, aes(Genotypes, FT_DAYS, color = Condition, fill = Condition)) +
  geom_bar(stat = 'identity', position = position_dodge(), alpha = .3, width = .5) +
  labs(y = "Flowering Time (Days)",
       title = "Comparing Flowering Time Between Individual Genotypes") +
  facet_wrap(~Generation, scales = "free_x") +
  theme(axis.text.x = element_text(size = 5)) +
  theme(axis.text.x = element_text(angle = 90))
ggsave("scripts/plotting/03ai_Bar_Graph_FT_Between_Conditions.png", width = 18, height = 10)

### 3aii_Bar_Graph_Avg_Fec_Between_Genotypes.R

ggplot(Average_Haplo_rep, aes(Genotypes, Fecundity, color = Condition, fill = Condition)) +
  geom_bar(stat = 'identity', position = position_dodge(), alpha = .3, width = .5) +
  labs(y = "Fecundity",
       title = "Comparing Fecundity Between Individual Genotypes") +
  facet_wrap(~Generation, scales = "free_x") +
  theme(axis.text.x = element_text(size = 5)) +
  theme(axis.text.x = element_text(angle = 90))
ggsave("scripts/plotting/03aii_Bar_Graph_Fec_Between_Conditions.png", width = 18, height = 10)

### 3aiii_Bar_Graph_Avg_Fit_Between_Genotypes.R

ggplot(Average_Haplo_rep, aes(Genotypes, Fitness, color = Condition, fill = Condition)) +
  geom_bar(stat = 'identity', position = position_dodge(), alpha = .3, width = .5) +
  labs(y = "Fitness",
       title = "Comparing Fitness Between Individual Genotypes") +
  facet_wrap(~Generation, scales = "free_x") +
  theme(axis.text.x = element_text(size = 5)) +
  theme(axis.text.x = element_text(angle = 90))
ggsave("scripts/plotting/03aiii_Bar_Graph_Fit_Between_Conditions.png",  width = 18, height = 10)

### 3b_Exp_Yield_Per_Plant_by_Condition.R

ggplot(Average_Haplo_rep, aes(Generation, Exp_TW_Per_Plant, color = Condition, add = "reg.line")) +
  geom_jitter() +
  geom_smooth(method = lm) +
  stat_regline_equation() +
  labs(y = "Expected Yield Per Plant (g)",
       title = "Generational Change in Expected Yield Per Plant")
ggsave("scripts/plotting/03b_Exp_Per_Plant_Yield_by_Condition.png")

### 3bi_Exp_Fec_Per_Plant_by_Condition.R

ggplot(Average_Haplo_rep, aes(Generation, Exp_Fec_Per_Plant, color = Condition, add = "reg.line")) +
  geom_jitter() +
  geom_smooth(method = lm) +
  stat_regline_equation() +
  labs(y = "Expected Fecundity Per Plant",
       title = "Generational Change in Expected Fecundity Per Plant")
ggsave("scripts/plotting/03bi_Exp_Per_Plant_Fec_by_Condition.png")

### 3bii_Exp_Fit_by_Condition.R

ggplot(Average_Haplo_rep, aes(Generation, Exp_Fit_Per_Plant, color = Condition, add = "reg.line")) +
  geom_jitter() +
  geom_smooth(method = lm) +
  stat_regline_equation() +
  labs(y = "Expected Fitness Per Plant",
       title = "Generational Change in Expected Fitness Per Plant")
ggsave("scripts/plotting/03bii_Exp_Per_Plant_Fit_by_Condition.png")

### 3c_Average_Yield_Over_Generations.R

Average_Haplo_rep$Generation <- as.numeric(Average_Haplo_rep$Generation)
fg <- ggplot(Average_Haplo_rep, aes(x = Generation, y = `Brown Bag Weight`, color = Condition, add = "reg.line")) +
  geom_jitter() +
  geom_smooth(method = lm) +
  stat_regline_equation() +
  labs(x = "Generation",
       y = "Average Yield (grams)")
ggsave("scripts/plotting/03c_Scatterplot_Avg_Yield_by_Condition.png")

### 3ci_Average_FT_Over_Generations.R

fa <- ggplot(Average_Haplo_rep, aes(x = Generation, y = FT_DAYS, color = Condition, add = "reg.line")) +
  geom_jitter() +
  geom_smooth(method = lm) +
  stat_regline_equation() +
  labs(x = "Generation",
       y = "Average Flowering Time (Days After Sowing)") +
  scale_y_continuous(breaks = seq(0, 140, 10))
ggsave("scripts/plotting/03ci_Scatterplot_Avg_FT_by_Condition.png")

### 3cii_Average_Fecundity_Over_Generations.R

fb <- ggplot(Average_Haplo_rep, aes(x = Generation, y = Fecundity, color = Condition, add = "reg.line")) +
  geom_jitter() +
  geom_smooth(method = lm) +
  stat_regline_equation() +
  labs(x = "Generation",
       y = "Average Fecundity")
ggsave("scripts/plotting/03cii_Scatterplot_Avg_Fec_by_Condition.png")

### 3ciii_Average_Fitness_Over_Generations.R

fv <- ggplot(Average_Haplo_rep, aes(x = Generation, y = Fitness, color = Condition, add = "reg.line")) +
  geom_jitter() +
  geom_smooth(method = lm) +
  stat_regline_equation() +
  labs(x = "Generation",
       y = "Average Fitness")
ggsave("scripts/plotting/03ciii_Scatterplot_Avg_Fit_by_Condition.png")

### 3ciiii_Combined_Mixed_v_Single_Scatterplots.R

g <- grid.arrange(fg, fa, fb, fv, top = "Comparing Evolution of our Four Measured Phenotypes Between Condition")
ggsave("scripts/plotting/03ciiii_Combined_Mixed_v_Single_Scatterplots.png", g, width = 14, height = 10)

### 3cc_T_test_Mixed_vs_Single_Yield

t.test(`Brown Bag Weight` ~ Condition, Average_Haplo_rep)

### 3cci_T_test_Mixed_vs_Single_Fitness

t.test(Fitness ~ Condition, Average_Haplo_rep)

### 3ccii_T_test_Mixed_vs_Single_Fecundity

t.test(Fecundity ~ Condition, Average_Haplo_rep)

### 3d_Seed_Yield_Per_Genotype

ggplot(Average_Haplo_rep, aes(x = reorder(Genotypes, +`Brown Bag Weight`), `Brown Bag Weight`, fill = Condition, group = Generation)) +
  geom_bar(stat = 'identity', position = position_dodge()) +
  labs(y = "Average Total Seed Weight",
       title = "Average Total Weight of Genotypes Across Generations") +
  facet_wrap(~Generation, scales = "free_x") +
  theme(axis.text.x = element_text(angle = 90))

### 3e_Intermediate_FT_reproductive_success.R

# FT vs. Fit
ggplot(Rep_Single, aes(FT_DAYS, Fitness, add = "reg.line")) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ x + I(x^2), size = 1) +
  facet_wrap(~Generation, scales = "free_x")
ggsave("scripts/plotting/03e_Int_FT_vs_Fit.png")

# FT vs. Fec
ggplot(Rep_Single, aes(FT_DAYS, Fecundity)) +
  geom_point() +
  geom_smooth() +
  facet_wrap(~Generation, scales = "free_x")
ggsave("scripts/plotting/03e_Int_FT_vs_fec.png")

# FT vs. Yield
ggplot(Rep_Single, aes(FT_DAYS, `Brown Bag Weight`)) +
  geom_point()+
  geom_smooth() +
  labs(y = "Yield") +
  facet_wrap(~Generation, scales = "free_x")
ggsave("scripts/plotting/03e_Int_FT_vs_Yield.png")






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
ggsave("scripts/plotting/Extra_3_Average_Yield_Comparisons.png", Compare_TW)

### Avg_FT_Over_Generations.R

fe <- ggplot(Average_Haplo_rep, aes(Generation, FT_DAYS, add = "reg.line")) +
  geom_jitter(alpha = .3) +
  geom_smooth(method = lm) +
  stat_regline_equation(label.y = 130) +
  labs(x = "Generation",
       y = "Average Flowering Time (Days)") +
  facet_wrap(~Condition)

ff <- ggplot(Average_Haplo_rep, aes(Generation, FT_DAYS, col = Condition)) +
  geom_jitter(alpha = .5) +
  geom_smooth(method = lm) +
  labs(x = "Generation",
       y = "Average Total Seed Weight (grams)") 

Compare_TW <- grid.arrange(fe,ff)
ggsave("scripts/plotting/Extra_3_Average_FT_Comparisons.png", Compare_TW)

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
ggsave("scripts/plotting/Extra_3_Average_Fecundity_Comparisons.png", Compare_TW)

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
