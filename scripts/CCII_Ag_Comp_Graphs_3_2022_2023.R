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

### 3a_Bar_Graph_to_Compare_Avg_Yield_Between_Genotypes.R

q <- ggplot(Average_Haplo_rep, aes(Genotypes, `Brown Bag Weight`, color = Condition, fill = Condition)) +
  geom_bar(stat = 'identity') +
  labs(y = "Average Total Weight (grams)",
       title = "Comparing Average Total Weight Among Genotypes") +
  coord_flip() +
  scale_y_discrete(breaks = Average_Haplo_rep$Genotypes,
                   labels = Average_Haplo_rep$Genotypes,
                   guide = guide_axis(n.dodge = 2)) 
q + theme(axis.text.y = element_text(angle = 0, vjust = 0, hjust =0))

### 3b_Exp_TW_Per_Plant_by_Condition.R

ggplot(Average_Haplo_rep, aes(Generation, Exp_TW_Per_Plant, color = Condition, add = "reg.line")) +
  geom_jitter() +
  geom_smooth(method = lm) +
  stat_regline_equation() +
  labs(y = "Expected Yield Per Plant (g)",
       title = "Generational Change in Expected Yield Per Plant")

### 3bi_Exp_Fec_Per_Plant_by_Condition.R

ggplot(Average_Haplo_rep, aes(Generation, Exp_Fec_Per_Plant, color = Condition, add = "reg.line")) +
  geom_jitter() +
  geom_smooth(method = lm) +
  stat_regline_equation() +
  labs(y = "Expected Fecundity Per Plant",
       title = "Generational Change in Expected Fecundity Per Plant")

### 3bii_Exp_Fit_by_Condition.R

ggplot(Average_Haplo_rep, aes(Generation, Exp_Fit_Per_Plant, color = Condition, add = "reg.line")) +
  geom_jitter() +
  geom_smooth(method = lm) +
  stat_regline_equation() +
  labs(y = "Expected Fitness Per Plant",
       title = "Generational Change in Expected Fitness Per Plant")

### 3c_Average_Total_Weight_Over_Generations.R

Average_Haplo_rep$Generation <- as.numeric(Average_Haplo_rep$Generation)
ggplot(Average_Haplo_rep, aes(x = Generation, y = `Brown Bag Weight`, color = Condition)) +
  geom_jitter() +
  geom_smooth(method = lm) +
  labs(x = "Generation",
       y = "Average Total Weight (grams)",
       title = "Average Total Weight Over Generations")

### 3ci_T_test_Mixed_vs_Single_Yield

t.test(`Brown Bag Weight` ~ Condition, Average_Haplo_rep)

### 3cii_T_test_Mixed_vs_Single_Fitness

t.test(Fitness ~ Condition, Average_Haplo_rep)

### 3ciii_T_test_Mixed_vs_Single_Fecundity

t.test(Fecundity ~ Condition, Average_Haplo_rep)

### 3d_Seed_Yield_Per_Genotype

Sorted_Genotypes <- Average_Haplo_rep %>% arrange(Generation, `Brown Bag Weight`)

ggplot(Sorted_Genotypes, aes(Genotypes, `Brown Bag Weight`, fill = Condition, group = Generation)) +
  geom_bar(stat = 'identity') +
  theme(axis.text.x = element_text(angle = 60)) +
  labs(y = "Average Total Seed Weight",
       title = "Average Total Weight of Genotypes Across Generations") 


### 3e_Intermediate_FT_reproductive_success.R

# FT vs. Fit
ggplot(Rep_Single, aes(FT_DAYS, Fitness)) +
  geom_point() +
  geom_smooth() +
  facet_wrap(~Generation, scales = "free_x")

# FT vs. Fec
ggplot(Rep_Single, aes(FT_DAYS, Fecundity)) +
  geom_point() +
  geom_smooth() +
  facet_wrap(~Generation, scales = "free_x")

# FT vs. Yield
ggplot(Rep_Single, aes(FT_DAYS, `Brown Bag Weight`)) +
  geom_point()+
  geom_smooth() +
  facet_wrap(~Generation, scales = "free_x")

