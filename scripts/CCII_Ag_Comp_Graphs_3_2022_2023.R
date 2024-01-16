library(tidyverse)
library(readr)
library(dplyr)
library(ggpubr)
library(ggplot2)
library(tidyr)
library(car)

### Load Data


Full_Data <- read_delim("~/Documents/GitHub/Ag-Competition/Full_Data")
Side_By_Side_Replicates <- read_delim("~/Documents/GitHub/Ag-Competition/Side_By_Side_Replicates")
Rep_Mixed <- read_delim("~/Documents/GitHub/Ag-Competition/Rep_Mixed")
Rep_Single <- read_delim("~/Documents/GitHub/Ag-Competition/Rep_Single")

### 3a_Bar_Graph_to_Compare_Yield_Between_Genotypes.R




ggplot(Side_By_Side_Replicates, aes(Genotypes, Avg_Total_Weight, color = Condition)) +
  geom_bar(stat = 'identity') +
  labs(y = "Average Total Weight (grams)",
       title = "Comparing Average Total Weight Among Genotypes") 

### 3b_Total_Seed_Weight_Exp_by_Condition.R
ggplot(Side_By_Side_Replicates, aes(Genotypes, Exp_TW_Per_Plant)) +
  geom_bar(stat = 'identity') 


  
### 3c_Average_Total_Weight_Over_Generations.R

Side_By_Side_Replicates$Generation <- as.numeric(Side_By_Side_Replicates$Generation)
ggplot(Side_By_Side_Replicates, aes(x = Generation, y = Avg_Total_Weight, color = Condition)) +
  geom_jitter() +
  geom_smooth(method = lm) +
  labs(x = "Generation",
       y = "Average Total Weight (grams)",
       title = "Average Total Weight Over Generations")

### 3cI_T_test_Mixed_vs_Single_Yield

t.test(Avg_Total_Weight ~ Condition, Side_By_Side_Replicates)

### 3d_Seed_Yield_Per_Genotype

Sorted_Genotypes <- Side_By_Side_Replicates %>% group_by(Genotypes) %>% arrange(Avg_Total_Weight, .by_group = T)


ggplot(Sorted_Genotypes, aes(Genotypes, Avg_Total_Weight, fill = Condition)) +
  geom_bar(stat = 'identity') +
  facet_grid(~Generation) +
  theme(axis.text.x = element_text(angle = 60)) +
  labs(y = "Average Total Seed Weight",
       title = "Average Total Weight of Genotypes Across Generations") +



  







