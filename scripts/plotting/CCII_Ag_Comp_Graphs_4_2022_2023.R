library(tidyverse)
library(readr)
library(dplyr)
library(ggpubr)
library(ggplot2)
library(googlesheets4)
library(tidyr)

### Load Data

Full_Data <- read_delim("~/Documents/GitHub/Ag-Competition/Full_Data")
Average_Haplo_rep <- read_delim("~/Documents/GitHub/Ag-Competition/Average_Haplo_rep")
Rep_Mixed <- read_delim("~/Documents/GitHub/Ag-Competition/Rep_Mixed")
Rep_Single <- read_delim("~/Documents/GitHub/Ag-Competition/Rep_Single")

### 4_Comparing_Yield_Between_Haplotypes.R

Rep_Single <- Rep_Single %>% filter(Haplotype != "NA")
Rep_Single$Haplotype <- as.factor(Rep_Single$Haplotype)
Haplo_graph <- Rep_Single %>% group_by(Generation, Haplotype) %>% summarise(Avg_TW = mean(`Brown Bag Weight`)) %>% ungroup()
Haplo_graph <- Haplo_graph %>% group_by(Generation) %>% mutate(GenAvg = mean(Avg_TW)) 

ggplot(Haplo_graph,aes(x= Haplotype, y = Avg_TW)) +
  geom_bar(stat = 'identity') +
  geom_hline(yintercept = mean(Haplo_graph$Avg_TW), color = "red") +
  geom_hline(aes(yintercept = GenAvg), color = 'blue', linetype = 2) +
  labs(y = "Average Yield",
       title = "Average Yield between Haplotypes") +
  facet_wrap(~Generation, scales = "free_x") + 
  theme(axis.text.x = element_text(angle = 90)) 
ggsave("scripts/plotting/Average_Yield_Haplotypes.png")

### 4a_Comparing_FT_Between_Haplotypes.R

Haplo_graph <- Rep_Single %>% group_by(Generation, Haplotype) %>% summarise(Avg_FT = mean(FT_DAYS)) %>% ungroup()
Haplo_graph <- Haplo_graph %>% group_by(Generation) %>% mutate(GenAvg = mean(Avg_FT)) 

ggplot(Haplo_graph,aes(x= Haplotype, y = Avg_FT)) +
  geom_bar(stat = 'identity') +
  geom_hline(yintercept = mean(Haplo_graph$Avg_FT), color = "red") +
  geom_hline(aes(yintercept = GenAvg), color = 'blue', linetype = 2) +
  labs(y = "Average Flowering Time (Days)",
       title = "Average Flowering Time between Haplotypes") +
  facet_wrap(~Generation, scales = "free_x") + 
  theme(axis.text.x = element_text(angle = 90)) 
ggsave("scripts/plotting/Average_FT_Haplotypes.png")

### 4b_Comparing_Fecundity_Between_Haplotypes.R

Haplo_graph <- Rep_Single %>% group_by(Generation, Haplotype) %>% summarise(Avg_Fec = mean(Fecundity)) %>% ungroup()
Haplo_graph <- Haplo_graph %>% group_by(Generation) %>% mutate(GenAvg = mean(Avg_Fec)) 

ggplot(Haplo_graph,aes(x= Haplotype, y = Avg_Fec)) +
  geom_bar(stat = 'identity') +
  geom_hline(yintercept = mean(Haplo_graph$Avg_Fec), color = "red") +
  geom_hline(aes(yintercept = GenAvg), color = 'blue', linetype = 2) +
  labs(y = "Average Fecundity",
       title = "Average Fecundity between Haplotypes") +
  facet_wrap(~Generation, scales = "free_x") + 
  theme(axis.text.x = element_text(angle = 90)) 
ggsave("scripts/plotting/Average_Fec_Haplotypes.png")

### 4c_Comparing_Fitness_Between_Haplotypes.R

Haplo_graph <- Rep_Single %>% group_by(Generation, Haplotype) %>% summarise(Avg_Fit = mean(Fitness)) %>% ungroup()
Haplo_graph <- Haplo_graph %>% group_by(Generation) %>% mutate(GenAvg = mean(Avg_Fit)) 

ggplot(Haplo_graph,aes(x= Haplotype, y = Avg_Fit)) +
  geom_bar(stat = 'identity') +
  geom_hline(yintercept = mean(Haplo_graph$Avg_Fit), color = "red") +
  geom_hline(aes(yintercept = GenAvg), color = 'blue', linetype = 2) +
  labs(y = "Average Fitness",
       title = "Average Fitness between Haplotypes") +
  facet_wrap(~Generation, scales = "free_x") + 
  theme(axis.text.x = element_text(angle = 90)) 
ggsave("scripts/plotting/Average_Fit_Haplotypes.png")


  
for (i in unique(T$Generation)){
  tmp <-  ggplot(subset(T, Generation == i), aes(Genotypes, TW)) + geom_bar(stat = 'identity') + 
    theme(axis.text.x = element_text(angle = 60)) + labs(y = "Total Weight (g)") +
    ggtitle(paste0("Generation ", i))
  ggsave(tmp, file = paste0("Bar_Plot_Total_Weight_Generation_", i, ".png"))
}


ggplot(Average_Haplo_rep, aes(x = Haplotype, y = `Brown Bag Weight`)) +
  geom_bar(stat = 'identity') +
  facet_wrap(~Generation) +
  labs(y = "Average Yield",
       title = "Average Yield between Haplotypes")

