library(tidyverse)
library(readr)
library(dplyr)
library(ggpubr)
library(ggplot2)
library(googlesheets4)
library(tidyr)

### 5a_Comparing_TW_Between_Haplotypes.R

Single_2021_2022 <- Single_2021_2022 %>% filter(Haplotype != "NA")
Single_2021_2022$Haplotype <- as.factor(Single_2021_2022$Haplotype)
Haplo_graph <- Single_2021_2022 %>% group_by(Generation, Haplotype) %>% summarise(Avg_TW = mean(total_seed_mass_g)) %>% ungroup()
Haplo_graph <- Haplo_graph %>% group_by(Generation) %>% mutate(GenAvg = mean(Avg_TW)) %>% ungroup()
Haplo_graph$`Avg Yield > Pop & Gen Avg Yield` <- ifelse(Haplo_graph$Avg_TW > Haplo_graph$GenAvg & Haplo_graph$Avg_TW > mean(Haplo_graph$Avg_TW), TRUE, FALSE)
f_AY <- Haplo_graph %>% filter(`Avg Yield > Pop & Gen Avg Yield` == TRUE) %>% filter(!duplicated(Haplotype)) %>% select(Haplotype)

ggplot(Haplo_graph,aes(x= Haplotype, y = Avg_TW, fill = `Avg Yield > Pop & Gen Avg Yield`)) +
  geom_bar(stat = 'identity') +
  geom_hline(yintercept = mean(Haplo_graph$Avg_TW), color = "green") +
  geom_hline(aes(yintercept = GenAvg), color = 'blue', linetype = 2) +
  labs(y = "Average Total Seed Weight",
       title = "Average Total Seed Weight between Haplotypes") +
  scale_y_continuous(breaks = seq(0, 200, 25)) +
  facet_wrap(~Generation, scales = "free_x") + 
  theme(axis.text.x = element_text(angle = 90)) +
  scale_fill_manual(values = c("black", "red")) 
ggsave("scripts/plotting/05a_Average_TW_Haplotypes_2021_2022.png", width = 16, height = 12)

### 5ai_Comparing_FT_Between_Haplotypes.R

Haplo_graph <- Single_2021_2022 %>% group_by(Generation, Haplotype) %>% summarise(Avg_FT = mean(Flowering_Date)) %>% ungroup()
Haplo_graph <- Haplo_graph %>% group_by(Generation) %>% mutate(GenAvg = mean(Avg_FT)) %>% ungroup()
Haplo_graph$`Avg FT > Pop & Gen Avg FT` <- ifelse(Haplo_graph$Avg_FT > Haplo_graph$GenAvg & Haplo_graph$Avg_FT > mean(Haplo_graph$Avg_FT), TRUE, FALSE)
f_FT <- Haplo_graph %>% filter(`Avg FT > Pop & Gen Avg FT` == TRUE) %>% filter(!duplicated(Haplotype)) %>% select(Haplotype)

ggplot(Haplo_graph,aes(x= Haplotype, y = Avg_FT, fill = `Avg FT > Pop & Gen Avg FT`)) +
  geom_bar(stat = 'identity') +
  geom_hline(yintercept = mean(Haplo_graph$Avg_FT), color = "green") +
  geom_hline(aes(yintercept = GenAvg), color = 'blue', linetype = 2) +
  labs(y = "Average Flowering Time (Days)",
       title = "Average Flowering Time between Haplotypes") +
  facet_wrap(~Generation, scales = "free_x") + 
  theme(axis.text.x = element_text(angle = 90)) +
  scale_fill_manual(values = c("black", "red"))
ggsave("scripts/plotting/05ai_Average_FT_Haplotypes_2021_2022.png", width = 16, height = 12)

### 5aii_Comparing_Fecundity_Between_Haplotypes.R

Haplo_graph <- Single_2021_2022 %>% group_by(Generation, Haplotype) %>% summarise(Avg_Fec = mean(Fecundity)) %>% ungroup()
Haplo_graph <- Haplo_graph %>% filter(Avg_Fec != "NA")
Haplo_graph <- Haplo_graph %>% group_by(Generation) %>% mutate(GenAvg = mean(Avg_Fec)) %>% ungroup()
Haplo_graph$`Avg Fec > Pop & Gen Avg Fec` <- ifelse(Haplo_graph$Avg_Fec > Haplo_graph$GenAvg & Haplo_graph$Avg_Fec > mean(Haplo_graph$Avg_Fec), TRUE, FALSE)
f_Fe <- Haplo_graph %>% filter(`Avg Fec > Pop & Gen Avg Fec` == TRUE) %>% filter(!duplicated(Haplotype)) %>% select(Haplotype)

ggplot(Haplo_graph,aes(x= Haplotype, y = Avg_Fec, fill = `Avg Fec > Pop & Gen Avg Fec`)) +
  geom_bar(stat = 'identity') +
  geom_hline(yintercept = mean(Haplo_graph$Avg_Fec), color = "green") +
  geom_hline(aes(yintercept = GenAvg), color = 'blue', linetype = 2) +
  labs(y = "Average Fecundity",
       title = "Average Fecundity between Haplotypes") +
  facet_wrap(~Generation, scales = "free_x") + 
  theme(axis.text.x = element_text(angle = 90)) +
  scale_fill_manual(values = c("black", "red"))
ggsave("scripts/plotting/05aii_Average_Fec_Haplotypes_2021_2022.png", width = 16, height = 12)

### 5aiii_Comparing_Fitness_Between_Haplotypes.R

Haplo_graph <- Single_2021_2022 %>% group_by(Generation, Haplotype) %>% summarise(Avg_Fit = mean(Fitness)) %>% ungroup()
Haplo_graph <- Haplo_graph %>% filter(Avg_Fit != "NA")
Haplo_graph <- Haplo_graph %>% group_by(Generation) %>% mutate(GenAvg = mean(Avg_Fit)) %>% ungroup()
Haplo_graph$`Avg Fit > Pop & Gen Avg Fit` <- ifelse(Haplo_graph$Avg_Fit > Haplo_graph$GenAvg & Haplo_graph$Avg_Fit > mean(Haplo_graph$Avg_Fit), TRUE, FALSE)
f_Fi <- Haplo_graph %>% filter(`Avg Fit > Pop & Gen Avg Fit` == TRUE) %>% filter(!duplicated(Haplotype)) %>% select(Haplotype)
ui <- Reduce(intersect, list(f_AY, f_Fe, f_Fi))

ggplot(Haplo_graph,aes(x= Haplotype, y = Avg_Fit, fill = `Avg Fit > Pop & Gen Avg Fit`)) +
  geom_bar(stat = 'identity') +
  geom_hline(yintercept = mean(Haplo_graph$Avg_Fit, na.rm = TRUE), color = "green") +
  geom_hline(aes(yintercept = GenAvg), color = 'blue', linetype = 2) +
  labs(y = "Average Fitness",
       title = "Average Fitness between Haplotypes") +
  facet_wrap(~Generation, scales = "free_x") + 
  theme(axis.text.x = element_text(angle = 90)) +
  scale_fill_manual(values = c("black", "red"))
ggsave("scripts/plotting/05aiii_Average_Fit_Haplotypes_2021_2022.png", width = 16, height = 12)



for (i in unique(T$Generation)){
  tmp <-  ggplot(subset(T, Generation == i), aes(Genotypes, TW)) + geom_bar(stat = 'identity') + 
    theme(axis.text.x = element_text(angle = 60)) + labs(y = "Total Weight (g)") +
    ggtitle(paste0("Generation ", i))
  ggsave(tmp, file = paste0("Bar_Plot_Total_Weight_Generation_", i, ".png"))
}
