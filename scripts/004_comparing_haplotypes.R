#!/usr/bin/env Rscript
#SBATCH --mem=30G
#SBATCH --time=02:00:00
#SBATCH --output=/rhome/jmarz001/bigdata/Ag-Competition/scripts/004_comparing_haplotypes.stdout
#SBATCH -p koeniglab

library(tidyverse)
library(ggpubr)
library(data.table)

setwd("/bigdata/koeniglab/jmarz001/Ag-Competition/data")
source("../scripts/CUSTOM_FNS.R")
hap <- fread("Competition_Lines_Haplotypes.csv")
df <- fread("FITNESS.tsv")


### Importing Haplotype Data and cleaning the dataframe
hap <- hap %>% select(c("Haplotype", "Generation", "Family"))
colnames(hap)[which(names(hap) == "Family")] <- "Genotype"
hap <- hap %>% filter(!is.na(Haplotype))
#hap$Genotype <- gsub("(\\d+_\\d+)_\\d", "\\1", hap$Genotype)

# cut 3-digit genotype codes down to 2 to more easily match families
df$Genotype <- gsub("(\\d+_\\d+)_\\d", "\\1", df$Genotype)

#setdiff(df$Genotype, hap$Genotype)
hap_join <- inner_join(hap, df, by = c("Genotype", "Generation"))
#anti_join(hap, df, by = c("Genotype", "Generation"))
write_delim(hap_join, "HAPLOTYPE_PHENOS.tsv")



# plot traits by haplotypes
df2 <- hap_join %>% filter(Condition == "single") %>% select(-Condition)
table(df2$Haplotype)
hap_list <- names(which(table(df2$Haplotype) > 1))
hap_df <- df2 %>% filter(Haplotype %in% hap_list)





### 5a_Comparing_TW_Between_Haplotypes.R

Single$Haplotype <- as.factor(Single$Haplotype)
Haplo_graph <- Single %>% group_by(Generation, Haplotype) %>% summarise(Avg_TW = mean(TOTAL_MASS, na.rm = T)) %>% ungroup()
Haplo_graph <- Haplo_graph %>% group_by(Generation) %>% mutate(GenAvg = mean(Avg_TW, na.rm = T)) %>% ungroup()
Haplo_graph$`Avg Yield > Pop & Gen Avg Yield` <- ifelse(Haplo_graph$Avg_TW > Haplo_graph$GenAvg & Haplo_graph$Avg_TW > mean(Haplo_graph$Avg_TW, na.rm = T), TRUE, FALSE)
f_AY <- Haplo_graph %>% filter(`Avg Yield > Pop & Gen Avg Yield` == TRUE) %>% filter(!duplicated(Haplotype)) %>% select(Haplotype)

ggplot(Haplo_graph,aes(x= Haplotype, y = Avg_TW, fill = `Avg Yield > Pop & Gen Avg Yield`)) +
  geom_bar(stat = "identity") +
  geom_hline(yintercept = mean(Haplo_graph$Avg_TW, na.rm = T), color = "green") +
  geom_hline(aes(yintercept = GenAvg), color = "blue", linetype = 2) +
  labs(y = "Average Total Seed Weight",
       title = "Average Total Seed Weight between Haplotypes") +
  scale_y_continuous(breaks = seq(0, 200, 25)) +
  facet_wrap(~Generation, scales = "free_x") +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_fill_manual(values = c("black", "red"))
ggsave("scripts/plotting/05a_Average_TW_Haplotypes.png", width = 16, height = 12)

### 5ai_Comparing_FT_Between_Haplotypes.R

Haplo_graph <- Single %>% group_by(Generation, Haplotype) %>% summarise(Avg_FT = mean(FT, na.rm = T)) %>% ungroup()
Haplo_graph <- Haplo_graph %>% group_by(Generation) %>% mutate(GenAvg = mean(Avg_FT, na.rm = T)) %>% ungroup()
Haplo_graph$`Avg FT > Pop & Gen Avg FT` <- ifelse(Haplo_graph$Avg_FT > Haplo_graph$GenAvg & Haplo_graph$Avg_FT > mean(Haplo_graph$Avg_FT, na.rm = T), TRUE, FALSE)
f_FT <- Haplo_graph %>% filter(`Avg FT > Pop & Gen Avg FT` == TRUE) %>% filter(!duplicated(Haplotype)) %>% select(Haplotype)

ggplot(Haplo_graph,aes(x= Haplotype, y = Avg_FT, fill = `Avg FT > Pop & Gen Avg FT`)) +
  geom_bar(stat = "identity") +
  geom_hline(yintercept = mean(Haplo_graph$Avg_FT, na.rm = T), color = "green") +
  geom_hline(aes(yintercept = GenAvg), color = "blue", linetype = 2) +
  labs(y = "Average Flowering Time (Days)",
       title = "Average Flowering Time between Haplotypes") +
  facet_wrap(~Generation, scales = "free_x") +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_fill_manual(values = c("black", "red"))
ggsave("scripts/plotting/05ai_Average_FT_Haplotypes.png", width = 16, height = 12)

### 5aii_Comparing_Fecundity_Between_Haplotypes.R

Haplo_graph <- Single %>% group_by(Generation, Haplotype) %>% summarise(Avg_Fec = mean(FECUNDITY, na.rm = T)) %>% ungroup()
Haplo_graph <- Haplo_graph %>% group_by(Generation) %>% mutate(GenAvg = mean(Avg_Fec, na.rm = T)) %>% ungroup()
Haplo_graph$`Avg Fec > Pop & Gen Avg Fec` <- ifelse(Haplo_graph$Avg_Fec > Haplo_graph$GenAvg & Haplo_graph$Avg_Fec > mean(Haplo_graph$Avg_Fec, na.rm = T), TRUE, FALSE)
f_Fe <- Haplo_graph %>% filter(`Avg Fec > Pop & Gen Avg Fec` == TRUE) %>% filter(!duplicated(Haplotype)) %>% select(Haplotype)

ggplot(Haplo_graph,aes(x= Haplotype, y = Avg_Fec, fill = `Avg Fec > Pop & Gen Avg Fec`)) +
  geom_bar(stat = "identity") +
  geom_hline(yintercept = mean(Haplo_graph$Avg_Fec, na.rm = T), color = "green") +
  geom_hline(aes(yintercept = GenAvg), color = "blue", linetype = 2) +
  labs(y = "Average Centered Fecundity",
       title = "Average Fecundity between Haplotypes") +
  facet_wrap(~Generation, scales = "free_x") +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_fill_manual(values = c("black", "red"))
ggsave("scripts/plotting/05aii_Average_Fec_Haplotypes.png", width = 16, height = 12)

### 5aiii_Comparing_Fitness_Between_Haplotypes.R

Haplo_graph <- Single %>% group_by(Generation, Haplotype) %>% summarise(Avg_Fit = mean(ABS_FITNESS, na.rm = T)) %>% ungroup()
Haplo_graph <- Haplo_graph %>% group_by(Generation) %>% mutate(GenAvg = mean(Avg_Fit, na.rm = T)) %>% ungroup()
Haplo_graph$`Avg Fit > Pop & Gen Avg Fit` <- ifelse(Haplo_graph$Avg_Fit > Haplo_graph$GenAvg & Haplo_graph$Avg_Fit > mean(Haplo_graph$Avg_Fit, na.rm = T), TRUE, FALSE)
f_Fi <- Haplo_graph %>% filter(`Avg Fit > Pop & Gen Avg Fit` == TRUE) %>% filter(!duplicated(Haplotype)) %>% select(Haplotype)
ui <- Reduce(intersect, list(f_AY, f_Fe, f_Fi))

ggplot(Haplo_graph,aes(x= Haplotype, y = Avg_Fit, fill = `Avg Fit > Pop & Gen Avg Fit`)) +
  geom_bar(stat = "identity") +
  geom_hline(yintercept = mean(Haplo_graph$Avg_Fit, na.rm = TRUE), color = "green") +
  geom_hline(aes(yintercept = GenAvg), color = "blue", linetype = 2) +
  labs(y = "Average Fitness",
       title = "Average Absolute Fitness between Haplotypes") +
  facet_wrap(~Generation, scales = "free_x") +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_fill_manual(values = c("black", "red"))
ggsave("scripts/plotting/05aiii_Average_Fit_Haplotypes.png", width = 16, height = 12)

### 5aiiii_Comparing_100SW_Between_Haplotypes.R

Haplo_graph <- Single %>% group_by(Generation, Haplotype) %>% summarise(Avg_100SW = mean(SEED_WEIGHT_100, na.rm = T)) %>% ungroup()
Haplo_graph <- Haplo_graph %>% group_by(Generation) %>% mutate(GenAvg = mean(Avg_100SW, na.rm = T)) %>% ungroup()
Haplo_graph$`Avg Fit > Pop & Gen Avg Fit` <- ifelse(Haplo_graph$Avg_100SW > Haplo_graph$GenAvg & Haplo_graph$Avg_100SW > mean(Haplo_graph$Avg_100SW, na.rm = T), TRUE, FALSE)
f_Fd <- Haplo_graph %>% filter(`Avg Fit > Pop & Gen Avg Fit` == TRUE) %>% filter(!duplicated(Haplotype)) %>% select(Haplotype)
ui <- Reduce(intersect, list(f_AY, f_Fe, f_Fi, f_Fd))

ggplot(Haplo_graph,aes(x= Haplotype, y = Avg_100SW, fill = `Avg Fit > Pop & Gen Avg Fit`)) +
  geom_bar(stat = "identity") +
  geom_hline(yintercept = mean(Haplo_graph$Avg_100SW, na.rm = TRUE), color = "green") +
  geom_hline(aes(yintercept = GenAvg), color = "blue", linetype = 2) +
  labs(y = "Average 100 Seed Weight (grams)",
       title = "Average 100 Seed Weigh between Haplotypes") +
  facet_wrap(~Generation, scales = "free_x") +
  theme(axis.text.x = element_text(angle = 90)) +
  scale_fill_manual(values = c("black", "red"))
ggsave("scripts/plotting/05aiii_Average_Fit_Haplotypes.png", width = 16, height = 12)














### Atlas

Average_Haplo_rep <- Average_Haplo_rep %>%  mutate(Atlas_Avg_Fec = 2363.51,
                  Atlas_Avg_Fitness = 21347.22,
                  Atlas_Avg_Total_Weight = 126.8267)
Average_Haplo_rep$Numbers <- ifelse(Average_Haplo_rep$Condition == "mixed", 1, 0)


### Adding Averaged Atlas values into the table and adding columns for centered data

Atlas_tbl <- sw_avg %>% filter(Genotype == "48_5") %>%
  mutate(Atlas_Avg_Fec = mean(FECUNDITY),
         Atlas_Avg_Fit = mean(FITNESS),
         Atlas_Avg_TW = mean(TOTAL_MASS))

sw_avg <- sw_avg %>% mutate(Atlas_Avg_Fec = 1413.799,
                                                              Atlas_Avg_Fit = 12880.63,
                                                              Atlas_Avg_TW = 64.55571)

sw_avg <- sw_avg %>% mutate(Centered_Fit = FITNESS - mean(FITNESS, na.rm = TRUE),
                                                  Centered_FT = FT - mean(FT, na.rm = TRUE),
                                                  Centered_Fec = FECUNDITY - mean(FECUNDITY, na.rm = TRUE),
                                                  Centered_TW = TOTAL_MASS - mean(TOTAL_MASS))

### Calculating Contribution of each seed to phenotype

### FECUNDITY
sw_avg$Exp_Fec_Per_Plant <- ifelse(sw_avg$Condition == "mixed",
                                              Exp_Fec_Mixed(sw_avg$FEC),
                                              Exp_Single(sw_avg$FEC))

### FITNESS
sw_avg$Exp_Fit_Per_Plant <- ifelse(sw_avg$Condition == "mixed",
                                              Exp_Fit_Mixed(sw_avg$FITNESS),
                                              Exp_Single(sw_avg$FITNESS))

### Total Weight
sw_avg$Exp_TW_Per_Plant <- ifelse(sw_avg$Condition == "mixed",
                                             Exp_TW_mix(sw_avg$TOTAL_MASS),
                                             Exp_Single(sw_avg$TOTAL_MASS))
