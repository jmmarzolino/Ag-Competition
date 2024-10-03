#!/usr/bin/env Rscript

#SBATCH --ntasks=1
#SBATCH --mem=50G
#SBATCH --time=02:00:00
#SBATCH --job-name="Ag-Competition"
#SBATCH --output=/rhome/jmarz001/bigdata/Ag-Competition/000_single_v_mixed_PerSeedWeight_correlation.stdout
#SBATCH -p short


setwd("/rhome/jmarz001/bigdata/Ag-Competition/data")
library(tidyverse)

# load phenotyping data
df <- read_delim("Ag-Comp Seed Weights - Sheet3.csv")

# find average seed weight and flowering time per genotype-condition-Replicate ## wait that shouldn't average anything then, get rid of Replicate grouping?
df_mean <- df %>% group_by(Genotype, Condition, Replicate) %>% summarise(a100seed = mean(SEED_WEIGHT_100, na.rm=T), ft= mean(as.numeric(FT), na.rm=T))

# filter data to only single-genotype plots (control condition)
df_single <- df_mean %>% filter(Condition=="single") %>% pivot_wider(names_from=Replicate, values_from=c(a100seed, ft)) %>% filter(!is.na(`a100seed_rep 1`)) %>% filter(!is.na(`a100seed_rep 2`))

# calculate correlation across single-genotype plots Replicates
cor(df_single$`a100seed_rep 1`, df_single$`a100seed_rep 2`)
cor(df_single$`ft_rep 1`, df_single$`ft_rep 2`)


# filter data to mixed-genotype plots (test condition)
df_mixed <- df_mean %>% filter(Condition=="mixed") %>% pivot_wider(names_from=Replicate, values_from=c(a100seed, ft)) %>% filter(!is.na(`a100seed_rep 1`)) %>% filter(!is.na(`a100seed_rep 2`))

# calculate correlation across mixed-genotype plot Replicates
cor(df_mixed$`a100seed_rep 1`, df_mixed$`a100seed_rep 2`)

df_mixed <- df_mixed %>% filter(!is.na(`ft_rep 1`)) %>% filter(!is.na(`ft_rep 2`))
cor(df_mixed$`ft_rep 1`, df_mixed$`ft_rep 2`)












## Add a Generation Col
df$Generation <- 0
df[grep("^1$", df$FAM_ID), which(colnames(df)=='Generation')] <- 18
df[grep("^2$", df$FAM_ID), which(colnames(df)=='Generation')] <- 28
df[grep("^3$", df$FAM_ID), which(colnames(df)=='Generation')] <- 50
df[grep("^7$", df$FAM_ID), which(colnames(df)=='Generation')] <- 58

#### generation means and vars
df %>% select(c(Generation, Condition, SEED_WEIGHT_100)) %>% group_by(Generation, Condition) %>% summarise(avg_seed_weight = mean(SEED_WEIGHT_100, na.rm=T), seed_weight_var = var(SEED_WEIGHT_100, na.rm=T)) %>% pivot_wider(names_from=Condition, values_from=c(avg_seed_weight, seed_weight_var))

df$FT <- as.numeric(df$FT)
df %>% select(c(Generation, Condition, FT)) %>% group_by(Generation, Condition) %>% summarise(avg_ft = mean(FT, na.rm=T), ft_var = var(FT, na.rm=T)) %>% pivot_wider(names_from=Condition, values_from=c(avg_ft, ft_var))



#### parent genotypes avg flowering date
report <- df %>% filter(Generation == 0) %>% filter(Condition=='single') %>% group_by(Genotype) %>% summarise(mean = mean(SEED_WEIGHT_100, na.rm=T), variance = var(SEED_WEIGHT_100, na.rm=T), n=n())
write_delim(report, "avg_100seed_weight_parents.tsv", "\t")
