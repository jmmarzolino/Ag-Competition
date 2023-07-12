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

df_mean <- df %>% group_by(Genotypes, Condition, replicate) %>% summarise(a100seed = mean(a100_seed_weight, na.rm=T), ft= mean(as.numeric(Flowering_Date), na.rm=T))

df_single <- df_mean %>% filter(Condition=="single") %>% pivot_wider(names_from=replicate, values_from=c(a100seed, ft)) %>% filter(!is.na(`a100seed_rep 1`)) %>% filter(!is.na(`a100seed_rep 2`))

cor(df_single$`a100seed_rep 1`, df_single$`a100seed_rep 2`)
cor(df_single$`ft_rep 1`, df_single$`ft_rep 2`)



df_mixed <- df_mean %>% filter(Condition=="mixed") %>% pivot_wider(names_from=replicate, values_from=c(a100seed, ft)) %>% filter(!is.na(`a100seed_rep 1`)) %>% filter(!is.na(`a100seed_rep 2`))

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
df %>% select(c(Generation, Condition, a100_seed_weight)) %>% group_by(Generation, Condition) %>% summarise(avg_seed_weight = mean(a100_seed_weight, na.rm=T), seed_weight_var = var(a100_seed_weight, na.rm=T)) %>% pivot_wider(names_from=Condition, values_from=c(avg_seed_weight, seed_weight_var))

df$Flowering_Date <- as.numeric(df$Flowering_Date)
df %>% select(c(Generation, Condition, Flowering_Date)) %>% group_by(Generation, Condition) %>% summarise(avg_ft = mean(Flowering_Date, na.rm=T), ft_var = var(Flowering_Date, na.rm=T)) %>% pivot_wider(names_from=Condition, values_from=c(avg_ft, ft_var))



#### parent genotypes avg flowering date
report <- df %>% filter(Generation == 0) %>% filter(Condition=='single') %>% group_by(Genotypes) %>% summarise(mean = mean(a100_seed_weight, na.rm=T), variance = var(a100_seed_weight, na.rm=T), n=n())
write_delim(report, "avg_100seed_weight_parents.tsv", "\t")
