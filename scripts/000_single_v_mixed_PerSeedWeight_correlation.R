#!/usr/bin/env Rscript

#SBATCH --ntasks=1
#SBATCH --mem=50G
#SBATCH --time=02:00:00
#SBATCH --job-name="Ag-Competition"
#SBATCH --output=/rhome/jmarz001/bigdata/Ag-Competition/000_single_v_mixed_PerSeedWeight_correlation.stdout
#SBATCH -p short


setwd("/rhome/jmarz001/bigdata/Ag-Competition")
library(tidyverse)

# load phenotyping data
df <- read_delim("Ag-Comp Seed Weights - Sheet3.tsv")



#### generation means and vars
df3 %>% group_by(Generation) %>% summarise(mean = mean(Flowering_Date, na.rm=T), variance = var(Flowering_Date, na.rm=T), n=n())
### conditions means and vars
df3 %>% group_by(Condition) %>% summarise(mean = mean(Flowering_Date, na.rm=T), variance = var(Flowering_Date, na.rm=T), n=n())
### replicates means and vars
df3 %>% group_by(replicate) %>% summarise(mean = mean(Flowering_Date, na.rm=T), variance = var(Flowering_Date, na.rm=T), n=n())



#### parent genotypes avg flowering date
report <- df3 %>% filter(Generation == 0) %>% filter(Condition=='single') %>% group_by(Genotypes) %>% summarise(mean = mean(Flowering_Date, na.rm=T), variance = var(Flowering_Date, na.rm=T), n=n())
write_delim(report, "avg_flowering_date_parents.tsv", "\t")
