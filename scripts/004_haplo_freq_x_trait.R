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

## goal --
## HAPLO  FREQ  TRAIT DATA
## haplotype    parents_freq    f18_freq    f28_freq    f50_freq    f58_freq    FT  FEC     SW


hap <- fread("hap_assign.txt")
df <- fread("BLUPs.tsv")

# set haplotypes to factors instead of numbers
hap$Haplotype <- as.factor(hap$Haplotype)

# cut 3-digit genotype codes down to 2 digits
# to match family lines
df$Genotype <- gsub("(\\d+_\\d+)_\\d", "\\1", df$Genotype)

df <- add_generation(df)


# average phenotype for each haplotype
hap_join <- inner_join(hap, df, by = "Genotype")

# how common is each haplotype per generation
hap_freq_per_gen <- table(hap_join$Haplotype, hap_join$Generation)
hap_freq_per_gen_progeny <- tibble("Haplotype"=rownames(hap_freq_per_gen), "F18_freq"= hap_freq_per_gen[,1], "F28_freq"= hap_freq_per_gen[,2], "F50_freq"= hap_freq_per_gen[,3], "F58_freq"= hap_freq_per_gen[,4])

# average phenotype across haplotype
hap_trait_avg <- hap_join %>% group_by(Haplotype) %>% summarise(across(where(is.numeric), mean)) 

# join trait-value-per-haplotype w haplotype frequency per generation
haplo_freq_x_trait <- full_join(hap_freq_per_gen_progeny, hap_trait_avg) %>% select(-Generation)
fwrite(haplo_freq_x_trait, "haplo_freq_x_trait.tsv")
