#!/usr/bin/env Rscript
#SBATCH -o /rhome/jmarz001/bigdata/Ag-Competition/scripts/plotting/5d_extract_sig_sites_across_traits.stdout
#SBATCH --ntasks=1
#SBATCH --mem=64gb
#SBATCH -t 01:00:00
#SBATCH -p koeniglab

library(tidyverse)
library(data.table)

setwd("/rhome/jmarz001/bigdata/Ag-Competition/results/gwas")
source("../../scripts/CUSTOM_FNS.R")


file_list <- fread("association_files_traits.txt")

for(i in 1:nrow(file_list)){

    # read in file
    df <- fread(file_list[[i, 3]])

    # extract the top x% of sites from gwas
    # define bonnferroni threshhold for significance with multiple testing
    bonn <- 0.05 / nrow(df)

    # keep only sites at or below our threshhold
    df <- df %>% filter(p_lrt <= bonn)

    # label associated trait
    df$associated_trait <- file_list[[i, 1]]

    # now save your resulting list of sites
    assign(gsub("(ASSOC_\\d+).assoc.txt", "\\1", file_list[[i, 3]]), df)
}



lm_files <- full_join(ASSOC_6, ASSOC_7)
lm_files <- full_join(lm_files, ASSOC_8)
lm_files <- full_join(lm_files, ASSOC_9)
lm_files <- full_join(lm_files, ASSOC_10)
lm_files <- full_join(lm_files, ASSOC_11)

lm_files$gwas_method <- "lm"


# repeat for lmm model files
# with file name ref col changed & assigned name change
for(i in 1:nrow(file_list)){

    # read in file
    df <- fread(file_list[[i, 4]])

    # extract the top x% of sites from gwas
    # define bonnferroni threshhold for significance with multiple testing
    bonn <- 0.05 / nrow(df)

    # keep only sites at or below our threshhold
    df <- df %>% filter(p_lrt <= bonn)

    # label associated trait
    df$associated_trait <- file_list[[i, 1]]

    # now save your resulting list of sites
    assign(gsub("(ASSOC_\\d+_lmm).assoc.txt", "\\1", file_list[[i, 4]]), df)
}

lmm_files <- full_join(ASSOC_6_lmm, ASSOC_7_lmm)
lmm_files <- full_join(lmm_files, ASSOC_8_lmm)
lmm_files <- full_join(lmm_files, ASSOC_9_lmm)
lmm_files <- full_join(lmm_files, ASSOC_10_lmm)
lmm_files <- full_join(lmm_files, ASSOC_11_lmm)

lmm_files$gwas_method <- "lmm"


df2 <- full_join(lm_files, lmm_files)
df2$associated_trait <- tidy_text_substitution(df2$associated_trait)

write_delim(df2, "all_gwas_sig_sites.tsv")