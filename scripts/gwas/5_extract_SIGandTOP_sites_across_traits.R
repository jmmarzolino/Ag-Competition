#!/usr/bin/env Rscript
#SBATCH -o /rhome/jmarz001/bigdata/Ag-Competition/scripts/plotting/5_extract_SIGandTOP_sites_across_traits.stdout
#SBATCH --ntasks=1
#SBATCH --mem=64gb
#SBATCH -t 01:30:00
#SBATCH -p koeniglab

library(tidyverse)
library(data.table)

setwd("/rhome/jmarz001/bigdata/Ag-Competition/results/gwas")
source("../../scripts/CUSTOM_FNS.R")

file_list <- fread("association_files_traits.txt")
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

lmm_files$gwas_method <- "lmm"

lmm_files$associated_trait <- tidy_text_substitution(lmm_files$associated_trait)

write_delim(lmm_files, "all_gwas_sig_sites.tsv")




#### stat checks
print(paste0("number of significant sites across traits and gwas models: ", length(lmm_files$rs)))
print(paste0("number of significant sites across traits: ", length(unique(lmm_files$rs))))

print("number of sig sites per trait")
print(table(lmm_files$associated_trait))

print("number of significant sites for more than one trait (if any. traits w the same sig site twice indicate significance with 2 gwas models): ",)
table(lmm_files$rs, lmm_files$associated_trait)[which(rowSums(table(lmm_files$rs, lmm_files$associated_trait)) > 1), ]


print("are there any sites overlapping trait/method?")
print("differing numbers indicate site overlap")
tmp <- lmm_files %>% select(c(chr, rs, ps, associated_trait, gwas_method))
length(tmp$chr)
length(unique(tmp$rs))


print("are there any sites that are sig for multiple traits?")
print("if numbers match, there are no overlaping sites")
tmp2 <- lmm_files %>% select(c(chr, rs, ps, associated_trait)) %>% unique
length(tmp2$rs)
length(unique(tmp2$rs))
