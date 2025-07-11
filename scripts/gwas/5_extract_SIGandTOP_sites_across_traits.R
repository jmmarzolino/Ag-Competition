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
files <- file_list$file

for(i in files){
    # read in file
    df <- fread(i)

    # extract the top x% of sites from gwas
    # define p-val based on quantiles to use as filter point
    p5 <- quantile(df$p_lrt, 0.05)
    p1 <- quantile(df$p_lrt, 0.01)
    p01 <- quantile(df$p_lrt, 0.001)
    bonn <- 0.05 / nrow(df)

    # keep only sites at or below our highest threshhold
    tmp <- df %>% filter(p_lrt <= p5)

    # label top 5% of associations
    tmp$top <- "5%"

    # for remaining sites, 
    # change top% col value if it meets higher threshhold
    # top 1% of assoc
    tmp[which(tmp$p_lrt<=p1), ncol(tmp)] <- "1%"

    # significant sites
    tmp[which(tmp$p_lrt<=bonn), ncol(tmp)] <- "Bonferroni"

    # top 0.1%
    tmp[which(tmp$p_lrt<=p01), ncol(tmp)] <- "0.1%"

    tmp <- tmp %>% select(c("chr", "ps", "beta", "p_lrt", "top"))
    # now save your resulting list of sites
    i2 <- gsub("(ASSOC_\\d+).assoc.txt", "\\1", i)
    write_delim(tmp, paste0(i2, "_top_sites.txt"))
}



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




#### stat checks
print(paste0("number of significant sites across traits and gwas models: ", length(df2$rs)))
print(paste0("number of significant sites across traits: ", length(unique(df2$rs))))

print("number of sig sites per trait")
print(table(df2$associated_trait))

print("number of significant sites for more than one trait (if any. traits w the same sig site twice indicate significance with 2 gwas models): ",)
table(df2$rs, df2$associated_trait)[which(rowSums(table(df2$rs, df2$associated_trait)) > 1), ]


print("are there any sites overlapping trait/method?")
print("differing numbers indicate site overlap")
tmp <- df2 %>% select(c(chr, rs, ps, associated_trait, gwas_method))
length(tmp$chr)
length(unique(tmp$rs))


print("are there any sites that are sig for multiple traits?")
print("if numbers match, there are no overlaping sites")
tmp2 <- df2 %>% select(c(chr, rs, ps, associated_trait)) %>% unique
length(tmp2$rs)
length(unique(tmp2$rs))
