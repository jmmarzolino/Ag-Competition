#!/usr/bin/env Rscript
#SBATCH -o /rhome/jmarz001/bigdata/Ag-Competition/scripts/plotting/5d_extract_leading_snps.stdout
#SBATCH --ntasks=1
#SBATCH --mem=64gb
#SBATCH -t 01:30:00
#SBATCH -p koeniglab

library(tidyverse)
library(data.table)

setwd("/rhome/jmarz001/bigdata/Ag-Competition/results/gwas")
source("../../scripts/CUSTOM_FNS.R")


file_list <- fread("association_files_traits.txt")
file_list <- file_list$file

# read in every association file

# list of paired association files

for(i in file_list){
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
