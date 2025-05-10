#!/usr/bin/env Rscript
#SBATCH -o /rhome/jmarz001/bigdata/Ag-Competition/scripts/gwas/2_joining_exp_phenotypes.stdout
#SBATCH --ntasks=1
#SBATCH --mem=24gb
#SBATCH -t 00:30:00
#SBATCH -p koeniglab

library(tidyverse)
library(data.table)

setwd("/rhome/jmarz001/bigdata/Ag-Competition/results/gwas")
source("../../scripts/CUSTOM_FNS.R")
## find the genotypes in common between Ag-Competition & CC II greenhouse experiments

# read in phenotypes file
fld <- fread("all_traits.fam")

gh <- fread("CCII_greenhouse_exp_gwas/CCII_Raw_phenotype_data.txt")
# remove parental lines
gh <- gh[-grep("^\\d+$", gh$Line), ]

# join the two trait data sets by their common genotypes
join <- inner_join(fld, gh, by=c("V1"="Line"))
write_delim(join, "CCII_greenhouse_exp_gwas/common_exp.fam")

# record traits and corresponding col number
#trait_n <- read_delim("trait_name_to_col_numbers.tsv")
trait_n2 <- tibble("trait_names"=colnames(gh))
trait_n2$trait_num <- (1:nrow(trait_n2))+5
trait_n2$file <- paste0("ASSOC_", trait_n2$trait_num, ".assoc.txt")
write_delim(trait_n2, "CCII_greenhouse_exp_gwas/CCII_GH_trait_file_nums.tsv")

# write out genotype list in order
#exp_common_genos <- intersect(gh$Line, fld$V1)
genos <- join %>% select(V1)
write_delim(genos, "CCII_greenhouse_exp_gwas/exp_common_genos", col_names=F)