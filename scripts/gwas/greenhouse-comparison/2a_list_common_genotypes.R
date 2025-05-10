#!/usr/bin/env Rscript
#SBATCH -o /rhome/jmarz001/bigdata/Ag-Competition/scripts/gwas/greenhouse-comparison/2_joining_exp_phenotypes.stdout
#SBATCH --ntasks=1
#SBATCH --mem=24gb
#SBATCH -t 00:30:00
#SBATCH -p koeniglab

library(tidyverse)
library(data.table)

setwd("/rhome/jmarz001/bigdata/Ag-Competition/results/gwas/CCII_greenhouse_exp_gwas")
source("../../../scripts/CUSTOM_FNS.R")
## find the genotypes in common between Ag-Competition & CC II greenhouse experiments

# read in phenotypes file
fld <- fread("../all_traits.fam") %>% select(V1)

gh <- fread("CCII_Raw_phenotype_data.txt") %>% select(Line)
# remove parental lines
gh <- gh[-grep("^\\d+$", gh$Line), ]

# join the two trait data sets by their common genotypes
join <- inner_join(fld, gh, by=c("V1"="Line"))
# make genotype list 2 columns (plink format)
join$V2 <- join$V1
# write out genotype list in order
write_delim(join, "exp_common_genos", delim=" ", col_names=F)
