#!/usr/bin/env Rscript
#SBATCH -o /rhome/jmarz001/bigdata/Ag-Competition/scripts/gwas/greenhouse-comparison/2b_joining_exp_phenotypes.stdout
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
fld <- fread("../all_traits.fam")

gh <- fread("CCII_Raw_phenotype_data.txt")
# remove parental lines
gh <- gh[-grep("^\\d+$", gh$Line), ]

# join the two trait data sets by their common genotypes
join <- inner_join(fld, gh, by=c("V1"="Line"))

# read in empty fam file & format
fam <- fread("gh.fam")
fam <- fam[,1:5]

# check for the same genotypes
ifelse(length(which(join$V1 != fam$V1)) != 0, "Problem with genotype match", "genotypes match")

# overwrite fam file w phenotype version
write_delim(join, "gh.fam", delim=" ", col_names=F)

# record traits and corresponding col number
#trait_n <- read_delim("trait_name_to_col_numbers.tsv")
trait_n2 <- tibble("trait_names"=colnames(gh)[2:ncol(gh)])
trait_n2$trait_num <- (1:nrow(trait_n2))+5
trait_n2$file <- paste0("ASSOC_", trait_n2$trait_num, ".assoc.txt")
trait_n2$file_lmm <- paste0("ASSOC_", trait_n2$trait_num, "_lmm", ".assoc.txt")
write_delim(trait_n2, "CCII_GH_trait_file_nums.tsv")