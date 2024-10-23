#!/usr/bin/env Rscript
#SBATCH -o /rhome/jmarz001/bigdata/Ag-Competition/scripts/gwas/2_phenotypes.stdout
#SBATCH --ntasks=1
#SBATCH --mem=64gb
#SBATCH -t 01:30:00
#SBATCH -p koeniglab

library(tidyverse)
library(data.table)

setwd("/rhome/jmarz001/bigdata/Ag-Competition/data")
# read in phenotypes file
pheno <- fread("FITNESS.tsv")

# remove parent lines that won't be in gwas
pheno <- pheno %>% filter(Generation != 0) 
# filter out 'mixed' condition plots
pheno <- pheno %>% filter(Condition != "mixed") %>% select(-Condition)

# record traits and corresponding col number
trait_names_list <- colnames(pheno)[3:ncol(pheno)]
trait_num_list <- 6:(length(trait_names_list)+5)
trait_n_list <- tibble(trait_names_list, trait_num_list)
write_delim(trait_n_list, "../results/trait_name_to_col_numbers.tsv")


## make genotypes match between genotype and phenotype files
plink_file <- fread("../results/all_traits.fam")
# remove 'no phenotype' tag from fam file
plink_file <- plink_file[,-6]

# find union of genotype numbers between the two data sets
# & keep only overlapping genotypes
joined <- inner_join(plink_file, pheno, by=c('V2'='Genotype'))
#left_join


#### over-write formatted genotype and phenotype .fam file
## format = famID / indvID / <phenotype / cols /...>
write_delim(joined, "../results/all_traits.fam", " ", col_names=F)
