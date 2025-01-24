#!/usr/bin/env Rscript
#SBATCH -o /rhome/jmarz001/bigdata/Ag-Competition/scripts/gwas/2_phenotypes.stdout
#SBATCH --ntasks=1
#SBATCH --mem=64gb
#SBATCH -t 01:30:00
#SBATCH -p koeniglab

library(tidyverse)
library(data.table)

setwd("/rhome/jmarz001/bigdata/Ag-Competition/results/gwas")
# read in phenotypes file
#pheno <- fread("../../data/FITNESS.tsv")
pheno <- fread("../../data/trait_BLUPs.tsv")

## filter trait cols with 0 variance
novar <- which(apply(pheno[2:ncol(pheno)], 2, var) == 0)
pheno <- pheno %>% select(-c(names(novar)))

# remove parent lines that won't be in gwas
#pheno <- pheno %>% filter(Generation != 0) %>% select(-Generation)
# filter out 'mixed' condition plots
#pheno <- pheno %>% filter(Condition != "mixed") %>% select(-Condition)

# scale the phenotypes
#pheno <- pheno %>% mutate(across(-c(Genotype), ~(scale(.) %>% as.vector)))


# record traits and corresponding col number
trait_names <- colnames(pheno)[2:ncol(pheno)]
trait_num <- 6:(length(trait_names)+5)
file <- paste0("ASSOC_", trait_num, ".assoc.txt")
trait_n <- tibble(trait_names, trait_num, file)
write_delim(trait_n, "trait_name_to_col_numbers.tsv")


## join phenotype and genotype file based on common genotypes
plink_file <- fread("all_traits.fam")
# remove any existing &/or 'no phenotype' tag from fam file
plink_file <- plink_file[,1:5]

# find union of genotype numbers between the two data sets
# & keep only overlapping genotypes
joined <- inner_join(plink_file, pheno, by=c('V2'='Genotype'))

### over-write formatted genotype and phenotype .fam file
## format = famID / indvID / <phenotype / cols /...>
write_delim(joined, "all_traits.fam", " ", col_names=F)
