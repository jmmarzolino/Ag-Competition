#!/usr/bin/env Rscript
#SBATCH -o /rhome/jmarz001/bigdata/Ag-Competition/scripts/gwas/2_phenotypes.stdout
#SBATCH --ntasks=1
#SBATCH --mem=64gb
#SBATCH -t 01:30:00
#SBATCH -p koeniglab

library(tidyverse)
library(data.table)

setwd("/rhome/jmarz001/bigdata/Ag-Competition/results/gwas")
source("../../scripts/CUSTOM_FNS.R")
# read in phenotypes file
pheno <- fread("../../data/trait_BLUPs.tsv")

## filter trait cols with 0 variance
#novar <- which(apply(pheno[2:ncol(pheno)], 2, var) == 0)
#pheno <- pheno %>% select(-c(names(novar)))

# shorten genotype codes to 2 generations (\\d_\\d), for consistency, IDing genotype families instead of lines
pheno$Genotype <- gsub("(\\d+_\\d+)_\\d+", "\\1", pheno$Genotype)
# remove parent lines that won't be in gwas
pheno <- add_generation(pheno)
pheno <- pheno %>% filter(Generation != 0) %>% select(-Generation)

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
joined <- left_join(plink_file, pheno, by=c('V2'='Genotype'))

# are fam & vcf genotypes in the same order?
#plink_file$V1 == pheno$Genotype
#plink_file$V1 == joined$V1
#plink_file$V1 == joined$V2
#joined[which(!(joined$V2 %in% plink_file$V2))]

### over-write formatted genotype and phenotype .fam file
## format = famID / indvID / <phenotype / cols /...>
write_delim(joined, "all_traits.fam", " ", col_names=F)
