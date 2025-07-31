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
# set up phenotypes in fam file for site in common between field and greenhouse experiments

# read in phenotype data
gh <- fread("CCII_Raw_phenotype_data.txt")
# select relevant columns
gh <- gh %>% select(c(Line, PlantHeight, Mass_100, Total_mass, Seed_Estimate, Flowering_days_2017, Flowering_2018_Median))
# which flowering time year has more data?
colSums(is.na(gh[,6:7]))

gh <- gh %>% group_by(Line) %>% mutate("avg_FT"=mean(c(Flowering_days_2017, Flowering_2018_Median), na.rm=T)) 


# remove parental lines
gh <- gh[-grep("^\\d+$", gh$Line), ]

# read in empty fam file & format
fam <- fread("gh.fam")
fam <- fam[,1:5]

# attach phenotypes to fam file by joining the trait data by common genotypes
join <- inner_join(fam, gh, by=c("V1"="Line"))

# overwrite fam file w phenotype version
write_delim(join, "gh.fam", delim=" ", col_names=F)

# record traits and corresponding col number
#trait_n <- read_delim("trait_name_to_col_numbers.tsv")
trait_n2 <- tibble("trait_names"=colnames(gh)[2:ncol(gh)])
trait_n2$trait_num <- (1:nrow(trait_n2))+5
trait_n2$file_lmm <- paste0("ASSOC_", trait_n2$trait_num, "_lmm", ".assoc.txt")
write_delim(trait_n2, "CCII_GH_trait_file_nums.tsv")