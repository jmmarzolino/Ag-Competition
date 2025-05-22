#!/usr/bin/env Rscript
#SBATCH -o /rhome/jmarz001/bigdata/Ag-Competition/scripts/plotting/4_allele_freqs_over_gens.stdout
#SBATCH --ntasks=1
#SBATCH --mem=64gb
#SBATCH -t 01:30:00
#SBATCH -p koeniglab

library(tidyverse)
library(data.table)

setwd("/rhome/jmarz001/bigdata/Ag-Competition/")
source("scripts/CUSTOM_FNS.R")

# read in file
df <- fread("")

# extract the top x% of sites from gwas
sites <- fread("")

"ASSOC_\\d+_top_sites.txt"