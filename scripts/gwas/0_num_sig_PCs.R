#!/usr/bin/env Rscript

#SBATCH --ntasks=1
#SBATCH --mem=200G
#SBATCH --time=12:00:00
#SBATCH --output=/rhome/jmarz001/bigdata/Ag-Competition/scripts/gwas/0_num_sig_PCs.stdout
#SBATCH -p koeniglab

### PCA of trait slopes and slope means

#######  LOAD ENVIRONMENT  #######
setwd("/bigdata/koeniglab/jmarz001/Ag-Competition/results/gwas")

library(tidyverse)
library(data.table)
#devtools::install_github("arleyc/PCAtest")
library(PCAtest)


#df <- fread("../../data/FITNESS.tsv")
df <- fread("all_traits.fam")

# select only numeric columns
df2 <- df[,c(6:11)]
# check for missing values
apply(df, 2, is.na) %>% rowSums
apply(df, 2, is.na) %>% colSums

# determine number of significant PCs
pdf("PCAtest_result.pdf")
result <- PCAtest(
              x = df2,
              nperm = 1000,
              nboot = 1000,
              alpha = 0.05,
              indload = TRUE,
              varcorr = T,
              counter = TRUE,
              plot = TRUE
            )

dev.off()

#PC 1 is significant and accounts for 48.9% of the total variation
#Variables 2, 5, and 6 have significant loadings on PC 1
