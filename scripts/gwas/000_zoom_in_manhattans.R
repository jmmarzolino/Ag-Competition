#!/usr/bin/env Rscript
#SBATCH -o /rhome/jmarz001/bigdata/Ag-Competition/scripts/gwas/000_zoom_in_manhattans.stdout
#SBATCH --ntasks=1
#SBATCH --mem=64gb
#SBATCH -t 01:30:00
#SBATCH -p short

library(tidyverse)
library(data.table)

setwd("/rhome/jmarz001/bigdata/Ag-Competition/results/gwas")
source("../../scripts/CUSTOM_FNS.R")
# read in phenotypes file
pheno <- fread("")

for(i in 1:5){
  df <- fread(paste0("plink_clump_lmm_r0.", i, ".clumped"))
  highs <- df$SNP

  png(paste0("clump_test_r", i, ".png"))

  manhattan(ft0, chr="chr", bp="ps", snp="rs", p="p_lrt", highlight=highs)

  dev.off()
}


for(i in 1:6){
  df <- fread(paste0("plink_clump_lmm_", i, "kb.clumped"))
  highs <- df$SNP

  png(paste0("clump_test_", i, ".png"))

  manhattan(ft0, chr="chr", bp="ps", snp="rs", p="p_lrt", highlight=highs)

  dev.off()
}





## 
manhattan(subset(ft0, chr==5), chr="chr", bp="ps", snp="rs", p="p_lrt", highlight=ft2$SNP)
manhattan(subset(ft0, chr==5), chr="chr", bp="ps", snp="rs", p="p_lrt", highlight=ft2$SNP, xlim=c(450000000, 550000000))

manhattan(subset(ft0, chr==4), chr="chr", bp="ps", snp="rs", p="p_lrt", highlight=ft2$SNP, xlim=c(588000000, 592000000))

annotatePval = 0.005, annotateTop = FALSE