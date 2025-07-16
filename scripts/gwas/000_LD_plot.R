#!/usr/bin/env Rscript
#SBATCH -o /rhome/jmarz001/bigdata/Ag-Competition/scripts/gwas/000_LD_plot.stdout
#SBATCH --ntasks=1
#SBATCH --mem=64gb
#SBATCH -t 01:30:00
#SBATCH -p short

library(tidyverse)
library(data.table)

setwd("/rhome/jmarz001/bigdata/Ag-Competition/results/gwas")
source("../../scripts/CUSTOM_FNS.R")

for(i in 1:7){

    # read in chr file 
    IN <- paste0("all_traits_chr", i, ".ld.gz")
    ld <- fread(IN)

    # print average LD over chromosome
    print(mean(ld$R2))
    print(quantile(ld$R2))

    # make heatmap plot
    tit <- paste0("CHR", i, " LD")
    g <- ggplot(ld, aes(x=SNP_A, y=SNP_B, fill=R2)) + geom_tile() + labs(title=tit)

    # write out
    ggsave(paste0("LD_chr", i, ".png"), g)

    # make distribution of LD values?
    ggplot(ld, aes(R2)) + geom_histogram()
}
