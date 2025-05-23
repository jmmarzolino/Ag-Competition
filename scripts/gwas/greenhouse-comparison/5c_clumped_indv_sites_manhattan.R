#!/usr/bin/env Rscript

#SBATCH --job-name=manhattan
#SBATCH -o /rhome/jmarz001/bigdata/Ag-Competition/scripts/gwas/greenhouse-comparison/5c_clumped_sites_manhattan.stdout
#SBATCH --ntasks=1
#SBATCH --mem=20gb
#SBATCH -t 00:30:00
#SBATCH -p short

setwd("/bigdata/koeniglab/jmarz001/Ag-Competition/results/gwas")
library(pacman)
p_load(tidyverse, data.table, ggsci, Cairo, qqman)

# make a list of the clumped-region gwas 
system("ls ASSOC_*.assoc.clumped > indv_clumped_trait_gwas.txt")
file_lst <- fread("indv_clumped_trait_gwas.txt", header=F)
lst <- file_lst$V1

addPlot <- function(FileName){
    # import data
    clump <- fread(unlist(FileName))
    clump2 <- clump[,1:11]
    clump2$NUM_CHR <- as.numeric(gsub("chr(\\d)H", "\\1", clump2$CHR))

    # cut assoc-file-name number for use in output file name
    i <- gsub("ASSOC_(.*).assoc.clumped", "\\1", FileName)
    # name output 
    OUTNAME <- paste0("indv_clumps_manhattan_", i, ".png")

    # print manhattan plot to output
    png(OUTNAME)
    manhattan(clump2, chr="NUM_CHR", bp="BP", p="P", col=c("dark blue", "dodger blue"))

    # close file
    dev.off()
}

# apply function to file list
lapply(lst, addPlot)
