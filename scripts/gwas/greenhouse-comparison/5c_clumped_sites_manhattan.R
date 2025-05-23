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

clump <- fread("plink_clump.clumped")
clump2 <- clump[,1:11]
clump2$NUM_CHR <- as.numeric(gsub("chr(\\d)H", "\\1", clump2$CHR))

png("all_clumps_manhattan.png")
manhattan(clump2, chr="NUM_CHR", bp="BP", p="P", col=c("dark blue", "dodger blue"))
dev.off()



########## lmm model results
clump <- fread("plink_clump_lmm.clumped")
clump2 <- clump[,1:11]
clump2$NUM_CHR <- as.numeric(gsub("chr(\\d)H", "\\1", clump2$CHR))

png("all_clumps_manhattan_lmm.png")
manhattan(clump2, chr="NUM_CHR", bp="BP", p="P", col=c("dark blue", "dodger blue"))
dev.off()

#suggestiveline: Where to draw a "suggestive" line. Default -log10(1e-5). Set to FALSE to disable.

#genomewideline: Where to draw a "genome-wide sigificant" line. Default -log10(5e-8). Set to FALSE to disable.

#highlight: A character vector of SNPs in your dataset to highlight.
#          These SNPs should all be in your dataset.
