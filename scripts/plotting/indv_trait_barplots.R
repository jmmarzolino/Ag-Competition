#!/usr/bin/env Rscript
#SBATCH --mem=30G
#SBATCH --time=02:00:00
#SBATCH --output=/rhome/jmarz001/bigdata/Ag-Competition/scripts/potting/indv_trait_barplots.stdout
#SBATCH -p short

library(tidyverse)
library(ggpubr)

setwd("/bigdata/koeniglab/jmarz001/Ag-Competition")
df <- fread("data/FITNESS.tsv")

pdf("results/barplot_trait_vals_per_genotype.pdf")


collist <- c('FT', 'TOTAL_MASS', 'SEED_WEIGHT_100', 'SURVIVAL', 'SEED_COUNT', 'FECUNDITY',  'FITNESS', 'RELATIVE_FITNESS', 'AT_REL_FITNESS')

for(trait in collist) {

    g <- ggplot(df, aes(x=reorder(Genotype, get(trait)), y=get(trait), fill=as.factor(Generation))) +
            geom_col() +
            facet_wrap(~Condition) +
            theme_bw() +
            labs(x="Genotype", y=trait, fill = "Generation")
    print(g)

}


dev.off()