#!/usr/bin/env Rscript
#SBATCH --mem=30G
#SBATCH --time=02:00:00
#SBATCH --output=/rhome/jmarz001/bigdata/Ag-Competition/scripts/potting/indv_trait_barplots.stdout
#SBATCH -p short

library(tidyverse)
library(ggpubr)
library(data.table)

setwd("/bigdata/koeniglab/jmarz001/Ag-Competition")
df <- fread("data/MASS_PER_PLANT.tsv")



collist <- c('FT', 'TOTAL_MASS', 'SEED_WEIGHT_100', 'Plants', 'SEED_COUNT', 'FECUNDITY',  'MASS_PER_PLANT', 'RELATIVE_MASS_PER_PLANT', 'AT_REL_MASS_PER_PLANT')

for(trait in collist) {

    g <- ggplot(df, aes(x=reorder(Genotype, get(trait)), y=get(trait), fill=as.factor(Generation))) +
            geom_col() +
            facet_wrap(~Condition) +
            theme_bw(base_size = 12) +
            theme(axis.text.x = element_text(angle = 90)) +
            labs(x="Genotype", y=trait, fill = "Generation")
    assign(paste0("g_", trait), g)

}

g_all <- ggarrange(g_FT, g_TOTAL_MASS, g_SEED_WEIGHT_100, g_Plants, g_SEED_COUNT, g_FECUNDITY, g_MASS_PER_PLANT, g_RELATIVE_MASS_PER_PLANT, g_AT_REL_MASS_PER_PLANT, ncol=3, nrow=3)

ggsave("results/barplot_trait_vals_per_genotype.png", g_all, height=20, width=20)
