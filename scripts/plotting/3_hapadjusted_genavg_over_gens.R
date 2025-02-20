#!/usr/bin/env Rscript

#SBATCH --ntasks=1
#SBATCH --mem=30G
#SBATCH --time=02:00:00
#SBATCH --output=/rhome/jmarz001/bigdata/Ag-Competition/scripts/plotting/3_hapadjusted_genavg_over_gens.stdout
#SBATCH -p koeniglab

library(tidyverse)

setwd("/bigdata/koeniglab/jmarz001/Ag-Competition/data")
source("../scripts/CUSTOM_FNS.R")
df <- read_delim("weighted_generation_trait_avgs.tsv")

## plot weighted trait averages over generations

# format data longways
df_lng <- df %>% pivot_longer(cols=colnames(df)[2:7], names_to="trait", values_to="value")
df_lng$trait <- tidy_text_substitution(df_lng$trait)

g <- ggplot(df_lng, aes(x= generation, y = value)) +
      geom_point() + geom_line() +
      facet_wrap(~trait, scales="free_y") + 
      labs(title = "Trait Average over Generations", subtitle="Haplotype-Frequency Weighted", x="Generation", y="") +
      theme_bw(base_size = 16) 
ggsave("../results/trait_avg_over_gens_hap_weighted.png", g, width = 16, height = 12)
