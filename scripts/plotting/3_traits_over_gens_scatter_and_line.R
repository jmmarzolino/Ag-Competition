#!/usr/bin/env Rscript
#SBATCH --ntasks=1
#SBATCH --mem=30G
#SBATCH --time=02:00:00
#SBATCH --output=/rhome/jmarz001/bigdata/Ag-Competition/scripts/plotting/3_traits_over_gens_scatter_and_line.stdout
#SBATCH -p koeniglab

# This script plots scatterplots with linear regressions for each trait in the single subpopulation, as well as distributions for these traits

library(tidyverse)
library(ggpubr)
library(ggplot2)
library(RColorBrewer)

setwd("/bigdata/koeniglab/jmarz001/Ag-Competition")
source("scripts/CUSTOM_FNS.R")
df <- read_delim("data/BLUPs.tsv")
df <- add_generation(df)

#comparing trait values over generations
## scatter plots 

# table for generation means
gen_men <- df %>%
            group_by(Generation) %>%
            summarise(across(where(is.numeric), \(x) mean(x, na.rm=T))) 

## strip plots for each trait
a <- ggplot() + 
      geom_jitter(data=df, aes(x=Generation, y=FT), alpha=0.7) + 
      geom_line(data=gen_men, aes(x=Generation, y=FT), color="darkblue", linewidth=1) + 
      theme_bw(base_size=16) +
      labs(x = "Generation", y = tidy_text_substitution("FT"), title="Flowering Time over Generations") 

c <- ggplot(df) + 
      geom_jitter(aes(x=Generation, y=SEED_WEIGHT_100), alpha=0.7) +
      geom_line(data=gen_men, aes(x=Generation, y=SEED_WEIGHT_100), color="darkblue", linewidth=1) + 
      theme_bw(base_size=16) +
      labs(x = "Generation", y = tidy_text_substitution("SEED_WEIGHT_100"), title="Mass per 100-seeds over Generations") 

f <- ggplot() + 
      geom_jitter(data=df, aes(x=Generation,y=FECUNDITY), alpha=0.7) + 
      geom_line(data=gen_men, aes(x=Generation, y=FECUNDITY), color="darkblue", linewidth=1) + 
      theme_bw(base_size=16) +
      labs(x = "Generation", y = tidy_text_substitution("FECUNDITY"), title="Genotype Fecundity over Generations") 

ggcombo <- ggarrange(a, c, f, ncol=3, nrow=1)
h <- (7 * 1) + 2
w <- (7.25 * 3) + 2

ggsave("results/traits_over_generations_scatterplots.png", ggcombo, width = w, height = h)
