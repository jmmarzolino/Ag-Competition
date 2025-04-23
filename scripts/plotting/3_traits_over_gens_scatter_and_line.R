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
df <- read_delim("data/trait_BLUPs.tsv")
df <- add_generation(df)

#comparing trait values over generations
## scatter- & box- plots 


# table for generation means
gen_men <- df %>%
            group_by(Generation) %>%
            summarise(across(where(is.numeric), \(x) mean(x, na.rm=T))) 


## strip plots for each trait
a <- ggplot() + 
      geom_jitter(data=df, aes(x=Generation, y=FT), alpha=0.7) + 
      geom_line(data=gen_men, aes(x=Generation, y=FT), color="dodgerblue3", linewidth=1) + 
      theme_bw(base_size=16) +
      labs(x = "Generation", y = tidy_text_substitution("FT"), title="Flowering Time over Generations") 

b <- ggplot(df) + 
      geom_jitter(aes(x=Generation, y=TOTAL_MASS), alpha=0.7) +
      geom_line(data=gen_men, aes(x=Generation, y=TOTAL_MASS), color="dodgerblue3", linewidth=1) + 
      theme_bw(base_size=16) +
      labs(x = "Generation", y = tidy_text_substitution("TOTAL_MASS"), title="Total Seed Weight over Generations") 

c <- ggplot() + 
      geom_jitter(data=df, aes(x=Generation, y=Germination), alpha=0.7) + 
      geom_line(data=gen_men, aes(x=Generation, y=Germination), color="dodgerblue3", linewidth=1) + 
      theme_bw(base_size=16) +
      labs(x = "Generation", y = tidy_text_substitution("Germination"), title="Seed Germination over Generations") 

d <- ggplot() + 
      geom_jitter(data=df, aes(x=Generation,y=FECUNDITY), alpha=0.7) + 
      geom_line(data=gen_men, aes(x=Generation, y=FECUNDITY), color="dodgerblue3", linewidth=1) + 
      theme_bw(base_size=16) +
      labs(x = "Generation", y = tidy_text_substitution("FECUNDITY"), title="Genotype Fecundity over Generations") 

e <- ggplot() + 
      geom_jitter(data=df, aes(x=Generation, y=MASS_PER_PLANT), alpha=0.7) + 
      geom_line(data=gen_men, aes(x=Generation, y=MASS_PER_PLANT), color="dodgerblue3", linewidth=1) + 
      theme_bw(base_size=16) +
      labs(x = "Generation", y = tidy_text_substitution("MASS_PER_PLANT"), title="Genotype Mass per Plant over Generations") 


ggcombo <- ggarrange(a, b, c, d, e, ncol=3, nrow=2)
h <- (7 * 2) + 2
w <- (7 * 3) + 2

ggsave("results/traits_over_generations_scatterplots.png", ggcombo, width = w, height = h)