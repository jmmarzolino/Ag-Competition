#!/usr/bin/env Rscript

#SBATCH --ntasks=1
#SBATCH --mem=30G
#SBATCH --time=02:00:00
#SBATCH --output=/rhome/jmarz001/bigdata/Ag-Competition/scripts/plotting/3_hapadjusted_traits_over_gens_scatter_and_line_and_box.stdout
#SBATCH -p koeniglab

library(tidyverse)
library(ggpubr)
library(RColorBrewer)

setwd("/bigdata/koeniglab/jmarz001/Ag-Competition/data")
source("../scripts/CUSTOM_FNS.R")
df <- read_delim("trait_BLUPs_HapRepPop.tsv")

### plot scatterplots with average line for each trait
# the lowest color value is too light, so adjust the color scale down one
#display.brewer.pal(6, "Blues")
adjusted_blues <- brewer.pal(7, "Blues")[3:7]

# calculate generation means to use as gen-avg line points
gen_men <- df %>%
            group_by(Generation) %>%
            summarise(across(where(is.numeric), \(x) mean(x, na.rm=T))) 

### make strip plots
a <- ggplot() + 
      geom_jitter(data=df, aes(x=Generation, group=as.factor(Generation), y=FT)) + 
      geom_line(data=gen_men, aes(x=Generation, y=FT), color="blue", linewidth=1) + 
      theme_bw(base_size=18) +
      labs(x = "Generation", y = tidy_text_substitution("FT"), title="Flowering Time over Generations") 

c <- ggplot() + 
      geom_jitter(data=df, aes(x=Generation, group=as.factor(Generation), y=SEED_WEIGHT_100)) + 
      geom_line(data=gen_men, aes(x=Generation, y=SEED_WEIGHT_100), color="blue", linewidth=1) + 
      theme_bw(base_size=18) +
      labs(x = "Generation", y = tidy_text_substitution("SEED_WEIGHT_100"), title="100-Seed Weight over Generations") 

f <- ggplot() + 
      geom_jitter(data=df, aes(x=Generation, group=as.factor(Generation), y=FECUNDITY)) + 
      geom_line(data=gen_men, aes(x=Generation, y=FECUNDITY), color="blue", linewidth=1) + 
      theme_bw(base_size=18) +
      labs(x = "Generation", y = tidy_text_substitution("FECUNDITY"), title="Fecundity over Generations") 

# combine and save plots
ggcombo <- ggarrange(a, c, f, nrow=1)
ggsave("../results/traits_over_generations_scatterplots_hap_weighted.png", ggcombo, width = (7*3), height = 7)
