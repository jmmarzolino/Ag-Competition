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
            summarise(across(where(is.numeric), mean)) 

### strip plots
a <- ggplot() + 
      geom_jitter(data=df, aes(x=Generation, group=as.factor(Generation), y=FT)) + 
      geom_line(data=gen_men, aes(x=Generation, y=FT), color="blue", linewidth=1) + 
      theme_bw(base_size=18) +
      labs(x = "Generation", y = tidy_text_substitution("FT"), title="Flowering Time over Generations") #, subtitle="Line connects generation averages") 

# scale total mass by a factor of 10^10
tmp <- df %>% mutate('TOTAL_MASS_scaled' = (TOTAL_MASS + 0.5484163055) * 10^11)
mean_tmp <- tmp %>% group_by(Generation) %>% summarise('totmas' = mean(TOTAL_MASS_scaled)) 
b <- ggplot(tmp) + 
      geom_jitter(aes(x=Generation, y=TOTAL_MASS_scaled)) +
      geom_line(data=mean_tmp, aes(x=Generation, y=totmas), color="blue", linewidth=1) + 
      theme_bw(base_size=18) +
      labs(x = "Generation", y = tidy_text_substitution("TOTAL_MASS"), title="Total Seed Weight over Generations") 

c <- ggplot() + 
      geom_jitter(data=df, aes(x=Generation, group=as.factor(Generation), y=Plants)) + 
      geom_line(data=gen_men, aes(x=Generation, y=Plants), color="blue", linewidth=1) + 
      theme_bw(base_size=18) +
      labs(x = "Generation", y = tidy_text_substitution("Plants"), title="Seed Germination over Generations") 

d <- ggplot() + 
      geom_jitter(data=df, aes(x=Generation, group=as.factor(Generation), y=FECUNDITY)) + 
      geom_line(data=gen_men, aes(x=Generation, y=FECUNDITY), color="blue", linewidth=1) + 
      theme_bw(base_size=18) +
      labs(x = "Generation", y = tidy_text_substitution("FECUNDITY"), title="Genotype Fecundity over Generations") 

e <- ggplot() + 
      geom_jitter(data=df, aes(x=Generation, group=as.factor(Generation), y=MASS_PER_PLANT)) + 
      geom_line(data=gen_men, aes(x=Generation, y=MASS_PER_PLANT), color="blue", linewidth=1) + 
      theme_bw(base_size=18) +
      labs(x = "Generation", y = tidy_text_substitution("MASS_PER_PLANT"), title="Genotype Mass per Plant over Generations") 

ggcombo <- ggarrange(a, b, c, d, e)
ggsave("results/traits_over_generations_scatterplots.png", ggcombo, width = 20, height = 10)