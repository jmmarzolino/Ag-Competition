#!/usr/bin/env Rscript
#SBATCH --ntasks=1
#SBATCH --mem=30G
#SBATCH --time=02:00:00
#SBATCH --output=/rhome/jmarz001/bigdata/Ag-Competition/scripts/plotting/3_trait_avg_and_var_over_gens.stdout
#SBATCH -p koeniglab

# This script makes dot/line plots for each trait 

library(tidyverse)
library(ggpubr)
library(RColorBrewer)
library(data.table)

setwd("/bigdata/koeniglab/jmarz001/Ag-Competition")
source("scripts/CUSTOM_FNS.R")

#df <- read_delim("data/trait_BLUPs.tsv")
df <- fread("data/generations_trait_avg_var.tsv")

# scale the phenotypes
#pheno <- pheno %>% mutate(across(-c(Genotype), ~(scale(.) %>% as.vector)))

### PLOTTING
## arrange data for facet plotting
df_long <- df %>% pivot_longer(cols=-c('Generation'), values_to="VALUE", names_to="TRAIT")
# substitute trait names w/ tidy text versions
df_long$TRAIT <- tidy_text_substitution(df_long$TRAIT)

# split data into mean and var
df_long_mean <- df_long[grep("mean", df_long$TRAIT), ]
df_long_var <- df_long[grep("sd", df_long$TRAIT), ]

df_long_mean$TRAIT <- gsub("_mean", "", df_long_mean$TRAIT)
df_long_var$TRAIT <- gsub("_var", "", df_long_var$TRAIT)


# the lowest color value is too light, so adjust the color scale down one
#display.brewer.pal(6, "Blues")
adjusted_blues <- brewer.pal(7, "Blues")[3:7]


a <- ggplot(df_long_mean, aes(y=VALUE, x=Generation)) +
      geom_point() + geom_line() + 
      facet_wrap(~TRAIT, scales="free") +
      theme_bw(base_size=20) +
      labs(x = "Generation", y = "") 
ggsave("results/trait_avgs_over_generations.png", a, width = 14, height = 10)

a <- ggplot(df_long_var, aes(y=VALUE, x=Generation)) +
      geom_point() + geom_line() + 
      facet_wrap(~TRAIT, scales="free") +
      theme_bw(base_size=20) +
      labs(x = "Generation", y = "") 
ggsave("results/trait_var_over_generations.png", a, width = 14, height = 10)

