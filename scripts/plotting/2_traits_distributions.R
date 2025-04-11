#!/usr/bin/env Rscript
#SBATCH --ntasks=1
#SBATCH --mem=30G
#SBATCH --time=02:00:00
#SBATCH --output=/rhome/jmarz001/bigdata/Ag-Competition/scripts/plotting/2_traits_distributions.stdout
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

# add column for parent or progeny marker
df$pgroup <- "Parents"
df[which(df$Generation != 0), which(colnames(df) == "pgroup")] <- "Progeny"

# the lowest color value is too light, so adjust the color scale down one
#display.brewer.pal(6, "Blues")
adjusted_blues <- brewer.pal(7, "Blues")[3:7]


### PLOTTING
## arrange data for facet plotting
df_long <- df %>% pivot_longer(cols=-c('Genotype', 'Generation', 'pgroup'), values_to="VALUE", names_to="TRAIT")

# substitute trait names w/ tidy text versions
df_long$TRAIT <- tidy_text_substitution(df_long$TRAIT)


# plot trait distributions
g <- ggplot(df_long, aes(VALUE)) +
  geom_density(linewidth=1) +
  facet_wrap(~TRAIT, scales="free") +  
  labs(x="", y="density") +
  theme_bw(base_size=20) +
  stat_summary(fun = median, geom = "vline", orientation = "y", aes(xintercept = after_stat(x), y = 0), color="dodgerblue3", linewidth=1) +
  stat_summary(fun = mean, geom = "vline", orientation = "y", aes(xintercept = after_stat(x), y = 0), color="dodgerblue", linewidth=1) 
ggsave("results/trait_distributions.png", g, width=12)


# one density line per generation
g <- ggplot(df_long, aes(VALUE, group=Generation, color=as.factor(Generation))) +
  geom_density(linewidth=1) +
  scale_color_manual(values=adjusted_blues, name="Generation") + 
  facet_wrap(~TRAIT, scales="free") +
  theme_bw(base_size=20) +
  labs(x="", y="density") +
  stat_summary(fun = mean, geom = "vline", orientation = "y", aes(xintercept = after_stat(x), y = 0), linewidth=1) 
ggsave("results/trait_distributions_Wgeneration.png", g, width=20)


# density lines for parents and progeny groups
g <- ggplot(df_long, aes(VALUE, group=pgroup, color=as.factor(pgroup))) +
  geom_density(linewidth=1) +
  scale_color_manual(values=adjusted_blues[c(1, 5)], name="Generation") + 
  facet_wrap(~TRAIT, scales="free") +
  theme_bw(base_size=20) +
  labs(x="", y="density") +
  stat_summary(fun = mean, geom = "vline", orientation = "y", aes(xintercept = after_stat(x), y = 0), linewidth=1) 
ggsave("results/trait_distributions_Wparentprogeny.png", g, width=20)

