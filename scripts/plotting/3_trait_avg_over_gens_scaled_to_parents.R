#!/usr/bin/env Rscript
#SBATCH --ntasks=1
#SBATCH --mem=30G
#SBATCH --time=02:00:00
#SBATCH --output=/rhome/jmarz001/bigdata/Ag-Competition/scripts/plotting/3_trait_avg_over_gens_scaled_to_parents.stdout
#SBATCH -p koeniglab

# This script makes dot/line plots for each trait 

library(tidyverse)
library(ggpubr)
library(RColorBrewer)
library(data.table)

setwd("/bigdata/koeniglab/jmarz001/Ag-Competition")
source("scripts/CUSTOM_FNS.R")

df <- read_delim("data/trait_BLUPs.tsv") 
df <- add_generation(df) %>% select(-Genotype)

# calculate parent generation mean and sd
parent_stats <- df %>%
  filter(Generation == 0) %>%
  select(-c(Generation)) %>%
  summarise(across(where(is.numeric),
                   list(mean = \(x) mean(x, na.rm=T), sd = \(y) sd(y, na.rm=T)),
                   .names = "{.col}_{.fn}"))

# list trait columns
trait_cols <- c("FT", "TOTAL_MASS", "SEED_WEIGHT_100", "MASS_PER_PLANT", "Germination", "FECUNDITY")

# z-scale the phenotypes with parent stats
parent_scaled <- df

for (t in trait_cols) {
  mean_val <- parent_stats[[paste0(t, "_mean")]]
  sd_val <- parent_stats[[paste0(t, "_sd")]]
  
  parent_scaled[[t]] <- (df[[t]] - mean_val) / sd_val
}


# the lowest color value is too light, so adjust the color scale down one
#display.brewer.pal(6, "Blues")
adjusted_blues <- brewer.pal(7, "Blues")[3:7]

## arrange data for facet plotting
df_long <- parent_scaled %>% pivot_longer(cols=-c('Generation'), values_to="VALUE", names_to="TRAIT")
# substitute trait names w/ tidy text versions
df_long$TRAIT <- tidy_text_substitution(df_long$TRAIT)

p_scale_avg <- parent_scaled %>% group_by(Generation) %>% summarise(across(where(is.numeric), \(x) mean(x, na.rm=T)))
p_scale_avg_lng <- p_scale_avg %>% pivot_longer(cols=-c('Generation'), values_to="VALUE", names_to="TRAIT")
p_scale_avg_lng$TRAIT <- tidy_text_substitution(p_scale_avg_lng$TRAIT)


a <- ggplot(df_long) +
      geom_point(aes(y=VALUE, x=Generation)) + 
      geom_line(data=p_scale_avg_lng, aes(x=Generation, y=VALUE)) +
      facet_wrap(~TRAIT) +
      theme_bw(base_size=20) +
      labs(x = "Generation", y = "") 
ggsave("results/trait_over_generations_parent_scaled.png", a, width = (7*3)+2, height = (7*2)+2)
