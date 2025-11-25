#!/usr/bin/env Rscript
#SBATCH --mem=10G
#SBATCH --time=02:00:00
#SBATCH --output=/rhome/jmarz001/bigdata/Ag-Competition/scripts/000_field_bed_and_row_boxplots.stdout
#SBATCH -p koeniglab

setwd("/rhome/jmarz001/bigdata/Ag-Competition/data")
library(tidyverse)
library(data.table)
library(ggpubr)

## plot phenotypes per field bed/row for both years & replicates
df <- fread("JOINED_PHENOTYPES.tsv") 
df <- df %>% select(c(Genotype, Condition, Replicate, BED, ROW, Exp_year, Germination, FT, TOTAL_MASS, SEED_WEIGHT_100))

tmp <- df %>% pivot_longer(cols=c("Germination", "FT", "TOTAL_MASS", "SEED_WEIGHT_100"), names_to="trait")

# plot phenotype per bed
gb <- ggplot(tmp, aes(x=BED, y=value, group=BED, color=as.factor(Replicate))) +
      geom_boxplot() +
      facet_wrap(Exp_year ~ trait, nrow=2, scales="free") +
      theme_minimal() +
      labs(title="Beds", color="Replicate")
ggsave("pheno_x_beds.png", gb, height=(7*2), width=(7*4), units="in")

# plot phenotype per row
tmp1 <- tmp %>% filter(Replicate == 1)
gr1 <- ggplot(tmp1, aes(x=ROW, y=value, group=ROW)) +
      geom_boxplot() +
      facet_wrap(Exp_year ~ trait, nrow=2, scales="free") +
      theme_minimal() +
      labs(title="Rows Replicate 1")
ggsave("pheno_x_rows1.png", gr1, height=(7*2), width=(7*4), units="in")

tmp2 <- tmp %>% filter(Replicate == 2)
gr2 <- ggplot(tmp2, aes(x=ROW, y=value, group=ROW)) +
      geom_boxplot() +
      facet_wrap(Exp_year ~ trait, nrow=2, scales="free") +
      theme_minimal() +
      labs(title="Rows Replicate 2")
ggsave("pheno_x_rows2.png", gr2, height=(7*2), width=(7*4), units="in")

super <- ggarrange(gb, gr1, gr2, nrow=3)
ggsave("pheno_x_beds_and_rows.png", super, height=(7*6), width=(7*4), units="in", dpi=300)