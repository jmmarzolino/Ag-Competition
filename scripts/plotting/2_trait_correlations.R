#!/usr/bin/env Rscript
#SBATCH --ntasks=1
#SBATCH --mem=30G
#SBATCH --time=02:00:00
#SBATCH --output=/rhome/jmarz001/bigdata/Ag-Competition/scripts/2_trait_correlations.stdout
#SBATCH -p koeniglab

library(tidyverse)
library(ggpubr)
library(data.table)
library(corrplot)
#library(lme4)
#library(car)

setwd("/rhome/jmarz001/bigdata/Ag-Competition")
source("scripts/CUSTOM_FNS.R")

df <- fread("data/FITNESS.tsv")
df <- df %>% filter(Condition == "single") %>% select(-Condition)



# check correlations between traits 
# remove highly correlated phenotypes from gwas
# fitness / atlas-fitness / 
traits_df <- df %>% select(-c(Genotype, Generation)) 
x <- cor(traits_df, use="na.or.complete", method="spearman")
png("data/trait_correlations.png")
corrplot(x, method="color", type="upper", order="original", title="", mar=c(0,0,4,0), addCoef.col = "black")
dev.off()

png("data/trait_correlations_filtered.png")
traits_df <- df %>% select(-c(Genotype, Generation, SEED_COUNT, RELATIVE_FITNESS, AT_REL_FITNESS)) 
#traits_df <- pheno %>% select(c(ends_with("_scaled"))) %>% select(-c(SEED_COUNT_scaled, RELATIVE_FITNESS_scaled, AT_REL_FITNESS_scaled)) 
x <- cor(traits_df, use="na.or.complete", method="spearman")
corrplot(x, method="color", type="upper", order="original", title="", mar=c(0,0,4,0), addCoef.col = "black")
dev.off()


pdf("data/trait_correlations_per_generation.pdf")
# check correlations between traits for each generation
for(i in c(0, 18, 28, 50, 58)) {

  traits_df <- df %>% filter(Generation == i) %>% select(-c(Genotype, Generation, SEED_COUNT, RELATIVE_FITNESS, AT_REL_FITNESS)) 
  x <- cor(traits_df, use="na.or.complete", method="spearman")
  corrplot(x, method="color", type="upper", order="original", title=paste0("Generation ", i), mar=c(0,0,4,0), addCoef.col = "black")

}
dev.off()



## plot trait relationships w scatterplots
png("all_trait_correlations.png", height=40, width=40, units="in")
plot(df, col=as.factor(df$Generation))
dev.off()


# to see more clearly if the relationship between traits are changing over generations
# plot trait v. trait, colored & approx line by generation
# and look to see if value range or relationship change much...
tmp <- df %>% select(-c(Genotype, SURVIVAL, RELATIVE_FITNESS, SEED_COUNT, AT_REL_FITNESS)) %>% pivot_longer(-c(FT, Generation), values_to = "VALUE", names_to="TRAIT")

ggplot(tmp, aes(x=FT, y=VALUE, group=Generation)) +
    geom_point(tmp, aes(color=Generation)) +
    geom_abline(aes()) +
    facet_wrap(~TRAIT)