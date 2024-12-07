#!/usr/bin/env Rscript
#SBATCH --ntasks=1
#SBATCH --mem=30G
#SBATCH --time=02:00:00
#SBATCH --output=/rhome/jmarz001/bigdata/Ag-Competition/scripts/2_trait_correlations.stdout
#SBATCH -p koeniglab

library(tidyverse)
library(ggpubr)
library(RColorBrewer)
library(data.table)
library(corrplot)
#library(lme4)
#library(car)

setwd("/rhome/jmarz001/bigdata/Ag-Competition/data")

df <- fread("DERIVED_PHENOTYPES.tsv")
df <- df %>% 
        filter(Condition == "single") %>% 
        select(-Condition) %>%
        group_by(Genotype) %>%
        summarise(across(where(is.numeric), mean)) %>%
        ungroup() %>%
        select(-c(Replicate, Exp_year, SEED_COUNT))


blup <- fread("trait_BLUPs.tsv")


# check correlations between traits 
# remove highly correlated phenotypes from gwas
# fitness / atlas-fitness / relative-fitness
traits_df <- df %>% select(-c(Genotype, Generation)) 
x <- cor(traits_df, use="na.or.complete", method="spearman")
png("trait_correlations.png")
corrplot(x, method="color", type="upper", order="original", title="", mar=c(0,0,4,0), addCoef.col = "black")
dev.off()


pdf("trait_correlations_per_generation.pdf")
# check correlations between traits for each generation
for(i in c(0, 18, 28, 50, 58)) {

  traits_df <- df %>% filter(Generation == i) %>% select(-c(Genotype, Generation)) 
  x <- cor(traits_df, use="na.or.complete", method="spearman")
  corrplot(x, method="color", type="upper", order="original", title=paste0("Generation ", i), mar=c(0,0,4,0), addCoef.col = "black")

}
dev.off()



## plot trait relationships w scatterplots
png("all_trait_correlations.png", height=40, width=40, units="in", res=300)
tmp <- df %>% select(-c(Genotype, Generation))
plot(tmp)
dev.off()


# to see more clearly if the relationship between traits are changing over generations
# plot trait v. trait, colored & approx line by generation
# and look to see if value range or relationship change much...


# the lowest color value is too light, so adjust the color scale down one
#display.brewer.pal(6, "Blues")
adjusted_blues <- brewer.pal(7, "Blues")[3:7]

gee <- tmp %>% 
        pivot_longer(-c(FT, Generation), values_to = "VALUE", names_to="TRAIT") %>% 
        ggplot(tmp, aes(x=FT, y=VALUE, group=Generation)) + 
          geom_point(aes(color=as.factor(Generation)), alpha=0.7) + 
          geom_smooth(aes(color=as.factor(Generation))) +
          scale_color_manual(values=adjusted_blues) + 
          facet_wrap(~TRAIT, scales="free_y") +
          labs(x="Flowering Time", y="", color="Generation") +
          theme_bw()
ggsave("traits_vs_FT.png", gee)



gee_wiz <- tmp %>% 
        pivot_longer(-c(FIT, Generation), values_to = "VALUE", names_to="TRAIT") %>% 
        ggplot(tmp, aes(x=FIT, y=VALUE, group=Generation)) + 
          geom_point(aes(color=as.factor(Generation)), alpha=0.7) + 
          geom_smooth(aes(color=as.factor(Generation))) +
          scale_color_manual(values=adjusted_blues) + 
          facet_wrap(~TRAIT, scales="free_y") +
          labs(x="Flowering Time", y="", color="Generation") +
          theme_bw()
ggsave("traits_vs_FIT.png", gee_wiz)