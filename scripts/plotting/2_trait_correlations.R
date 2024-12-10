#!/usr/bin/env Rscript
#SBATCH --ntasks=1
#SBATCH --mem=30G
#SBATCH --time=02:00:00
#SBATCH --output=/rhome/jmarz001/bigdata/Ag-Competition/scripts/plotting/2_trait_correlations.stdout
#SBATCH -p koeniglab

library(tidyverse)
library(ggpubr)
library(RColorBrewer)
library(data.table)
library(corrplot)

setwd("/rhome/jmarz001/bigdata/Ag-Competition/data")
source("../scripts/CUSTOM_FNS.R")

df <- fread("DERIVED_PHENOTYPES.tsv")
df <- df %>% 
        filter(Condition == "single") %>% 
        select(-Condition) %>%
        group_by(Genotype) %>%
        summarise(across(where(is.numeric), mean)) %>%
        ungroup() %>%
        select(-c(Replicate, Genotype, Exp_year, SEED_COUNT))


# check correlations between traits 
# remove highly correlated phenotypes from gwas
# fitness / atlas-fitness / relative-fitness
traits_df <- df %>% select(-c(Generation)) 
x <- cor(traits_df, use="na.or.complete", method="spearman")
png("trait_correlations.png")
corrplot(x, method="color", type="upper", order="original", title="", mar=c(0,0,4,0), addCoef.col = "black")
dev.off()


pdf("trait_correlations_per_generation.pdf")
# check correlations between traits for each generation
for(i in c(0, 18, 28, 50, 58)) {

  traits_df <- df %>% filter(Generation == i) %>% select(-c(Generation)) 
  x <- cor(traits_df, use="na.or.complete", method="spearman")
  corrplot(x, method="color", type="upper", order="original", title=paste0("Generation ", i), mar=c(0,0,4,0), addCoef.col = "black")

}
dev.off()



## plot trait relationships w scatterplots
png("all_trait_correlations.png", height=40, width=40, units="in", res=300)
tmp <- df %>% select(-c(Generation))
plot(tmp)
dev.off()


# to see more clearly if the relationship between traits are changing over generations
# plot trait v. trait, colored & approx line by generation
# and look to see if value range or relationship change much...

# the lowest color value is too light, so adjust the color scale down one
#display.brewer.pal(6, "Blues")
adjusted_blues <- brewer.pal(7, "Blues")[3:7]

# flowering time plots
df_ft <- df %>% 
        pivot_longer(-c(FT, Generation), values_to = "VALUE", names_to="TRAIT")
df_ft$TRAIT <- tidy_text_substitution(df_ft$TRAIT)

gee <- df_ft %>% 
        ggplot(aes(x=FT, y=VALUE,)) + 
          geom_point(aes(color=as.factor(Generation)), alpha=0.7) + 
          geom_smooth() +
          scale_color_manual(values=adjusted_blues) + 
          facet_wrap(~TRAIT, scales="free_y") +
          labs(x="Flowering Time", y="", color="Generation") +
          theme_bw()
ggsave("traits_vs_FT.png", gee, width=14)

gee <- df_ft %>% 
        ggplot(aes(x=FT, y=VALUE, group=as.factor(Generation))) + 
          geom_point(aes(color=as.factor(Generation)), alpha=0.7) + 
          geom_smooth(aes(color=as.factor(Generation))) +
          scale_color_manual(values=adjusted_blues) + 
          facet_wrap(~TRAIT, scales="free_y") +
          labs(x="Flowering Time", y="", color="Generation") +
          theme_bw()
ggsave("traits_vs_FT_by_generation.png", gee, width=14)


# fitness plots
df_fit <- df %>% 
        pivot_longer(-c(FITNESS, Generation), values_to = "VALUE", names_to="TRAIT")
df_fit$TRAIT <- tidy_text_substitution(df_fit$TRAIT)

gee_wiz <- df_fit %>% 
        ggplot(aes(x=FITNESS, y=VALUE)) + 
          geom_point(aes(color=as.factor(Generation)), alpha=0.7) + 
          geom_smooth() +
          scale_color_manual(values=adjusted_blues) + 
          facet_wrap(~TRAIT, scales="free_y") +
          labs(x="Fitness", y="", color="Generation") +
          theme_bw()
ggsave("traits_vs_FIT.png", gee_wiz, width=14)


gee_wiz <- df_fit %>% 
        ggplot(aes(x=FITNESS, y=VALUE, group=Generation)) + 
          geom_point(aes(color=as.factor(Generation)), alpha=0.7) + 
          geom_smooth(aes(color=as.factor(Generation))) +
          scale_color_manual(values=adjusted_blues) + 
          facet_wrap(~TRAIT, scales="free_y") +
          labs(x="Fitness", y="", color="Generation") +
          theme_bw()
ggsave("traits_vs_FIT_by_generation.png", gee_wiz, width=14)