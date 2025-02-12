#!/usr/bin/env Rscript
#SBATCH --ntasks=1
#SBATCH --mem=30G
#SBATCH --time=02:00:00
#SBATCH --output=/rhome/jmarz001/bigdata/Ag-Competition/scripts/plotting/2_blup_correlations.stdout
#SBATCH -p koeniglab

library(tidyverse)
library(ggpubr)
library(RColorBrewer)
library(data.table)
library(corrplot)

setwd("/rhome/jmarz001/bigdata/Ag-Competition/data")
source("../scripts/CUSTOM_FNS.R")

df <- fread("trait_BLUPs.tsv")
df <- add_generation(df) %>% select(-c(Genotype))


# check correlations between traits 
# remove highly correlated phenotypes from gwas
# fitness / atlas-fitness / relative-fitness
traits_df <- df %>% select(-c(Generation)) 
x <- cor(traits_df, use="na.or.complete", method="spearman")
png("blup_correlations.png")
corrplot(x, method="color", type="upper", order="original", title="", mar=c(0,0,4,0), addCoef.col = "black")
dev.off()


pdf("blup_correlations_per_generation.pdf")
# check correlations between traits for each generation
for(i in c(0, 18, 28, 50, 58)) {

  traits_df <- df %>% filter(Generation == i) %>% select(-c(Generation)) 
  x <- cor(traits_df, use="na.or.complete", method="spearman")
  corrplot(x, method="color", type="upper", order="original", title=paste0("Generation ", i), mar=c(0,0,4,0), addCoef.col = "black")

}
dev.off()



## plot trait relationships w scatterplots
png("all_blup_correlations.png", height=40, width=40, units="in", res=300)
plot(traits_df)
dev.off()


# to see more clearly if the relationship between traits are changing over generations
# plot trait v. trait, colored & approx line by generation
# and look to see if value range or relationship change much...

# the lowest color value is too light, so adjust the color scale down one
#display.brewer.pal(6, "Blues")
adjusted_blues <- brewer.pal(7, "Blues")[3:7]

# flowering time plots
df_ft <- df %>% 
        pivot_longer(-c(FT_blup, Generation), values_to = "VALUE", names_to="TRAIT")
df_ft$TRAIT <- tidy_text_substitution(df_ft$TRAIT)

gee <- df_ft %>% 
        ggplot(aes(x=FT_blup, y=VALUE,)) + 
          geom_point(aes(color=as.factor(Generation)), alpha=0.7) + 
          geom_smooth() +
          scale_color_manual(values=adjusted_blues) + 
          facet_wrap(~TRAIT, scales="free_y") +
          labs(x="Flowering Time BLUP", y="", color="Generation") +
          theme_bw()
ggsave("blups_vs_FT.png", gee, width=14)

gee <- df_ft %>% 
        ggplot(aes(x=FT_blup, y=VALUE, group=as.factor(Generation))) + 
          geom_point(aes(color=as.factor(Generation)), alpha=0.7) + 
          geom_smooth(aes(color=as.factor(Generation))) +
          scale_color_manual(values=adjusted_blues) + 
          facet_wrap(~TRAIT, scales="free_y") +
          labs(x="Flowering Time BLUP", y="", color="Generation") +
          theme_bw()
ggsave("blups_vs_FT_by_generation.png", gee, width=14)


# mass per plant plots
df_massper <- df %>% 
        pivot_longer(-c(MASS_PER_PLANT_blup, Generation), values_to = "VALUE", names_to="TRAIT")
df_massper$TRAIT <- tidy_text_substitution(df_massper$TRAIT)

gee_wiz <- df_massper %>% 
        ggplot(aes(x=MASS_PER_PLANT_blup, y=VALUE)) + 
          geom_point(aes(color=as.factor(Generation)), alpha=0.7) + 
          geom_smooth() +
          scale_color_manual(values=adjusted_blues) + 
          facet_wrap(~TRAIT, scales="free_y") +
          labs(x="Mass Per Plant BLUP", y="", color="Generation") +
          theme_bw()
ggsave("blups_vs_MASS_PER_PLANT.png", gee_wiz, width=14)


gee_wiz <- df_massper %>% 
        ggplot(aes(x=MASS_PER_PLANT_blup, y=VALUE, group=Generation)) + 
          geom_point(aes(color=as.factor(Generation)), alpha=0.7) + 
          geom_smooth(aes(color=as.factor(Generation))) +
          scale_color_manual(values=adjusted_blues) + 
          facet_wrap(~TRAIT, scales="free_y") +
          labs(x="Mass Per Plant  BLUP", y="", color="Generation") +
          theme_bw()
ggsave("blups_vs_MASS_PER_PLANT_by_generation.png", gee_wiz, width=14)



# seed-weight & total-mass plots
df_tm <- df %>% 
        pivot_longer(-c(TOTAL_MASS_blup, Generation), values_to = "VALUE", names_to="TRAIT")
df_tm$TRAIT <- tidy_text_substitution(df_tm$TRAIT)
golly <- df_tm %>% 
        ggplot(aes(x=TOTAL_MASS_blup, y=VALUE)) + 
          geom_point(aes(color=as.factor(Generation)), alpha=0.7) + 
          geom_smooth() +
          scale_color_manual(values=adjusted_blues) + 
          facet_wrap(~TRAIT, scales="free_y") +
          labs(x="Total Mass BLUP", y="", color="Generation") +
          theme_bw()
ggsave("blups_vs_TOTALMASS.png", golly, width=14)

golly <- df_tm %>% 
        ggplot(aes(x=TOTAL_MASS_blup, y=VALUE, group=Generation)) + 
          geom_point(aes(color=as.factor(Generation)), alpha=0.7) + 
          geom_smooth(aes(color=as.factor(Generation))) +
          scale_color_manual(values=adjusted_blues) + 
          facet_wrap(~TRAIT, scales="free_y") +
          labs(x="Total Mass BLUP", y="", color="Generation") +
          theme_bw()
ggsave("blups_vs_TOTALMASS_by_generation.png", golly, width=14)


# seed-weight & total-mass plots
df_sw <- df %>% 
        pivot_longer(-c(SEED_WEIGHT_100_blup, Generation), values_to = "VALUE", names_to="TRAIT")
df_sw$TRAIT <- tidy_text_substitution(df_sw$TRAIT)
golly <- df_sw %>% 
        ggplot(aes(x=SEED_WEIGHT_100_blup, y=VALUE)) + 
          geom_point(aes(color=as.factor(Generation)), alpha=0.7) + 
          geom_smooth() +
          scale_color_manual(values=adjusted_blues) + 
          facet_wrap(~TRAIT, scales="free_y") +
          labs(x="100-Seed Weight BLUP", y="", color="Generation") +
          theme_bw()
ggsave("blups_vs_100sw.png", golly, width=14)

golly <- df_sw %>% 
        ggplot(aes(x=SEED_WEIGHT_100_blup, y=VALUE, group=Generation)) + 
          geom_point(aes(color=as.factor(Generation)), alpha=0.7) + 
          geom_smooth(aes(color=as.factor(Generation)), alpha=0.5) +
          scale_color_manual(values=adjusted_blues) + 
          facet_wrap(~TRAIT, scales="free_y") +
          labs(x="100-Seed Weight BLUP", y="", color="Generation") +
          theme_bw()
ggsave("blups_vs_100sw_by_generation.png", golly, width=14)
