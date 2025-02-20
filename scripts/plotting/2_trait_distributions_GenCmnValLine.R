#!/usr/bin/env Rscript
#SBATCH --ntasks=1
#SBATCH --mem=30G
#SBATCH --time=02:00:00
#SBATCH --output=/rhome/jmarz001/bigdata/Ag-Competition/scripts/plotting/2_trait_distributions_GenCmnValLine.stdout
#SBATCH -p koeniglab

# plot trait distributions, with each generation 5 most common values for generations 15 & 58 pointed out

library(tidyverse)
library(ggpubr)
library(ggplot2)
library(RColorBrewer)

setwd("/bigdata/koeniglab/jmarz001/Ag-Competition")
source("scripts/CUSTOM_FNS.R")
df <- read_delim("data/trait_BLUPs.tsv")
df <- add_generation(df)

# filter trait w/out any variance
df <- df %>% select(-SEED_WEIGHT_100)

# the lowest color value is too light, so adjust the color scale down one
#display.brewer.pal(6, "Blues")
adjusted_blues <- brewer.pal(7, "Blues")[3:7]


### PLOTTING

# round trait blups as a way of binning to common values
df_rnd <- df
df_rnd$FT <- round(df$FT, 2)
df_rnd$Plants <- round(df$Plants, 3)
df_rnd$FECUNDITY <- round(df$FECUNDITY)
df_rnd$MASS_PER_PLANT <- round(df$MASS_PER_PLANT)


### find the most common phenotypes for early and late generations
# ft
F18 <- as.numeric(names(sort(table(df_rnd[which(df_rnd$Generation==18), which(colnames(df_rnd)=="FT")]), decreasing=T)[1:5]))
F58 <- as.numeric(names(sort(table(df_rnd[which(df_rnd$Generation==58), which(colnames(df_rnd)=="FT")]), decreasing=T)[1:5]))
FT_cmn <- tibble(F18, F58) %>% pivot_longer(names_to="Generation", cols=everything())

# germination
F18 <- as.numeric(names(sort(table(df_rnd[which(df_rnd$Generation==18), which(colnames(df_rnd)=="Plants")]), decreasing=T)[1:5]))
F58 <- as.numeric(names(sort(table(df_rnd[which(df_rnd$Generation==58), which(colnames(df_rnd)=="Plants")]), decreasing=T)[1:5]))
GM_cmn <- tibble(F18, F58) %>% pivot_longer(names_to="Generation", cols=everything())

# fecundity
F18 <- as.numeric(names(sort(table(df_rnd[which(df_rnd$Generation==18), which(colnames(df_rnd)=="FECUNDITY")]), decreasing=T)[1:5]))
F58 <- as.numeric(names(sort(table(df_rnd[which(df_rnd$Generation==58), which(colnames(df_rnd)=="FECUNDITY")]), decreasing=T)[1:5]))
FEC_cmn <- tibble(F18, F58) %>% pivot_longer(names_to="Generation", cols=everything())

# mass per plant
F18 <- as.numeric(names(sort(table(df_rnd[which(df_rnd$Generation==18), which(colnames(df_rnd)=="MASS_PER_PLANT")]), decreasing=T)[1:5]))
F58 <- as.numeric(names(sort(table(df_rnd[which(df_rnd$Generation==58), which(colnames(df_rnd)=="MASS_PER_PLANT")]), decreasing=T)[1:5]))
MASSPER_cmn <- tibble(F18, F58) %>% pivot_longer(names_to="Generation", cols=everything())



## plot trait distributions
a <- 
ggplot() +
  geom_density(data=df, aes(FT), linewidth=1) +
  geom_vline(data=FT_cmn, aes(xintercept=value, color=Generation), linewidth=1) +
  scale_color_manual(values=c("#6BAED6", "#084594")) + 
  labs(x="", y="density", title="Flowering Time Distribution", subtitle="with most common values colored") +
  theme_bw(base_size=18) +
  guides(color = "none")


b <- 
ggplot() +
  geom_density(data=df, aes(Plants), linewidth=1) +
  geom_vline(data=GM_cmn, aes(xintercept=value, color=Generation), linewidth=1) +
  scale_color_manual(values=c("#6BAED6", "#084594")) + 
  labs(x="", y="density", title="Germination Distribution", subtitle="with most common values colored") +
  theme_bw(base_size=18)


c <-
ggplot() +
  geom_density(data=df, aes(FECUNDITY), linewidth=1) +
  geom_vline(data=FEC_cmn, aes(xintercept=value, color=Generation), linewidth=1) +
  scale_color_manual(values=c("#6BAED6", "#084594")) + 
  labs(x="", y="density", title="Fecundity Distribution", subtitle="with most common values colored") +
  theme_bw(base_size=18) +
  guides(color = "none")


d <- 
ggplot() +
  geom_density(data=df, aes(MASS_PER_PLANT), linewidth=1) +
  geom_vline(data=MASSPER_cmn, aes(xintercept=value, color=Generation), linewidth=1) +
  scale_color_manual(values=c("#6BAED6", "#084594")) + 
  labs(x="", y="density", title="Mass Per Plant Distribution", subtitle="with most common values colored") +
  theme_bw(base_size=18) +
  guides(color = "none")


plots <- ggarrange(a, b, c, d)
ggsave("results/trait_distributions_GenCmnValLine.png", plots)
