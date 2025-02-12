#!/usr/bin/env Rscript
#SBATCH --ntasks=1
#SBATCH --mem=30G
#SBATCH --time=02:00:00
#SBATCH --output=/rhome/jmarz001/bigdata/Ag-Competition/scripts/plotting/2_traits_over_generations.stdout
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

# filter trait w/out any variance
df <- df %>% select(-SEED_WEIGHT_100)

# add column for parent or progeny marker
df$pgroup <- "Parents"
df[which(df$Generation != 0), which(colnames(df) == "pgroup")] <- "Progeny"
#df <- df %>% filter(Generation != 50)
#df$Generation <- as.factor(df$Generation)

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
  geom_density(color="#084594", linewidth=1) +
  facet_wrap(~TRAIT, scales="free") +  
  labs(x="", y="density") +
  theme_bw(base_size=20) #+
  #stat_summary(fun = median, geom = "vline", orientation = "y", aes(xintercept = after_stat(x), y = 0), color="#eca50b", linewidth=1) 
ggsave("results/trait_distributions.png", g, width=12)


# one density line per generation
g <- ggplot(df_long, aes(VALUE, group=Generation, color=as.factor(Generation))) +
  geom_density(linewidth=1) +
  scale_color_manual(values=adjusted_blues, name="Generation") + 
  facet_wrap(~TRAIT, scales="free") +
  theme_bw(base_size=20) +
  labs(x="", y="density") +
  stat_summary(fun = median, geom = "vline", orientation = "y", aes(xintercept = after_stat(x), y = 0), linewidth=1) 
ggsave("results/trait_distributions_Wgeneration.png", g, width=16)


# density lines for parents and progeny groups
g <- ggplot(df_long, aes(VALUE, group=pgroup, color=as.factor(pgroup))) +
  geom_density(linewidth=1) +
  scale_color_manual(values=adjusted_blues[c(1, 5)], name="Generation") + 
  facet_wrap(~TRAIT, scales="free") +
  theme_bw(base_size=20) +
  labs(x="", y="density") +
  stat_summary(fun = median, geom = "vline", orientation = "y", aes(xintercept = after_stat(x), y = 0), linewidth=1) 
ggsave("results/trait_distributions_Wparentprogeny.png", g, width=16)














## boxplots comparing traits over generations
#a <- ggplot(df_long, aes(y=VALUE, x=Generation, group=Generation)) +
#      facet_wrap(~TRAIT, scales="free") +
#      geom_boxplot(width=0.5, outliers=F) +
#      geom_jitter(width=0.2) +
#      theme_bw(base_size=20) +
#      labs(x = "Generation", y = "") 
#ggsave("results/traits_over_generations_scatterplot.png", a, width = 14, height = 10)

# when using boxplot & jitter ....      
# some points are doubled b/c they're plotted as part of jitter & boxplot-outliers
# you can color boxplot outliers, or not include them (outliers=T/F, outlier.color/shape/size/alpha...)


# dot plot with generation-fit line
#t <- ggplot(df_long, aes(y=VALUE, x=Generation)) +
#      geom_jitter() +
#      geom_smooth(method="lm") + stat_regline_equation() +
#      theme_bw(base_size=20) +
#      facet_wrap(~TRAIT, scales="free")+
#      labs(x = "Generation", y = "") 
#ggsave("results/traits_over_generations_scatterplot_w_trendline.png", t)





# line for generation means...
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
      geom_jitter(data=df, aes(x=Generation, group=as.factor(Generation), y=FITNESS)) + 
      geom_line(data=gen_men, aes(x=Generation, y=FITNESS), color="blue", linewidth=1) + 
      theme_bw(base_size=18) +
      labs(x = "Generation", y = tidy_text_substitution("FITNESS"), title="Genotype Fitness over Generations") 

ggcombo <- ggarrange(a, b, c, d, e)
ggsave("results/traits_over_generations_scatterplots.png", ggcombo, width = 20, height = 10)