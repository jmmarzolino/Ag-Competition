#!/usr/bin/env Rscript
#SBATCH --ntasks=1
#SBATCH --mem=30G
#SBATCH --time=02:00:00
#SBATCH --output=/rhome/jmarz001/bigdata/Ag-Competition/scripts/plotting/3_traits_over_gens_scatter_and_line.stdout
#SBATCH -p koeniglab

# This script plots scatterplots with linear regressions for each trait in the single subpopulation, as well as distributions for these traits

library(tidyverse)
library(ggpubr)
library(ggplot2)
library(ggsci)
library(RColorBrewer)

setwd("/bigdata/koeniglab/jmarz001/Ag-Competition")
source("scripts/CUSTOM_FNS.R")
df <- read_delim("data/combined_CCII.tsv")
df <- add_generation(df)
df$Generation <- as.factor(df$Generation)
#comparing trait values over generations
## scatter- & box- plots 
#I. scatter plots comparing growth environments, for all traits
#    1. all samples
#    2. all progeny
#    3. each generation
#II. paired box plots (of scaled phenotypes?) for field and greenhouse experiments, trait values over generations



###    Scatter Plots

## All samples
# parent lines aren't included because they were identified by name in CC II greenhouse data and I don't want to manually match them...


## All progeny
g <- ggplot(df, aes(y=FT, x=days_to_heading)) +
    geom_point() +
    geom_smooth(method="lm", color="black", se=F) + 
    stat_cor(method="spearman", cor.coef.name="rho") +
    theme_bw(base_size=16) +
    labs(x = "Days to Heading (greenhouse)", y = paste0(tidy_text_substitution("FT"), " (field)"), title="Comparing Flowering Time across Environments") 


## ggpub syntax
#ggscatter(df, x = "days_to_heading", y = "FT", add = "reg.line") + stat_regline_equation(label.y=120)

###
#cor(df$days_to_heading, df$FT, use="pairwise.complete.obs", method="spearman")
# r2 is literally just r*r, r^2


#cor.coef.name = "rho" # rho is for spearman coef
# label.x.npc = "left",
# label.y.npc = "top",
#c('right', 'left', 'center', 'centre', 'middle') 
#c('bottom', 'top', 'center', 'centre', 'middle')

###   total mass
m1 <- ggplot(df, aes(y=TOTAL_MASS, x=seed_mass)) +
    geom_point() +
    geom_smooth(method="lm", color="black", se=F) + 
    stat_cor(method="spearman", cor.coef.name="rho") +
    theme_bw(base_size=16) +
    labs(x = "Seed mass (greenhouse)", y = paste0(tidy_text_substitution("TOTAL_MASS"), " (field)"), title="Comparing Yield across Environments") 


m2 <- ggplot(df, aes(y=MASS_PER_PLANT, x=seed_mass)) +
    geom_point() +
    geom_smooth(method="lm", color="black", se=F) + 
    stat_cor(method="spearman", cor.coef.name="rho") +
    theme_bw(base_size=16) +
    labs(x = "Seed mass (greenhouse)", y = paste0(tidy_text_substitution("MASS_PER_PLANT"), " (field)"), title="Comparing Yield across Environments") 


###   fecundity
f <- ggplot(df, aes(y=FECUNDITY, x=seed_estimate)) +
    geom_point() +
    geom_smooth(method="lm", color="black", se=F) + 
    stat_cor(method="spearman", cor.coef.name="rho") +
    theme_bw(base_size=16) +
    labs(x = "Seed estimate (greenhouse)", y = paste0(tidy_text_substitution("FECUNDITY"), " (field)"), title="Comparing Seed Production across Environments") 

###   100 seed mass
s <- ggplot(df, aes(y=SEED_WEIGHT_100, x=X100_seed_mass)) +
    geom_point() +
    geom_smooth(method="lm", color="black", se=F) + 
    stat_cor(method="spearman", cor.coef.name="rho") +
    theme_bw(base_size=16) +
    labs(x = "100-Seed weight (greenhouse)", y = paste0(tidy_text_substitution("SEED_WEIGHT_100"), " (field)"), title="Comparing Yield across Environments") 

### print all together
zzz <- ggarrange(g, m1, m2, f, s)
ggsave(filename="results/cross_exp_scatter_cor.png", plot=zzz, height=18, width=24, units="in")






## Each generation
for(i in c(18, 28, 50, 58)) {

    tmp <- df %>% filter(Generation == i) 

    ## flowering time
    g <- ggplot(tmp, aes(y=FT, x=days_to_heading)) +
        geom_point() +
        geom_smooth(method="lm", color="black", se=F) + 
        stat_cor(method="spearman", cor.coef.name="rho") +
        theme_bw(base_size=16) +
        labs(x = "Days to Heading (greenhouse)", y = paste0(tidy_text_substitution("FT"), " (field)"), title="Comparing Flowering Time across Environments", subtitle=paste0("Individuals from Generation ", i)) 

    ###   total mass
    m1 <- ggplot(tmp, aes(y=TOTAL_MASS, x=seed_mass)) +
        geom_point() +
        geom_smooth(method="lm", color="black", se=F) + 
        stat_cor(method="spearman", cor.coef.name="rho") +
        theme_bw(base_size=16) +
        labs(x = "Seed mass (greenhouse)", y = paste0(tidy_text_substitution("TOTAL_MASS"), " (field)"), title="Comparing Yield across Environments", subtitle=paste0("Individuals from Generation ", i)) 

    m2 <- ggplot(tmp, aes(y=MASS_PER_PLANT, x=seed_mass)) +
        geom_point() +
        geom_smooth(method="lm", color="black", se=F) + 
        stat_cor(method="spearman", cor.coef.name="rho") +
        theme_bw(base_size=16) +
        labs(x = "Seed mass (greenhouse)", y = paste0(tidy_text_substitution("MASS_PER_PLANT"), " (field)"), title="Comparing Yield across Environments", subtitle=paste0("Individuals from Generation ", i))  

    ###   fecundity
    f <- ggplot(tmp, aes(y=FECUNDITY, x=seed_estimate)) +
        geom_point() +
        geom_smooth(method="lm", color="black", se=F) + 
        stat_cor(method="spearman", cor.coef.name="rho") +
        theme_bw(base_size=16) +
        labs(x = "Seed estimate (greenhouse)", y = paste0(tidy_text_substitution("FECUNDITY"), " (field)"), title="Comparing Seed Production across Environments", subtitle=paste0("Individuals from Generation ", i))  


    ###   100 seed mass
    s <- ggplot(tmp, aes(y=SEED_WEIGHT_100, x=X100_seed_mass)) +
        geom_point() +
        geom_smooth(method="lm", color="black", se=F) + 
        stat_cor(method="spearman", cor.coef.name="rho") +
        theme_bw(base_size=16) +
        labs(x = "100-Seed weight (greenhouse)", y = paste0(tidy_text_substitution("SEED_WEIGHT_100"), " (field)"), title="Comparing Yield across Environments", subtitle=paste0("Individuals from Generation ", i)) 

    zzz_gen <- ggarrange(g, m1, m2, f, s)
    ggsave(filename=paste0("results/cross_exp_scatter_cor_gen", i, ".png"), plot=zzz_gen, height=18, width=24, units="in")
}






####    Boxplots
#II. paired box plots (of scaled phenotypes?) for field and greenhouse experiments, trait values over generations
adjusted_blues <- brewer.pal(7, "Blues")[3:7]

###   Flowering Time
tmp <- df %>% 
        select(c(FT, days_to_heading, Generation)) %>% 
        pivot_longer(-Generation, names_to="experiment")

tmp$experiment <- gsub("FT", "field", tmp$experiment)
tmp$experiment <- gsub("days_to_heading", "greenhouse", tmp$experiment)

a <- ggplot(tmp, aes(y=value, color=experiment, x=Generation)) +
    geom_boxplot() +
    theme_bw(base_size=16) +
    labs(title="Comparing Experiments", subtitle="Flowering Time", y="Days to heading")

b <- ggplot(tmp, aes(y=value, x=experiment, color=Generation)) +
    geom_boxplot() +
    scale_color_manual(values=adjusted_blues) +
    theme_bw(base_size=16) +
    labs(title="Comparing Experiments", subtitle="Flowering Time", y="Days to heading")

# scale data for comparison across experiments
tmp2 <- df %>% select(c(FT, days_to_heading, Generation)) %>% 
    mutate(across(-c(Generation), ~(scale(.) %>% as.vector))) %>% 
    pivot_longer(-Generation, names_to="experiment")

tmp2$experiment <- gsub("FT", "field", tmp2$experiment)
tmp2$experiment <- gsub("days_to_heading", "greenhouse", tmp2$experiment)

c <- ggplot(tmp2, aes(y=value, color=experiment, x=Generation)) +
    geom_boxplot() +
    theme_bw(base_size=16) +
    labs(title="Comparing Experiments", subtitle="Flowering Time", y="Scaled Days to heading") 
d <- ggplot(tmp2, aes(y=value, x=experiment, color=Generation)) +
    geom_boxplot() +
    scale_color_manual(values=adjusted_blues) +
    theme_bw(base_size=16) +
    labs(title="Comparing Experiments", subtitle="Flowering Time", y="Scaled Days to heading") 

ggcombo <- ggarrange(a, b, c, d)
ggsave("results/cross_exp_boxplot_FT.png", ggcombo, width = 20, height = 10)






###   Seed Weight

tmp <- df %>% select(c(TOTAL_MASS, seed_mass, Generation)) %>% pivot_longer(-Generation, names_to="experiment")

tmp <- df %>% 
        select(c(TOTAL_MASS, seed_mass, Generation)) %>% 
        pivot_longer(-Generation, names_to="experiment")

tmp$experiment <- gsub("TOTAL_MASS", "field", tmp$experiment)
tmp$experiment <- gsub("seed_mass", "greenhouse", tmp$experiment)

a <- ggplot(tmp, aes(y=value, color=experiment, x=Generation)) +
    geom_boxplot() +
    theme_bw(base_size=16) +
    labs(title="Comparing Experiments", subtitle="Total Seed Mass", y="mass (g)")

b <- ggplot(tmp, aes(y=value, x=experiment, color=Generation)) +
    geom_boxplot() +
    scale_color_manual(values=adjusted_blues) +
    theme_bw(base_size=16) +
    labs(title="Comparing Experiments", subtitle="Total Seed Mass", y="mass (g)")

# scale data for comparison across experiments
tmp2 <- df %>% select(c(TOTAL_MASS, seed_mass, Generation)) %>% 
    mutate(across(-c(Generation), ~(scale(.) %>% as.vector))) %>% 
    pivot_longer(-Generation, names_to="experiment")

tmp2$experiment <- gsub("TOTAL_MASS", "field", tmp2$experiment)
tmp2$experiment <- gsub("seed_mass", "greenhouse", tmp2$experiment)

c <- ggplot(tmp2, aes(y=value, color=experiment, x=Generation)) +
    geom_boxplot() +
    theme_bw(base_size=16) +
    labs(title="Comparing Experiments", subtitle="Total Seed Mass", y="Scaled mass (g)") 
d <- ggplot(tmp2, aes(y=value, x=experiment, color=Generation)) +
    geom_boxplot() +
    scale_color_manual(values=adjusted_blues) +
    theme_bw(base_size=16) +
    labs(title="Comparing Experiments", subtitle="Total Seed Mass", y="Scaled mass (g)") 

ggcombo <- ggarrange(a, b, c, d)
ggsave("results/cross_exp_boxplot_MASS.png", ggcombo, width = 20, height = 10)


##  Mass v. 2
tmp <- df %>% select(c(MASS_PER_PLANT, seed_mass, Generation)) %>% pivot_longer(-Generation, names_to="experiment")

tmp <- df %>% 
        select(c(MASS_PER_PLANT, seed_mass, Generation)) %>% 
        pivot_longer(-Generation, names_to="experiment")

tmp$experiment <- gsub("MASS_PER_PLANT", "field", tmp$experiment)
tmp$experiment <- gsub("seed_mass", "greenhouse", tmp$experiment)

a <- ggplot(tmp, aes(y=value, color=experiment, x=Generation)) +
    geom_boxplot() +
    theme_bw(base_size=16) +
    labs(title="Comparing Experiments", subtitle="Per-plant Seed Mass", y="mass (g)")

b <- ggplot(tmp, aes(y=value, x=experiment, color=Generation)) +
    geom_boxplot() +
    scale_color_manual(values=adjusted_blues) +
    theme_bw(base_size=16) +
    labs(title="Comparing Experiments", subtitle="Per-plant Seed Mass", y="mass (g)")

# scale data for comparison across experiments
tmp2 <- df %>% select(c(MASS_PER_PLANT, seed_mass, Generation)) %>% 
    mutate(across(-c(Generation), ~(scale(.) %>% as.vector))) %>% 
    pivot_longer(-Generation, names_to="experiment")

tmp2$experiment <- gsub("MASS_PER_PLANT", "field", tmp2$experiment)
tmp2$experiment <- gsub("seed_mass", "greenhouse", tmp2$experiment)

c <- ggplot(tmp2, aes(y=value, color=experiment, x=Generation)) +
    geom_boxplot() +
    theme_bw(base_size=16) +
    labs(title="Comparing Experiments", subtitle="Per-plant Seed Mass", y="Scaled mass (g)") 
d <- ggplot(tmp2, aes(y=value, x=experiment, color=Generation)) +
    geom_boxplot() +
    scale_color_manual(values=adjusted_blues) +
    theme_bw(base_size=16) +
    labs(title="Comparing Experiments", subtitle="Per-plant Seed Mass", y="Scaled mass (g)") 

ggcombo <- ggarrange(a, b, c, d)
ggsave("results/cross_exp_boxplot_MASS2.png", ggcombo, width = 20, height = 10)



###   100 Seed weight

tmp <- df %>% select(c(SEED_WEIGHT_100, X100_seed_mass, Generation)) %>% pivot_longer(-Generation, names_to="experiment")

tmp$experiment <- gsub("SEED_WEIGHT_100", "field", tmp$experiment)
tmp$experiment <- gsub("X100_seed_mass", "greenhouse", tmp$experiment)

a <- ggplot(tmp, aes(y=value, color=experiment, x=Generation)) +
    geom_boxplot() +
    theme_bw(base_size=16) +
    labs(title="Comparing Experiments", subtitle="100-Seed Mass", y="mass for 100 seeds (g)")

b <- ggplot(tmp, aes(y=value, x=experiment, color=Generation)) +
    geom_boxplot() +
    scale_color_manual(values=adjusted_blues) +
    theme_bw(base_size=16) +
    labs(title="Comparing Experiments", subtitle="100-Seed Mass", y="mass for 100 seeds (g)")


# scale data for comparison across experiments
tmp2 <- df %>% select(c(SEED_WEIGHT_100, X100_seed_mass, Generation)) %>% 
    mutate(across(-c(Generation), ~(scale(.) %>% as.vector))) %>% 
    pivot_longer(-Generation, names_to="experiment")

tmp2$experiment <- gsub("SEED_WEIGHT_100", "field", tmp2$experiment)
tmp2$experiment <- gsub("X100_seed_mass", "greenhouse", tmp2$experiment)

c <- ggplot(tmp2, aes(y=value, color=experiment, x=Generation)) +
    geom_boxplot() +
    theme_bw(base_size=16) +
    labs(title="Comparing Experiments", subtitle="100-Seed Mass", y="Scaled mass for 100 seeds (g)") 
d <- ggplot(tmp2, aes(y=value, x=experiment, color=Generation)) +
    geom_boxplot() +
    scale_color_manual(values=adjusted_blues) +
    theme_bw(base_size=16) +
    labs(title="Comparing Experiments", subtitle="100-Seed Mass", y="Scaled mass for 100 seeds (g)") 

ggcombo <- ggarrange(a, b, c, d)
ggsave("results/cross_exp_boxplot_100sm.png", ggcombo, width = 20, height = 10)





###   Fecundity / number of seeds
tmp <- df %>% select(c(FECUNDITY, seed_estimate, Generation)) %>% pivot_longer(-Generation, names_to="experiment")

tmp$experiment <- gsub("FECUNDITY", "field", tmp$experiment)
tmp$experiment <- gsub("seed_estimate", "greenhouse", tmp$experiment)

a <- ggplot(tmp, aes(y=value, color=experiment, x=Generation)) +
    geom_boxplot() +
    theme_bw(base_size=16) +
    labs(title="Comparing Experiments", subtitle="Estimated seed number", y="Seeds")

b <- ggplot(tmp, aes(y=value, x=experiment, color=Generation)) +
    geom_boxplot() +
    scale_color_manual(values=adjusted_blues) +
    theme_bw(base_size=16) +
    labs(title="Comparing Experiments", subtitle="Estimated seed number", y="Seeds")


# scale data for comparison across experiments
tmp2 <- df %>% select(c(FECUNDITY, seed_estimate, Generation)) %>% 
    mutate(across(-c(Generation), ~(scale(.) %>% as.vector))) %>% 
    pivot_longer(-Generation, names_to="experiment")

tmp2$experiment <- gsub("FECUNDITY", "field", tmp2$experiment)
tmp2$experiment <- gsub("seed_estimate", "greenhouse", tmp2$experiment)

c <- ggplot(tmp2, aes(y=value, color=experiment, x=Generation)) +
    geom_boxplot() +
    theme_bw(base_size=16) +
    labs(title="Comparing Experiments", subtitle="Estimated seed number", y="Scaled Seeds") 
d <- ggplot(tmp2, aes(y=value, x=experiment, color=Generation)) +
    geom_boxplot() +
    scale_color_manual(values=adjusted_blues) +
    theme_bw(base_size=16) +
    labs(title="Comparing Experiments", subtitle="Estimated seed number", y="Scaled Seeds") 

ggcombo <- ggarrange(a, b, c, d)
ggsave("results/cross_exp_boxplot_fec.png", ggcombo, width = 20, height = 10)
