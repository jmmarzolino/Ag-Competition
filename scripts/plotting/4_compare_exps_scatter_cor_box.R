#!/usr/bin/env Rscript
#SBATCH --ntasks=1
#SBATCH --mem=30G
#SBATCH --time=02:00:00
#SBATCH --output=/rhome/jmarz001/bigdata/Ag-Competition/scripts/plotting/4_compare_exps_scatter_cor_box.stdout
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
g <- df %>% filter(days_to_heading < 120) %>% #remove excessive values that don't get plotted anyways
    ggplot(aes(y=FT, x=days_to_heading)) +
        geom_point() +
        geom_smooth(method="lm", color="black", se=T) + 
        stat_cor(method="spearman", cor.coef.name="rho") +
        theme_bw(base_size=16) +
        labs(x = "greenhouse", y = "field", title="Comparing Flowering Time across Environments") 

###   100 seed mass
s <- ggplot(df, aes(y=SEED_WEIGHT_100, x=X100_seed_mass)) +
    geom_point() +
    geom_smooth(method="lm", color="black", se=T) + 
    stat_cor(method="spearman", cor.coef.name="rho") +
    theme_bw(base_size=16) +
    labs(x = "greenhouse", y = "field", title="Comparing 100-Seed Mass across Environments") 

###   fecundity
f <- ggplot(df, aes(y=FECUNDITY, x=seed_estimate)) +
    geom_point() +
    geom_smooth(method="lm", color="black", se=T) + 
    stat_cor(method="spearman", cor.coef.name="rho") +
    theme_bw(base_size=16) +
    labs(x = "greenhouse", y = "field", title="Comparing Seed Production across Environments") 

### print all together
zzz <- ggarrange(g, s, f, ncol=3)
ggsave(filename="results/cross_exp_scatter_cor.png", plot=zzz, height=7, width=(7*3)+1, units="in")





## 4 specific cross-experiment comparison plots
### field flowering time vs. 
### 100 seed weight, fecundity
### greenhouse: 100 seed weight, fecundity
###   100 seed mass
f_sw <- df %>% 
    ggplot(aes(y=SEED_WEIGHT_100, x=FT)) +
        geom_point() +
        geom_smooth(method="lm", color="black", se=T) + 
        stat_cor(method="spearman", cor.coef.name="rho") +
        theme_bw(base_size=16) +
        labs(x = "field flowering time", y = "field 100-seed weight") 
###   fecundity
f_fe <- df %>% 
    ggplot(aes(y=FECUNDITY, x=FT)) +
        geom_point() +
        geom_smooth(method="lm", color="black", se=T) + 
        stat_cor(method="spearman", cor.coef.name="rho") +
        theme_bw(base_size=16) +
        labs(x = "field flowering time", y = "field fecundity") 

###   GH 100 seed mass
f_ghsw <- df %>% 
    ggplot(aes(y=X100_seed_mass, x=FT)) +
        geom_point() +
        geom_smooth(method="lm", color="black", se=T) + 
        stat_cor(method="spearman", cor.coef.name="rho") +
        theme_bw(base_size=16) +
        labs(x = "field flowering time", y = "greenhouse 100-seed weight") 

###   GH fecundity
f_ghfe <- df %>% 
    ggplot(aes(y=seed_estimate, x=FT)) +
        geom_point() +
        geom_smooth(method="lm", color="black", se=T) + 
        stat_cor(method="spearman", cor.coef.name="rho") +
        theme_bw(base_size=16) +
        labs(x = "field flowering time", y = "greenhouse fecundity") 

### print all together
zzz <- ggarrange(f_sw, f_fe, f_ghsw, f_ghfe, ncol=2, nrow=2)
ggsave(filename="results/cross_exp_scatter_cor_fig2_4panels.png", plot=zzz, height=7*2, width=7*2, units="in")



## All samples colored by generation
g <- df %>% filter(days_to_heading < 120) %>% #remove excessive values that don't get plotted anyways
    ggplot(aes(y=FT, x=days_to_heading, color=Generation)) +
        geom_point() +
        geom_smooth(method="lm", color="black", se=T) + 
        stat_cor(method="spearman", cor.coef.name="rho") +
        theme_bw(base_size=16) +
        labs(x = "greenhouse", y = "field", title="Comparing Flowering Time across Environments") 

###   fecundity
f <- ggplot(df, aes(y=FECUNDITY, x=seed_estimate, color=Generation)) +
    geom_point() +
    geom_smooth(method="lm", color="black", se=T) + 
    stat_cor(method="spearman", cor.coef.name="rho") +
    theme_bw(base_size=16) +
    labs(x = "greenhouse", y = "field", title="Comparing Seed Production across Environments") 

###   100 seed mass
s <- ggplot(df, aes(y=SEED_WEIGHT_100, x=X100_seed_mass, color=Generation)) +
    geom_point() +
    geom_smooth(method="lm", color="black", se=T) + 
    stat_cor(method="spearman", cor.coef.name="rho") +
    theme_bw(base_size=16) +
    labs(x = "greenhouse", y = "field", title="Comparing Yield across Environments") 

### print all together
zzz <- ggarrange(g, s, f, ncol=3)
ggsave(filename="results/cross_exp_scatter_cor_color_gens.png", plot=zzz, height=7, width=(7*3)+1, units="in")




## All parents
parent <- df %>% filter(Generation == 0)
g <- parent %>% filter(days_to_heading < 120) %>% #remove excessive values that don't get plotted anyways
    ggplot(aes(y=FT, x=days_to_heading)) +
        geom_point() +
        geom_smooth(method="lm", color="black", se=T) + 
        stat_cor(method="spearman", cor.coef.name="rho") +
        theme_bw(base_size=16) +
        labs(x = "greenhouse", y = "field", title="Comparing Flowering Time across Environments") 

###   fecundity
f <- ggplot(parent, aes(y=FECUNDITY, x=seed_estimate)) +
    geom_point() +
    geom_smooth(method="lm", color="black", se=T) + 
    stat_cor(method="spearman", cor.coef.name="rho") +
    theme_bw(base_size=16) +
    labs(x = "greenhouse", y = "field", title="Comparing Seed Production across Environments") 

###   100 seed mass
s <- ggplot(parent, aes(y=SEED_WEIGHT_100, x=X100_seed_mass)) +
    geom_point() +
    geom_smooth(method="lm", color="black", se=T) + 
    stat_cor(method="spearman", cor.coef.name="rho") +
    theme_bw(base_size=16) +
    labs(x = "greenhouse", y = "field", title="Comparing Yield across Environments") 

### print all together
zzz <- ggarrange(g, s, f, ncol=3)
ggsave(filename="results/cross_exp_scatter_cor_parents.png", plot=zzz, height=7, width=(7*3)+1, units="in")






## All progeny
prog <- df %>% filter(Generation != 0)

g <- ggplot(prog, aes(y=FT, x=days_to_heading)) +
    geom_point() +
    geom_smooth(method="lm", color="black", se=T) + 
    stat_cor(method="spearman", cor.coef.name="rho") +
    theme_bw(base_size=16) +
    labs(x = "greenhouse", y = "field", title="Comparing Flowering Time across Environments", subtitle="Progeny") 


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

###   fecundity
f <- ggplot(prog, aes(y=FECUNDITY, x=seed_estimate)) +
    geom_point() +
    geom_smooth(method="lm", color="black", se=T) + 
    stat_cor(method="spearman", cor.coef.name="rho") +
    theme_bw(base_size=16) +
    labs(x = "greenhouse", y = "field", title="Comparing Seed Production across Environments", subtitle="Progeny") 

###   100 seed mass
s <- ggplot(prog, aes(y=SEED_WEIGHT_100, x=X100_seed_mass)) +
    geom_point() +
    geom_smooth(method="lm", color="black", se=T) + 
    stat_cor(method="spearman", cor.coef.name="rho") +
    theme_bw(base_size=16) +
    labs(x = "greenhouse", y = "field", title="Comparing Yield across Environments", subtitle="Progeny") 

### print all together
zzz <- ggarrange(g, s, f, ncol=3)
ggsave(filename="results/cross_exp_scatter_cor_progeny.png", plot=zzz, height=7, width=(7*3)+1, units="in")





####    Boxplots
#II. paired box plots (of scaled phenotypes?) for field and greenhouse experiments, trait values over generations
adjusted_blues <- brewer.pal(7, "Blues")[3:7]

###   Flowering Time
# scale data for comparison across experiments
tmp2 <- df %>% select(c(FT, days_to_heading, Generation)) %>% 
    mutate(across(-c(Generation), ~(scale(.) %>% as.vector))) %>% 
    pivot_longer(-Generation, names_to="experiment")

tmp2$experiment <- gsub("FT", "field", tmp2$experiment)
tmp2$experiment <- gsub("days_to_heading", "greenhouse", tmp2$experiment)


#c <- ggplot(tmp2, aes(y=value, color=experiment, x=Generation)) +
#    geom_boxplot() +
#    theme_bw(base_size=16) +
#    labs(title="Comparing Experiments", subtitle="Flowering Time", y="Scaled Days to heading") 

# switch the order of greenhouse and field datasets in plot
d <- tmp2 %>% 
    mutate(experiment = factor(experiment, levels=c("greenhouse", "field"))) %>%
        ggplot(aes(y=value, x=experiment, color=Generation)) +
            geom_boxplot() +
            scale_color_manual(values=adjusted_blues) +
            theme_bw(base_size=16) +
            labs(title="Comparing Experiments", subtitle="Flowering Time", y="Scaled Days to heading") 

#ggcombo <- ggarrange(c, d)
ggsave("results/cross_exp_boxplot_FT.png", d, width = 7, height = 7)

###   100 Seed weight
# scale data for comparison across experiments
tmp2 <- df %>% select(c(SEED_WEIGHT_100, X100_seed_mass, Generation)) %>% 
    mutate(across(-c(Generation), ~(scale(.) %>% as.vector))) %>% 
    pivot_longer(-Generation, names_to="experiment")

tmp2$experiment <- gsub("SEED_WEIGHT_100", "field", tmp2$experiment)
tmp2$experiment <- gsub("X100_seed_mass", "greenhouse", tmp2$experiment)

#c <- ggplot(tmp2, aes(y=value, color=experiment, x=Generation)) +
 #   geom_boxplot() +
  #  theme_bw(base_size=16) +
   # labs(title="Comparing Experiments", subtitle="100-Seed Mass", y="Scaled mass for 100 seeds (g)") 
d <- tmp2 %>% 
    mutate(experiment = factor(experiment, levels=c("greenhouse", "field"))) %>%
        ggplot(aes(y=value, x=experiment, color=Generation)) +
            geom_boxplot() +
            scale_color_manual(values=adjusted_blues) +
            theme_bw(base_size=16) +
            labs(title="Comparing Experiments", subtitle="100-Seed Mass", y="Scaled mass for 100 seeds (g)") 

#ggcombo <- ggarrange(c, d)
ggsave("results/cross_exp_boxplot_100sm.png", d, width = 7, height = 7)



###   Fecundity / number of seeds
# scale data for comparison across experiments
tmp2 <- df %>% select(c(FECUNDITY, seed_estimate, Generation)) %>% 
    mutate(across(-c(Generation), ~(scale(.) %>% as.vector))) %>% 
    pivot_longer(-Generation, names_to="experiment")

tmp2$experiment <- gsub("FECUNDITY", "field", tmp2$experiment)
tmp2$experiment <- gsub("seed_estimate", "greenhouse", tmp2$experiment)

#c <- ggplot(tmp2, aes(y=value, color=experiment, x=Generation)) +
 #   geom_boxplot() +
  #  theme_bw(base_size=16) +
   # labs(title="Comparing Experiments", subtitle="Estimated seed number", y="Scaled Seeds") 
d <- tmp2 %>% 
    mutate(experiment = factor(experiment, levels=c("greenhouse", "field"))) %>%
        ggplot(aes(y=value, x=experiment, color=Generation)) +
            geom_boxplot() +
            scale_color_manual(values=adjusted_blues) +
            theme_bw(base_size=16) +
            labs(title="Comparing Experiments", subtitle="Estimated seed number", y="Scaled Seeds") 

#ggcombo <- ggarrange(c, d)
ggsave("results/cross_exp_boxplot_fec.png", d, width = 7, height = 7)
