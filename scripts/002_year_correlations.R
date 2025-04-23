#!/usr/bin/env Rscript
#SBATCH --mem=50G
#SBATCH --time=02:00:00
#SBATCH --output=/rhome/jmarz001/bigdata/Ag-Competition/scripts/002_year_correlations.stdout
#SBATCH -p short

library(tidyverse)
library(ggpubr)
library(data.table)
library(ggrepel)
library(Cairo)

setwd("/rhome/jmarz001/bigdata/Ag-Competition/data")
df <- fread("JOINED_PHENOTYPES.tsv")

# average phenotypes across replicates 
df_mean <- df %>% 
    mutate("MASS_PER_PLANT"=TOTAL_MASS/Germination) %>% 
    select(-c(Germination, BED, ROW, Replicate)) %>%
    group_by(Genotype, Condition, Exp_year) %>% 
    summarise(across(where(is.numeric), \(x) mean(x, na.rm=T))) %>%
    ungroup()

df_rep <- df_mean %>% 
    pivot_longer(cols=c("FT", "TOTAL_MASS", "MASS_PER_PLANT", "SEED_WEIGHT_100"), names_to="PHENOTYPES", values_to="VALUE") %>% 
    pivot_wider(names_from=Exp_year, names_prefix="YEAR_", values_from="VALUE")
 
df_rep %>%
    group_by(Condition, PHENOTYPES) %>%
    summarise("COR" = cor(YEAR_2022, YEAR_2023, use="pairwise.complete.obs"))



# plot correlation between years
g <- df_rep %>%
    ggplot(aes(x=YEAR_2022, y=YEAR_2023), add = "reg.line") +
    geom_jitter() +
    geom_abline(slope = 1, intercept = 0, color = "red") +
    geom_smooth(method = "lm") +
    #stat_cor() +
    labs(x = "Year 2022",
        y = "Year 2023",
        title = paste("Correlation of", sep = " ", y)) +
    theme_bw() +
    facet_wrap(Condition~PHENOTYPES, scales="free", nrow=2,ncol=4)


ggsave("../results/annual_correlations.png", g, width=24, height=24, dpi=300)
