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
    select(-c(Plants, BED, ROW, Replicate)) %>%
    group_by(Genotype, Condition, Exp_year) %>% 
    summarise(across(where(is.numeric), \(x) mean(x, na.rm=T))) %>%
    ungroup()

# scale data
df_mean$FT <- scale(df_mean$FT)[,1]
df_mean$TOTAL_MASS <- scale(df_mean$TOTAL_MASS)[,1]
df_mean$SEED_WEIGHT_100 <- scale(df_mean$SEED_WEIGHT_100)[,1]

# functions to print phenotypes yearly correlation
# x = data frame
# y = phenotype, string

annual_correlation <- function(x, y){

    # correlation between years
    x2 <- x %>% 
        select(c(Genotype, Condition, Exp_year, all_of(y))) %>% 
        pivot_wider(names_from = Exp_year, names_prefix = "YEAR_", values_from = y) %>%
        select(-Genotype) %>%
        na.omit() 
    
    x2 %>%
        group_by(Condition) %>%
        summarise("COR" = cor(YEAR_2022, YEAR_2023, use="pairwise.complete.obs")) %>%
        print

    # plot correlation between years
    g1 <- x2 %>%
            ggplot(aes(YEAR_2022, YEAR_2023), add = "reg.line") +
            geom_jitter() +
            geom_abline(slope = 1, intercept = 0, color = "red") +
            geom_smooth(method = "lm") +
            #stat_cor() +
            labs(x = "Year 2022",
                y = "Year 2023",
                title = paste("Correlation of", sep = " ", y)) +
            theme_bw() +
            facet_wrap(~Condition)

    return(g1)
}         


# loop over each phenotype
# calculate trait correlation btwn years & graph
y=c('TOTAL_MASS', 'FT', 'SEED_WEIGHT_100')

p1 <- annual_correlation(df_mean, y=y[1])
p2 <- annual_correlation(df_mean, y=y[2])
p3 <- annual_correlation(df_mean, y=y[3])

ppp <- ggarrange(p1, p2, p3)
ggsave("../results/annual_correlations.png", ppp) 
#, width=20, height=20, dpi=300)



# labelling genotype source of outlier points
# geom_text_repel(label = ifelse(sc_22$RES > 2 | sc_22$RES < -2,
#                                scaled$Genotype, ""), 
#                                size = 3, hjust = 1, max.overlaps = 30) +

