#!/usr/bin/env Rscript
#SBATCH --mem=50G
#SBATCH --time=02:00:00
#SBATCH --output=/rhome/jmarz001/bigdata/Ag-Competition/scripts/002_replicate_correlations.stdout
#SBATCH -p short

library(tidyverse)
library(ggpubr)
library(data.table)
library(ggrepel)


setwd("/rhome/jmarz001/bigdata/Ag-Competition/data")
# load data
df <- fread("JOINED_PHENOTYPES.tsv")
df <- df %>%
  select(-c(BED, ROW)) %>%
  mutate("MASS_PER_PLANT"=TOTAL_MASS/Germination) 

# correlation options
#method=c("pearson", "spearman", "kendall")
#use=c( ‘"everything"’, ‘"all.obs"’, ‘"complete.obs"’, ‘"na.or.complete"’, ‘"pairwise.complete.obs"’)

replicate_df <- df %>%
  pivot_longer(cols=c(FT, TOTAL_MASS, SEED_WEIGHT_100, MASS_PER_PLANT), names_to='PHENOTYPE', values_to="VALUE") %>%
  pivot_wider(names_from=Replicate, names_prefix="REP_", values_from='VALUE')

# another way of summarizing trait replicate correlations by year, condition, and phenotype
replicate_df %>% group_by(Exp_year, Condition, PHENOTYPE) %>% summarise('correlation'=cor(REP_1, REP_2, use="pairwise.complete.obs", method="pearson"))


#replicate_df %>% group_by(Exp_year, Condition, PHENOTYPE) %>% summarise('correlation'=cor(REP_1, REP_2, use="pairwise.complete.obs", method="spearman"))
#replicate_df %>% group_by(Exp_year, Condition, PHENOTYPE) %>% summarise('correlation'=cor(REP_1, REP_2, use="pairwise.complete.obs", method="kendall"))






# Creating a function to graph correlations between replicates
# x = dataframe
# y = phenotype, string

graph_correlation <- function(x=df, y='TOTAL_MASS'){
 x %>%
  filter(Germination>0 & !is.na(TOTAL_MASS) & !is.na(SEED_WEIGHT_100)) %>% 
  select(-c(Germination)) %>%
  select(c('Genotype', 'Condition', 'Replicate', 'Exp_year', all_of(y))) %>% 
  pivot_wider(names_from=Replicate, names_prefix="REP_", values_from=y) %>%
  ggplot(aes(REP_1, REP_2), add = "reg.line") +
    geom_jitter() +
    geom_abline(slope = 1, intercept = 0, color = "red") +
    geom_smooth(method = "lm") +
    labs(x = "Replicate 1",
         y = "Replicate 2",
         title = paste("Correlation of", sep = " ", y, "Replicates")) +
    theme_bw() +
    facet_wrap(vars(Condition, Exp_year))
}


c1 <- graph_correlation(df) 
c1a <- graph_correlation(df, y='MASS_PER_PLANT') 
c2 <- graph_correlation(df, y='SEED_WEIGHT_100') 
c3 <- graph_correlation(df, y='FT') 

ggsave("../results/replicate_correlations_total_mass.png", c1, width = 12, height = 12, dpi=300)
ggsave("../results/replicate_correlations_mass_per_plant.png", c1a, width = 12, height = 12, dpi=300)

ggsave("../results/replicate_correlations_seed_weight_100.png", c2, width = 12, height = 12, dpi=300)
ggsave("../results/replicate_correlations_FT.png", c3, width = 12, height = 12, dpi=300)

gc <- ggarrange(c1, c1a, c2, c3)
ggsave("../results/replicate_correlations_all.png", gc, width=24, height=24)






# Residual Plotting

## scale data
scaled <- df %>% select(-c(Germination)) 
#scaled$FT <- scale(scaled$FT)[,1]
#scaled$TOTAL_MASS <- scale(scaled$TOTAL_MASS)[,1]
#scaled$SEED_WEIGHT_100 <- scale(scaled$SEED_WEIGHT_100)[,1]


plot_residuals <- function(scaled, y='TOTAL_MASS', condition="single", year=2022) {
    scaled <- na.omit(scaled)
    sc2 <- scaled %>% 
        filter(Condition == condition) %>%
        select(c(Genotype, Replicate, Exp_year, all_of(y))) %>%
        pivot_wider(names_from=Replicate, names_prefix="REP_", values_from=y)

    sc_22 <- sc2 %>% filter(Exp_year == year) %>% select(-Exp_year) %>% na.omit()

    d <- sc_22 %>% reframe("RES" = resid(lm(REP_2 ~ REP_1))) 
    sc_22$RES <- d$RES / sd(d$RES)
    print(sc_22 %>% filter(RES > 2 | RES < -2))


    g1 <- ggplot(sc_22, aes(x = REP_1, y = RES)) +
      geom_jitter() +
      geom_hline(yintercept = 0, color = "red") +
      geom_hline(yintercept = 2, color = "blue") +
      geom_hline(yintercept = -2, color = "blue") +
      scale_y_continuous(breaks = seq(-4, 4, 1)) +
      ylim(c(-4,4)) +
      geom_text_repel(label = ifelse(sc_22$RES > 2 | sc_22$RES < -2,
                                scaled$Genotype, ""), 
                                size = 3, hjust = 1, max.overlaps = 30) +
      labs(x = "Replicate 1", y = "Residuals", title = paste0("Residual Plot - ", year, ", ", condition), subtitle = y)

    # g2 is just a qq plot
    #g2 <- ggplot(sc_22, aes(x = REP_2, y = RES)) +
     # geom_jitter() +
     # geom_hline(yintercept = 0, color = "red") +
     # geom_hline(yintercept = 2, color = "blue") +
     # geom_hline(yintercept = -2, color = "blue") +
     # scale_y_continuous(breaks = seq(-4, 4, 1)) +
     # ylim(c(-4,4)) +
     # geom_text_repel(label = ifelse(sc_22$RES > 2 | sc_22$RES < -2,
     #                           scaled$Genotype, ""), 
     #                           size = 3, hjust = 1, max.overlaps = 30) +
     # labs(x = "Replicate 2", y = "Residuals", title = paste0("Residual Plot - ", year, ", ", condition), subtitle = y)

    #g <- ggarrange(g1, g2)
    return(g1)
}


y=c('TOTAL_MASS', 'MASS_PER_PLANT', 'FT', 'SEED_WEIGHT_100')

for(i in 1:3){
    p1 <- plot_residuals(scaled, y=y[i], condition="single", year=2022)
    p2 <- plot_residuals(scaled, y=y[i], condition="mixed", year=2022)

    p3 <- plot_residuals(scaled, y=y[i], condition="single", year=2023)
    p4 <- plot_residuals(scaled, y=y[i], condition="mixed", year=2023)

    ppp <- ggarrange(p1, p2, p3, p4, ncol=2, nrow=2)
    ggsave(paste0("../results/residuals_", y[i], ".png"), ppp, width=20, height=20, dpi=300)
}
