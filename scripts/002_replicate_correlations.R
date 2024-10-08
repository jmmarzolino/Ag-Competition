#!/usr/bin/env Rscript
#SBATCH --ntasks=1
#SBATCH --mem=50G
#SBATCH --time=02:00:00
#SBATCH --job-name="Ag-Competition"
#SBATCH --output=/rhome/jmarz001/bigdata/Ag-Competition/002_replicate_correlations.stdout
#SBATCH -p short

library(tidyverse)
library(ggpubr)
library(data.table)

library(gridExtra)
library(car)
library(ggrepel)


setwd("/rhome/jmarz001/bigdata/Ag-Competition/data")
# load data
df <- fread("JOINED_PHENOTYPES.tsv")
df1 <- df %>% filter(Exp_year==2022)
df2 <- df %>% filter(Exp_year==2023)

# subset data
single1 <- df1 %>%
  filter(Condition=="single") %>%
  pivot_wider(id_cols=Genotype, names_from=Replicate, values_from=c(TOTAL_MASS, SEED_WEIGHT_100, FT))
  
single2 <- df2 %>%
  filter(Condition=="single") %>%
  pivot_wider(id_cols=Genotype, names_from=Replicate, values_from=c(TOTAL_MASS, SEED_WEIGHT_100, FT)) 

mix1 <- df1 %>%
  filter(Condition=="mixed") %>%
  pivot_wider(id_cols=Genotype, names_from=Replicate, values_from=c(TOTAL_MASS, SEED_WEIGHT_100, FT))
  
mix2 <- df2 %>%
  filter(Condition=="mixed") %>%
  pivot_wider(id_cols=Genotype, names_from=Replicate, values_from=c(TOTAL_MASS, SEED_WEIGHT_100, FT)) 




# calculate correlation across replicates
print("Trait Correlations Per Year")
for(i in c('single1', 'single2', 'mix1', 'mix2')){
  print(i)
  i <- get(i)
  print("Total Mass")
  print(cor(i$TOTAL_MASS_1, i$TOTAL_MASS_2, use="pairwise.complete.obs"))
  print("100 Seed Weight")
  print(cor(i$SEED_WEIGHT_100_1, i$SEED_WEIGHT_100_2, use="pairwise.complete.obs"))
  print("Flowering Time")
  print(cor(i$FT_1, i$FT_2, use="pairwise.complete.obs"))
}

# correlation options
#method=c("pearson", "spearman")
#use=c( ‘"everything"’,
#          ‘"all.obs"’, ‘"complete.obs"’, ‘"na.or.complete"’, or
#          ‘"pairwise.complete.obs"’.
#)

replicate_df <- df %>%
  select(-c(Plants, BED, ROW)) %>%
  #filter(Condition=="single") %>%
  pivot_longer(cols=c(FT, TOTAL_MASS, SEED_WEIGHT_100), names_to='PHENOTYPE', values_to="VALUE") %>%
  pivot_wider(names_from=Replicate, names_prefix="REP_", values_from='VALUE')




tm <- df %>%
  filter(Plants>0 & !is.na(TOTAL_MASS) & !is.na(SEED_WEIGHT_100)) %>% 
  select(-c(Plants, BED, ROW)) %>%
  #filter(Condition=="single") %>%
  #pivot_longer(cols=c(FT, TOTAL_MASS, SEED_WEIGHT_100), names_to='PHENOTYPE', values_to="VALUE") %>%
  select(c(Genotype, Condition, Replicate, Exp_year, TOTAL_MASS)) %>% 
  pivot_wider(names_from=Replicate, names_prefix="REP_", values_from='TOTAL_MASS')








# Creating a function to graph correlations between replicates
# x = dataframe
# y = phenotype, string
# z = year, string

graph_correlation <- function(x=df, y='TOTAL_MASS', z=""){
 x %>%
  filter(Plants>0 & !is.na(TOTAL_MASS) & !is.na(SEED_WEIGHT_100)) %>% 
  select(-c(Plants, BED, ROW)) %>%
  select(c('Genotype', 'Condition', 'Replicate', 'Exp_year', all_of(y))) %>% 
  pivot_wider(names_from=Replicate, names_prefix="REP_", values_from=y) %>%
  ggplot(aes(REP_1, REP_2), add = "reg.line") +
    geom_jitter() +
    geom_abline(slope = 1, intercept = 0, color = "red") +
    geom_smooth(method = "lm") +
    labs(x = "Replicate 1",
         y = "Replicate 2",
         title = paste("Correlation of", sep = " ", y, "Replicates", z)) +
    theme_bw() +
    facet_wrap(vars(Condition, Exp_year))
}


c1 <- graph_correlation(df) 
c2 <- graph_correlation(df, y='SEED_WEIGHT_100') 
c3 <- graph_correlation(df, y='FT') 

ggsave(c1, "../results/replicate_correlations_total_mass.png", width = 12, height = 12, dpi=300)
ggsave(c2, "../results/replicate_correlations_seed_weight_100.png", width = 12, height = 12, dpi=300)
ggsave(c3, "../results/replicate_correlations_FT.png", width = 12, height = 12, dpi=300)

png("../results/replicate_correlations_all.png")
ggarrange(c1, c2, c3)
dev.off()




## temporary code to identify extreme seed weight values
df %>% filter(SEED_WEIGHT_100 >30)








#p1 <- graph_correlation(cmp, "Total Weight (g)", "2021-2022") +
#  stat_cor(label.y = 170)
#p2 <- graph_correlation(cmp, "Relative Fecundity", "2021-2022") +
#  stat_cor(label.y = 2.5, label.x = 6.5) +
#  geom_text_repel(label = ifelse(cmp$REP_1 > 3 | cmp$REP_2 > 3,
 #                                cmp$Genotype,
  #                               ""), size = 3, hjust =1, max.overlaps = 30)



upper <- median(mix1$TOTAL_MASS, na.rm = T) + (2 * IQR(mix1$TOTAL_MASS, na.rm = T))
lower <- median(mix1$TOTAL_MASS, na.rm = T) - (2 * IQR(mix1$TOTAL_MASS, na.rm = T))
outlier_data <- mix1 %>% filter(TOTAL_MASS > upper | TOTAL_MASS < lower)





# Residual Plotting
cmp <- PHENO_MIXED_2022 %>% group_by(Genotype, Generation, Replicate) %>% summarise(TOTAL_MASS) %>% spread(key = Replicate, value = TOTAL_MASS) %>% ungroup()


plot_residuals <- function(cmp, lab_x="", lab_title="Residual Plot") {
    cmp <- na.omit(cmp)
    residual_df <- cmp %>% reframe(Residuals = resid(lm(REP_2 ~ REP_1)))
    cmp$Residuals <- residual_df$Residuals
    Standard_dev <- sd(cmp$Residuals)
    cmp <- cmp %>% mutate(Residuals = Residuals/Standard_dev)
    tmp_outlier <- cmp %>% filter(Residuals > 2 | Residuals < -2)

    g <- ggplot(cmp, aes(x = REP_1, y = Residuals)) +
      geom_jitter() +
      geom_hline(yintercept = 0, color = "red", line_type="dashed") +
      geom_hline(yintercept = 2, color = "blue") +
      geom_hline(yintercept = -2, color = "blue") +
      scale_y_continuous(breaks = seq(-4, 4, 1)) +
      geom_text_repel(label = ifelse(cmp$Residuals > 2 | cmp$Residuals < -2,
                                    cmp$Genotype,
                                    ""), size = 3, hjust = 1, max.overlaps = 30) +
      labs(x = lab_x,
          title = lab_title)

    return(g)
}

a <- plot_residuals(cmp, lab_x="Centered Fecundity", lab_title="Residual Plot - Mixed, Centered Fecundity, 2022")
a <- plot_residuals(cmp, lab_x="Total Weight (grams)", lab_title="Residual Plot - Mixed, Total Weight, 2022")






y <- arrangeGrob(a,b,c,d,e,f,g,h,i,j, nrow = 2, ncol = 5)
ggsave("../results/residuals.png", y, width = 26, height =14)
