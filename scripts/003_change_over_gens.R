#!/usr/bin/env Rscript
#SBATCH --mem=50G
#SBATCH --time=02:00:00
#SBATCH --output=/rhome/jmarz001/bigdata/Ag-Competition/scripts/003_change_over_gens.stdout
#SBATCH -p koeniglab

setwd("/rhome/jmarz001/bigdata/Ag-Competition")
source("scripts/CUSTOM_FNS.R")

library(tidyverse)
library(data.table)
library(corrplot)
library(lme4)
#install_packages("lmeTest")
#library(lmeTest)


# read in data
df <- fread("data/FITNESS.tsv")
df <- df %>% filter(Condition == "single")



# check correlations trait between traits 
# remove highly correlated phenotypes from gwas
# fitness / atlas-fitness / 
traits_df <- pheno %>% select(c(ends_with("_scaled"))) 
#traits_df <- pheno %>% select(c(ends_with("_scaled"))) %>% select(-c(SEED_COUNT_scaled, RELATIVE_FITNESS_scaled, AT_REL_FITNESS_scaled)) 

x <- cor(traits_df, use="na.or.complete", method="spearman")
corrplot(x, method="color", type="upper", order="original", title="", mar=c(0,0,4,0), addCoef.col = "black")





## BASE STATISTICS
# summarise mean & variance
df %>% 
    group_by(Generation, Condition) %>% 
    summarise(across(where(is.numeric), list(mean=mean, var=var), .names="{.col}_{.fn}")) -> x
write_delim(x, "data/generations_trait_avg_var.tsv")
# write table out with generation/condition trait averages / summary statistics



# set up for normality and variance equity tests
single <- df %>% filter(Condition == "single") 
collist <- c('FT', 'TOTAL_MASS', 'SEED_WEIGHT_100', 'SURVIVAL', 'SEED_COUNT', 'FECUNDITY', 'FITNESS', 'RELATIVE_FITNESS', 'AT_REL_FITNESS')
generationlist <- c(0, 18, 28, 50, 58)


pdf("results/normality_test.pdf")

# test for normal distributions 
# w shapiro test
# and visualize w qq-plot
for(i in collist) {

  print(i)
  tmp <- single %>% select(c(Genotype, Generation, Condition, all_of(i))) %>% tibble

  for(g in generationlist) {

    print(g)
    tmp_gen <- tmp %>% filter(Generation == g)
    s <- shapiro.test(tmp_gen[,4][[1]])
    print(s)

    qqnorm(tmp_gen[,4][[1]], main = paste0("Generation ", i))
    qqline(tmp_gen[,4][[1]])
  }
}

dev.off()




# test for equality of variance between groups before anova
# w Levene test
for(i in collist){
  print(i)
  leveneTest(get(i) ~ as.factor(Generation), single) %>% print
}
## only flowering time has unequal variance between generations


### Not normal distributions &| not equal variance btwn trait-groups
# seed-weight-100, survival, fecundity, ft
### Kruskal Wallis and Dunn Tests
print('100 seed weight')
kruskal.test(SEED_WEIGHT_100 ~ Generation, single)
dunn.test(single$SEED_WEIGHT_100, single$Generation)

print('survival')
kruskal.test(SURVIVAL ~ Generation, single)
dunn.test(single$SURVIVAL, single$Generation)

print('flowering time')
kruskal.test(FT ~ Generation, single)
dunn.test(single$FT, single$Generation)


# normally distributed traits & equal trait-group variance
# total mass, seed-count, fitness, relative-fitness, at-rel-fitness
### AVOVA/Tukey Post-hoc
print('total mass')
lil <- aov(TOTAL_MASS ~ as.factor(Generation), single)
summary(lil)
TukeyHSD(lil)

print('seed count')
lil <- aov(SEED_COUNT ~ as.factor(Generation), single)
summary(lil)
TukeyHSD(lil)

print('fecundity')
lil <- aov(FECUNDITY ~ as.factor(Generation), single)
summary(lil)
TukeyHSD(lil)

print('fitness')
lil <- aov(FITNESS ~ as.factor(Generation), single)
summary(lil)
TukeyHSD(lil)

print('relative fitness')
lil <- aov(RELATIVE_FITNESS ~ as.factor(Generation), single)
summary(lil)
TukeyHSD(lil)

print('Atlas-relative fitness')
lil <- aov(AT_REL_FITNESS ~ as.factor(Generation), single)
summary(lil)
TukeyHSD(lil)

















### Plotting
## boxplots comparing conditions
trait_df <- df %>% pivot_longer(cols=c('FT', 'TOTAL_MASS', 'SEED_WEIGHT_100', 'SURVIVAL', 'SEED_COUNT', 'FECUNDITY', 'FITNESS', 'RELATIVE_FITNESS', 'AT_REL_FITNESS'), values_to="VALUE", names_to="trait")

ggplot(trait_df, aes(y=VALUE, x=Condition)) +
    geom_boxplot() +
    theme_minimal() +
    facet_wrap(~trait, scales="free")


## boxplots comparing traits over generations
ggplot(single, aes(y=FT, x=as.factor(Generation), group=as.factor(Generation), fill=as.factor(Generation))) +
    geom_boxplot() +
    theme_minimal() +
    facet_wrap(~trait, scales="free")

ggplot(single, aes(y=FT, x=as.factor(Generation), group=as.factor(Generation), fill=as.factor(Generation))) +
    geom_boxplot() +
    theme_minimal() +
    facet_wrap(~trait, scales="free")



