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
library(car)
library(dunn.test)
#install_packages("lmeTest")
#library(lmeTest)


# read in data
df <- fread("data/FITNESS.tsv")
df <- df %>% filter(Condition == "single") %>% select(-Condition)




## BASE STATISTICS
# summarise mean & variance
x <- df %>% 
    group_by(Generation) %>% 
    summarise(across(where(is.numeric), list(mean=mean, var=var), .names="{.col}_{.fn}")) 
write_delim(x, "data/generations_trait_avg_var.tsv")
# write table out with generation/condition trait averages / summary statistics



# set up for normality and variance equity tests
collist <- c('FT', 'TOTAL_MASS', 'SEED_WEIGHT_100', 'SURVIVAL', 'SEED_COUNT', 'FECUNDITY', 'FITNESS', 'RELATIVE_FITNESS', 'AT_REL_FITNESS')
generationlist <- c(0, 18, 28, 50, 58)


pdf("results/normality_test.pdf")

# test for normal distributions 
# w shapiro test
# and visualize w qq-plot
for(i in collist) {

  print(i)
  tmp <- df %>% select(c(Genotype, Generation, all_of(i))) %>% tibble

  for(g in generationlist) {

    print(g)
    tmp_gen <- tmp %>% filter(Generation == g)
    s <- shapiro.test(tmp_gen[,3][[1]])
    print(s)

    qqnorm(tmp_gen[,3][[1]], main = paste0("Generation ", i))
    qqline(tmp_gen[,3][[1]])
  }
}

dev.off()




# test for equality of variance between groups before anova
# w Levene test
for(i in collist){
  print(i)
  leveneTest(get(i) ~ as.factor(Generation), df) %>% print
}
## only flowering time has unequal variance between generations


### Not normal distributions &| not equal variance btwn trait-groups
# seed-weight-100, survival, fecundity, ft
### Kruskal Wallis and Dunn Tests
print('100 seed weight')
kruskal.test(SEED_WEIGHT_100 ~ Generation, df)
dunn.test(df$SEED_WEIGHT_100, df$Generation)

print('survival')
kruskal.test(SURVIVAL ~ Generation, df)
dunn.test(df$SURVIVAL, df$Generation)

print('flowering time')
kruskal.test(FT ~ Generation, df)
dunn.test(df$FT, df$Generation)


# normally distributed traits & equal trait-group variance
# total mass, seed-count, fitness, relative-fitness, at-rel-fitness
### AVOVA/Tukey Post-hoc
print('total mass')
lil <- aov(TOTAL_MASS ~ as.factor(Generation), df)
summary(lil)
TukeyHSD(lil)

print('seed count')
lil <- aov(SEED_COUNT ~ as.factor(Generation), df)
summary(lil)
TukeyHSD(lil)

print('fecundity')
lil <- aov(FECUNDITY ~ as.factor(Generation), df)
summary(lil)
TukeyHSD(lil)

print('fitness')
lil <- aov(FITNESS ~ as.factor(Generation), df)
summary(lil)
TukeyHSD(lil)

print('relative fitness')
lil <- aov(RELATIVE_FITNESS ~ as.factor(Generation), df)
summary(lil)
TukeyHSD(lil)

print('Atlas-relative fitness')
lil <- aov(AT_REL_FITNESS ~ as.factor(Generation), df)
summary(lil)
TukeyHSD(lil)

