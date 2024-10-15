#!/usr/bin/env Rscript
#SBATCH --mem=30G
#SBATCH --time=02:00:00
#SBATCH --output=/rhome/jmarz001/bigdata/Ag-Competition/scripts/002_QC.stdout
#SBATCH -p short

library(tidyverse)
library(ggpubr)
library(data.table)
library(car)
library(gridExtra)
library(dunn.test)

setwd("/bigdata/koeniglab/jmarz001/Ag-Competition")
df <- fread("data/FITNESS.tsv")

#  Centered data 
#sw$FEC <- as.vector(scale(sw$FEC, center = TRUE, scale =TRUE))

# arrange data for facet plotting
df_long <- df %>%
  pivot_longer(cols=-c(Genotype, Condition, Generation), names_to='PHENOTYPE', values_to="VALUE")

#"#eca50b"


# check normality & plot trait distributions
ggplot(df_long, aes(VALUE)) +
  geom_density() +
  #geom_vline(aes(xintercept = mean(df2$Plants)), color = "#0c820c") +
  facet_wrap(~PHENOTYPE, scales="free") +
  theme_bw() +
  stat_summary(fun = mean, geom = "vline", orientation = "y", 
    aes(xintercept = after_stat(x), y = 0)) #+
  #labs(title = "Distribution of plot survival")
ggsave("results/trait_distributions.png")
# traits look like normal enough distribution
# survival has a bit of skew but it's not a trait I expect to be normal
# 100-seed-weight has overly high values





# filter for outlier values
summary(sw$TOTAL_MASS)

upper <- median(mix1$TOTAL_MASS, na.rm = T) + (2 * IQR(mix1$TOTAL_MASS, na.rm = T))
lower <- median(mix1$TOTAL_MASS, na.rm = T) - (2 * IQR(mix1$TOTAL_MASS, na.rm = T))

ggplot(X, aes(SEED_WEIGHT_100)) + geom_histogram()
  df[which(df$SEED_WEIGHT_100 > 30),]




## plot all trait distributions w various facets

ggplot(df_long, aes(VALUE, color=Condition)) +
  geom_density() +
  facet_wrap(~PHENOTYPE, scales="free") +
  theme_bw()+
  stat_summary(fun = mean, geom = "vline", orientation = "y", 
    aes(xintercept = after_stat(x), y = 0)) 
ggsave("results/trait_distributions_Wcondition.png", width=12)

ggplot(df_long, aes(VALUE, group=Generation, color=as.factor(Generation))) +
  geom_density(linewidth=0.75) +
  facet_wrap(~PHENOTYPE, scales="free") +
  theme_bw() +
  stat_summary(fun = mean, geom = "vline", orientation = "y", 
    aes(xintercept = after_stat(x), y = 0)) 
ggsave("results/trait_distributions_Wgeneration.png", width=14)

ggplot(df_long, aes(VALUE, color=as.factor(Generation))) +
  geom_density(aes(linetype=Condition), linewidth=0.75) +
  facet_wrap(~PHENOTYPE, scales="free") +
  theme_bw() +
  stat_summary(fun = mean, geom = "vline", orientation = "y", 
    aes(xintercept = after_stat(x), y = 0)) 
ggsave("results/trait_distributions_Wgeneration_Wcondition.png", width=14)




# set up for normality and variance equity tests
single <- df %>% filter(Condition == "single") 
collist <- c('FT', 'TOTAL_MASS', 'SEED_WEIGHT_100', 'SURVIVAL', 'SEED_COUNT', 'FECUNDITY', 'FITNESS', 'RELATIVE_FITNESS', 'AT_REL_FITNESS')
generationlist <- c(0, 18, 28, 50, 58)
par(mfrow = c(2, 3))

# test for normal distributions 
# w shapiro test
# and visualize w qq-plot
for(i in collist) {

  print(i)
  tmp <- df %>% select(c(Genotype, Generation, Condition, all_of(i)))

  for(g in generationlist) {

    print(g)
    tmp_gen <- tmp %>% filter(Generation == g)
    s <- shapiro.test(tmp_gen$FITNESS)
    print(s)

    qqnorm(p$TOTAL_MASS, main = paste0("Generation ", i))
    qqline(p$TOTAL_MASS)
  }
}

#TOTAL_MASS (SOME GROUPS NOT NORMAL)
# FT (SOME GROUPS NOT NORMAL)
# Fec (NOT NORMAL)
# Fitness (NOT ALL NORMAL)
# 100 seed weight (NORMAL)




# test for equality of variance between groups before anova
# w Levene test
for(i in collist){
  print(i)
  leveneTest(get(i) ~ as.factor(Generation), single) %>% print
}
## only flowering time has unequal variance between generations


### Not normal distributions &| not equal variance btwn trait-groups
### Kruskal Wallis and Dunn Tests
kruskal.test(TOTAL_MASS ~ Generation, single)
dunn.test(single$TOTAL_MASS, single$Generation)

kruskal.test(FECUNDITY ~ Generation, single)
dunn.test(single$FEC, single$Generation)

kruskal.test(ABS_FITNESS ~ Generation, single)
dunn.test(single$ABS_FITNESS, single$Generation)


# normally distributed traits & equal trait-group variance
### AVOVA/Tukey Post-hoc
ANOVA_Fec <- aov(Fecundity ~ as.factor(Generation), single)
summary(ANOVA_Fec)
TukeyHSD(ANOVA_Fec)

ANOVA_100 <- aov(SEED_WEIGHT_100 ~ as.factor(Generation), single)
summary(ANOVA_100)
TukeyHSD(ANOVA_100)
