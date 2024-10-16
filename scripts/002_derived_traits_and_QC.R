#!/usr/bin/env Rscript
#SBATCH --mem=30G
#SBATCH --time=02:00:00
#SBATCH --output=/rhome/jmarz001/bigdata/Ag-Competition/scripts/002_QC.stdout
#SBATCH -p short

library(tidyverse)
library(ggpubr)
library(data.table)
#library(car)
#library(gridExtra)
#library(dunn.test)

setwd("/bigdata/koeniglab/jmarz001/Ag-Competition")

df <- read_delim("data/JOINED_PHENOTYPES.tsv")
# remove any rows without genotype
df <- df %>% 
    filter(!is.na(Genotype)) %>%
    select(-c(BED, ROW, Replicate)) %>%
    group_by(Genotype, Condition, Exp_year) %>%
    summarise(across(where(is.numeric), \(x) mean(x, na.rm=T))) 

# fn adds generation col based on 'Genotype' col
df <- add_generation(df)


# calculate derived phenotypes

## survival = (plants/10) 
over_sown <- df %>% filter(Plants > 10)
# for plots have 11 or 12 seeds planted, 
# assume max number planted is the same
over_sown$SURVIVAL <- 1
# plants / 11 or 12 will be 100% germination

df2 <- df %>% filter(Plants <= 10)
df2 <- df2 %>% mutate(SURVIVAL = Plants / 10)

df3 <- full_join(df2, over_sown, by=c('Genotype', 'Condition', 'Exp_year', 'Plants', 'FT', 'TOTAL_MASS', 'SEED_WEIGHT_100', 'Generation', 'SURVIVAL'))
df3 <- df3 %>% group_by(Genotype, Condition) %>% summarise(across(where(is.numeric), mean)) %>% select(-Exp_year)


# fecundity
## fecundity = seed produced per plot
## total seed weight / seed-weight-100/100
df3 <- df3 %>% 
    mutate(PER_SEED_WEIGHT = SEED_WEIGHT_100/100) %>%
    mutate(SEED_COUNT = TOTAL_MASS / PER_SEED_WEIGHT) %>%
    mutate(FECUNDITY = SEED_COUNT / Plants) %>% 
    select(-c(PER_SEED_WEIGHT, Plants))

## fitness = survival * fecundity
df3 <- df3 %>% mutate(FITNESS = SURVIVAL * FECUNDITY) 
df3$RELATIVE_FITNESS <- df3$FITNESS / max(df3$FITNESS)


# fitness relative to Atlas (parent #48)
AT <- df3 %>% filter(Genotype == "48_5") %>% summarise(across(where(is.numeric), mean))
AT$Condition <- 'single'

df4 <- df3 %>% filter(Genotype != "48_5")
df4 <- full_join(df4, AT)

df4$AT_REL_FITNESS <- df4$FITNESS / AT$FITNESS
# write out data frame w derived phenotpyes
write_delim(df4, "data/FITNESS.tsv")



##########
df <- fread("data/FITNESS.tsv")

#  Centered data 
#df$FEC <- as.vector(scale(df$FEC, center = TRUE, scale =TRUE))

# arrange data for facet plotting
df_long <- df %>%
  pivot_longer(cols=-c(Genotype, Condition, Generation), names_to='PHENOTYPE', values_to="VALUE")

#"#eca50b"


# check normality & plot trait distributions
g <- ggplot(df_long, aes(VALUE)) +
  geom_density() +
  #geom_vline(aes(xintercept = mean(df2$Plants)), color = "#0c820c") +
  facet_wrap(~PHENOTYPE, scales="free") +
  theme_bw() +
  stat_summary(fun = mean, geom = "vline", orientation = "y", 
    aes(xintercept = after_stat(x), y = 0)) #+
  #labs(title = "Distribution of plot survival")
ggsave("results/trait_distributions.png", g)
# traits look like normal enough distribution
# survival has a bit of skew but it's not a trait I expect to be normal
# 100-seed-weight has overly high values





# filter for outlier values
print('SEED_WEIGHT_100')
summary(df$SEED_WEIGHT_100)

upper <- median(df$SEED_WEIGHT_100) + (3 * IQR(df$SEED_WEIGHT_100))
lower <- median(df$SEED_WEIGHT_100) - (3 * IQR(df$SEED_WEIGHT_100))

g <- ggplot(df, aes(SEED_WEIGHT_100)) + geom_histogram() + 
  geom_vline(aes(xintercept=upper, color = "#F31919")) + 
  geom_vline(aes(xintercept=lower), color = "#F31919") + 
  geom_vline(aes(xintercept=median(df$SEED_WEIGHT_100), color = "#F31919"), linetype="dashed") + 
  theme_bw()

ggsave("results/seed_weight_outlier_distribution.png", g)
# filter outlier 100 seed weight
df[which(df$SEED_WEIGHT_100 > 8),]
df <- df %>% filter(SEED_WEIGHT_100 < 8)



## plot all trait distributions w various facets

g <- ggplot(df_long, aes(VALUE, color=Condition)) +
  geom_density() +
  facet_wrap(~PHENOTYPE, scales="free") +
  theme_bw()+
  stat_summary(fun = mean, geom = "vline", orientation = "y", 
    aes(xintercept = after_stat(x), y = 0)) 
ggsave("results/trait_distributions_Wcondition.png", g, width=12)

g <- ggplot(df_long, aes(VALUE, group=Generation, color=as.factor(Generation))) +
  geom_density(linewidth=0.75) +
  facet_wrap(~PHENOTYPE, scales="free") +
  theme_bw() +
  stat_summary(fun = mean, geom = "vline", orientation = "y", 
    aes(xintercept = after_stat(x), y = 0)) 
ggsave("results/trait_distributions_Wgeneration.png", g, width=14)

g <- ggplot(df_long, aes(VALUE, color=as.factor(Generation))) +
  geom_density(aes(linetype=Condition), linewidth=0.75) +
  facet_wrap(~PHENOTYPE, scales="free") +
  theme_bw() +
  stat_summary(fun = mean, geom = "vline", orientation = "y", 
    aes(xintercept = after_stat(x), y = 0)) 
ggsave("results/trait_distributions_Wgeneration_Wcondition.png", g, width=14)




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
