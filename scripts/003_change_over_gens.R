#!/usr/bin/env Rscript
#SBATCH --mem=50G
#SBATCH --time=02:00:00
#SBATCH --output=/rhome/jmarz001/bigdata/Ag-Competition/scripts/003_change_over_gens.stdout
#SBATCH -p short

setwd("/rhome/jmarz001/bigdata/Ag-Competition/data")
source("../scripts/CUSTOM_FNS.R")

library(tidyverse)
#library(lme4)
#install_packages("lmeTest")
#library(lmeTest)

# read in data
df <- read_delim("JOINED_PHENOTYPES.tsv")
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
df4 <- df3 %>% filter(Genotype != "48_5")
df4 <- full_join(df4, AT)


df4$AT_REL_FITNESS <- df4$FITNESS / AT$FITNESS


# explore traits' mean and variance by different factors
for(i in 5:ncol(df)){
    y <- colnames(df[,i])

    df %>% 
        group_by(Generation) %>% 
        summarise(across(where(is.numeric), \(x) mean(x, na.rm=T), \(x) var(x, na.rm=T), n=n())) %>%
        print()

    # parent genotypes avg flowering date
    df %>% 
        filter(Generation == 0) %>% 
        filter(Condition == "single") %>% 
        group_by(Genotype) %>% 
        summarise(across(where(is.numeric), \(x) mean(x, na.rm=T), \(x) var(x, na.rm=T), n=n())) %>%
        print()

}

### Atlas is parent 48


#    pivot_wider(names_from=Condition, values_from=c(average, variance))


write_delim(report, "avg_FT_parents.tsv", "\t")














### SIGNIFICANCE TESTS
#anova(lmer()) # anova for fixed effects
#summary(aov())
#ranova(lmer()) # for random effects
#t.test(extra ~ group, data = sleep)

#x <-aov(lmer(FT ~ Condition + (1|Generation:Genotype) + (1|Replicate), df))
#x <-aov(lmer(FT ~ Condition + Generation + (1 + Generation|Genotype) + (1|Replicate)), df) # try to plot this line!
# Fixed: Condition, Generation
# Random: Genotype, Replicate

summary(aov(FT ~ Replicate, df))
summary(aov(FT ~ Condition, df))
summary(aov(FT ~ Replicate + Condition, df))
summary(aov(FT ~ Generation, df))

# combinations of factors
summary(aov(FT ~ Condition + Replicate + Generation + Plot_Survival, df))
summary(aov(FT ~ Condition*Generation, df))
summary(aov(FT ~ Condition + Replicate + Generation + Plot_Survival, df))
summary(aov(FT ~ Replicate + Generation + Condition + Plot_Survival + Generation*Plot_Survival + Replicate*Plot_Survival + Generation*Plot_Survival + Condition*Generation, df))
















### Plotting
## boxplot of experimental conditions
df %>% ggplot( aes(y=FT, x=as.factor(Generation), group=as.factor(Generation), fill=as.factor(Generation))) +
df %>% filter(Condition == "single") %>% ggplot( aes(y=FT, x=as.factor(Generation), group=as.factor(Generation), fill=as.factor(Generation))) +

save_name <- paste0("boxplot", factor, ".png")
png(save_name)



ggplot(df, aes(y=FT, x=Replicate, group=Replicate, fill=Replicate)) +


dev.off()