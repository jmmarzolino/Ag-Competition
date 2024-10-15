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
# make list of trait columns
collist <- c('FT', 'TOTAL_MASS', 'SEED_WEIGHT_100', 'SURVIVIAL', 'SEED_COUNT', 'FECUNDITY', 'FITNESS', 'RELATIVE_FITNESS', 'AT_REL_FITNESS')

# loop over each trait column
for(i in collist){
    y <- df4 %>% select(c(Genotype, Generation, Condition, all_of(i)))

    summary(aov(unlist(y[,4]) ~ Condition, y))
    summary(aov(unlist(y[,4]) ~ Generation, y))
    summary(aov(unlist(y[,4]) ~ Condition*Generation, y))




    # combinations of factors
    summary(aov(unlist(y2[,4]) ~ Condition + Generation + SURVIVAL, y2))
    summary(aov(unlist(y[,4]) ~ Condition + Generation + SURVIVAL, y))

}


# check whether combinations of factors are significant / relevant
collist2 <- c('FT', 'TOTAL_MASS', 'SEED_WEIGHT_100', 'SEED_COUNT', 'FECUNDITY', 'FITNESS', 'RELATIVE_FITNESS', 'AT_REL_FITNESS')
y2 <- df4 %>% select(c(Genotype, Generation, Condition, SURVIVAL, all_of(collist2)))

for(i in 5:ncol(y2)) {
    print(colnames(y2[,i]))
    # test factor combinations w anova & AIC
    x <- aov(unlist(y2[,i]) ~ Generation + Condition + SURVIVAL + Generation*SURVIVAL, y2)
    summary(x) %>% print
    AIC(x) %>% print

    y <- aov(unlist(y2[,i]) ~ Generation + Condition + SURVIVAL, y2)
    summary(y) %>% print
    AIC(y) %>% print

    aov(unlist(y2[,i]) ~ Generation + Condition, y2) %>% AIC %>% print
    aov(unlist(y2[,i]) ~ Generation + SURVIVAL, y2) %>% AIC %>% print
    aov(unlist(y2[,i]) ~ SURVIVAL + Condition, y2) %>% AIC %>% print
    aov(unlist(y2[,i]) ~ Generation, y2) %>% AIC %>% print
    aov(unlist(y2[,i]) ~ Condition, y2) %>% AIC %>% print
    aov(unlist(y2[,i]) ~ SURVIVAL, y2) %>% AIC %>% print
    }

SEED_WEIGHT_100 ~ Generation



lm(AT_REL_FITNESS ~ Generation + Condition + SURVIVAL + Generation*SURVIVAL, y2)
# all 4 factors/combos get best AIC, only indv factors significant



#anova(lmer()) # anova for fixed effects
#ranova(lmer()) # for random effects
#t.test(y ~ gS, data = x)

# Fixed: Condition, Generation
# Random: Genotype, Replicate
lmer(unlist(y[,4]) ~ Condition + (1|Generation:Genotype), df4) %>% AIC
lmer(unlist(y[,4]) ~ Condition + Generation + (1 |Genotype), df4) %>% AIC










### Plotting
## boxplots comparing conditions
trait_df <- df4 %>% pivot_longer(cols=c('FT', 'TOTAL_MASS', 'SEED_WEIGHT_100', 'SURVIVAL', 'SEED_COUNT', 'FECUNDITY', 'FITNESS', 'RELATIVE_FITNESS', 'AT_REL_FITNESS'), values_to="VALUE", names_to="trait")

ggplot(trait_df, aes(y=VALUE, x=Condition)) +
geom_boxplot() +
theme_minimal() +
facet_wrap(~trait, scales="free")






ggplot( aes(y=FT, x=as.factor(Generation), group=as.factor(Generation), fill=as.factor(Generation))) +
df %>% filter(Condition == "single") %>% ggplot( aes(y=FT, x=as.factor(Generation), group=as.factor(Generation), fill=as.factor(Generation))) +

save_name <- paste0("boxplot", factor, ".png")
png(save_name)



ggplot(df, aes(y=FT, x=Replicate, group=Replicate, fill=Replicate)) +


dev.off()