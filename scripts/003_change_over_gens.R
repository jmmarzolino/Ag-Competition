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


### BASE STATISTICS
# parent genotype averages
df4 %>% 
    filter(Generation == 0) %>% 
    filter(Condition == "single") 

# summarise mean & variance
df4 %>% 
    group_by(Generation, Condition) %>% 
    summarise(across(where(is.numeric), mean))
df4 %>% 
    group_by(Generation, Condition) %>% 
    summarise(across(where(is.numeric), var))








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