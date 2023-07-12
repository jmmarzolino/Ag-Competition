#!/usr/bin/env Rscript

#SBATCH --ntasks=1
#SBATCH --mem=50G
#SBATCH --time=02:00:00
#SBATCH --job-name="Plot FT"
#SBATCH --output=/rhome/jmarz001/bigdata/Ag-Competition/scripts/001_compare_ft_over_years.stdout
#SBATCH -p short


setwd("/rhome/jmarz001/bigdata/Ag-Competition/data")
library(tidyverse)

ft22 <- read_delim("FT_2022.tsv", "\t")
ft22 <- ft22 %>% select(-c('2021BED', '2021ROW'))
ft22$number_of_plants <- as.numeric(ft22$number_of_plants)
ft22$Exp_year <- 2022

ft23 <- read_delim("FT_2023.tsv", "\t")
ft23 <- ft23 %>% select(-c('Bed_2022', 'Row_2022', 'PLOT_ID'))
ft23$Exp_year <- 2023

bright_eyes <- full_join(ft22, ft23, by=c('Genotypes', 'Condition', 'replicate', 'Flowering_Date'='FT_DAYS', 'Generation', 'number_of_plants'='Plot_Survival', 'Exp_year')) #%>% select(-'number_of_plants')

write_delim(bright_eyes, "FT_BOTH_YEARS.tsv", "\t")


### SIGNIFICANCE TESTS
# test relationship between flowering date (days between planting and spike emergence) and...
summary(aov(Flowering_Date ~ Exp_year, bright_eyes))
summary(aov(Flowering_Date ~ Condition + replicate + Generation + number_of_plants + Exp_year + Condition*Generation, bright_eyes))
summary(aov(Flowering_Date ~ Condition*Exp_year + replicate*Exp_year + Generation*Exp_year + number_of_plants*Exp_year + Exp_year + Condition*Generation, bright_eyes))


### CORRELATION
exp1 <- bright_eyes %>% filter(Exp_year==2022) %>% filter(!is.na(Flowering_Date)) %>% select(-c(number_of_plants, Generation, Exp_year)) %>% group_by(Genotypes, Condition) %>% summarise(geno_ft = mean(Flowering_Date)) %>% pivot_wider(names_from=c(Condition), values_from=geno_ft)
exp2 <- bright_eyes %>% filter(Exp_year==2023) %>% filter(!is.na(Flowering_Date)) %>% select(-c(number_of_plants, Generation, Exp_year)) %>% group_by(Genotypes, Condition) %>% summarise(geno_ft = mean(Flowering_Date)) %>% pivot_wider(names_from=c(Condition), values_from=geno_ft)

cor(exp1$mixed, exp1$single, use="complete.obs")
cor(exp2$mixed, exp2$single, use="complete.obs")
cor(exp1$single, exp2$single, use="complete.obs")
cor(exp1$mixed, exp2$mixed, use="complete.obs")




### Plotting
setwd("../results/")
## distribution (all)
png("distribution.png")
ggplot(bright_eyes, aes(x=Flowering_Date)) +
geom_density(alpha=0.5) +
theme_minimal() +
facet_wrap(~Exp_year, nrow=2)
dev.off()


## boxplot of experimental conditions
png("boxplot_condition.png")
ggplot(bright_eyes, aes(y=Flowering_Date, x=Condition, fill=Condition)) +
geom_boxplot() +
theme_minimal() +
facet_wrap(~Exp_year)
dev.off()

## distribution of experimental conditions
png("distribution_condition.png")
ggplot(bright_eyes, aes(x=Flowering_Date, group=Condition, color=Condition, fill=Condition)) +
geom_density(alpha=0.5) +
theme_minimal() +
facet_wrap(~Exp_year, nrow=2)
dev.off()



## boxplot of Generations
png("boxplot_generation.png")
bright_eyes %>% ggplot( aes(y=Flowering_Date, x=as.factor(Generation), group=as.factor(Generation), fill=as.factor(Generation))) +
geom_boxplot() +
theme_minimal() +
facet_wrap(~Exp_year)
dev.off()

## distribution of Generations
png("distribution_generation.png")
ggplot(bright_eyes, aes(x=Flowering_Date, group=as.factor(Generation), color=as.factor(Generation), fill=as.factor(Generation))) +
geom_density(alpha=0.5) +
theme_minimal() +
facet_wrap(~Exp_year, nrow=2)
dev.off()

## distribution of Generations - single (non-mixed) plots only
png("distribution_generation_single_condition.png")
bright_eyes %>% filter(Condition == "single") %>% ggplot(aes(x=Flowering_Date, group=as.factor(Generation), color=as.factor(Generation), fill=as.factor(Generation))) +
geom_density(alpha=0.5) +
theme_minimal() +
facet_wrap(~Exp_year, nrow=2)
dev.off()

## boxplot of Generations - single plots only
png("boxplot_generation_single_condition.png")
bright_eyes %>% filter(Condition == "single") %>% ggplot( aes(y=Flowering_Date, x=as.factor(Generation), group=as.factor(Generation), fill=as.factor(Generation))) +
geom_boxplot() +
theme_minimal() +
facet_wrap(~Exp_year)
dev.off()



## distribution of Generations - mixed plots only
png("distribution_generation_mixed_condition.png")
bright_eyes %>% filter(Condition == "mixed") %>% ggplot(aes(x=Flowering_Date, group=as.factor(Generation), color=as.factor(Generation), fill=as.factor(Generation))) +
geom_density(alpha=0.5) +
theme_minimal() +
facet_wrap(~Exp_year, nrow=2)
dev.off()

## boxplot of Generations - mixed plots only
png("boxplot_generation_mixed_condition.png")
bright_eyes %>% filter(Condition == "mixed") %>% ggplot( aes(y=Flowering_Date, x=as.factor(Generation), group=as.factor(Generation), fill=as.factor(Generation))) +
geom_boxplot() +
theme_minimal() +
facet_wrap(~Exp_year)
dev.off()



## boxplot of replicates
png("boxplot_replicates.png")
ggplot(bright_eyes, aes(y=Flowering_Date, x=replicate, group=replicate, fill=replicate)) +
geom_boxplot() +
theme_minimal() +
facet_wrap(~Exp_year)
dev.off()

## distribution of replicates
png("distribution_replicates.png")
ggplot(bright_eyes, aes(x=Flowering_Date, group=replicate, color=replicate, fill=replicate)) +
geom_density(alpha=0.5) +
theme_minimal() +
facet_wrap(~Exp_year, nrow=2)
dev.off()
