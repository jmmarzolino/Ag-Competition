#!/usr/bin/env Rscript

#SBATCH --ntasks=1
#SBATCH --mem=50G
#SBATCH --time=02:00:00
#SBATCH --job-name="Plot FT"
#SBATCH --output=/rhome/jmarz001/bigdata/Ag-Competition/scripts/001_compare_ft_over_years.stdout
#SBATCH -p short

# set up workspace
setwd("/rhome/jmarz001/bigdata/Ag-Competition/data")
library(tidyverse)

# load and format data files
ft22 <- read_delim("FT_2021_2022.tsv", "\t")
ft22 <- ft22 %>% select(-c('BED_2021', 'ROW_2021'))
ft22$Plants <- as.numeric(ft22$Plants)
ft22$Exp_year <- 2022

ft23 <- read_delim("FT_2023.tsv", "\t")
ft23 <- ft23 %>% select(-c('BED_2022', 'ROW_2022', 'PLOT_ID'))
ft23$Exp_year <- 2023

# join data for experiment years
joined_ft <- full_join(ft22, ft23, by=c('Genotype', 'Condition', 'Replicate', 'FT'='FT', 'Generation', 'Plants'='Plot_Survival', 'Exp_year')) #%>% select(-'Plants')
# write out a copy of joined data for ease
write_delim(joined_ft, "FT_BOTH_YEARS.tsv", "\t")


### SIGNIFICANCE TESTS
# test relationship between flowering date (days between planting and spike emergence) and other variables
summary(aov(FT ~ Exp_year, joined_ft))
summary(aov(FT ~ Condition + Replicate + Generation + Plants + Exp_year + Condition*Generation, joined_ft))
summary(aov(FT ~ Condition*Exp_year + Replicate*Exp_year + Generation*Exp_year + Plants*Exp_year + Exp_year + Condition*Generation, joined_ft))


### CORRELATION
exp1 <- joined_ft %>% filter(Exp_year==2022) %>% filter(!is.na(FT)) %>% select(-c(Plants, Generation, Exp_year)) %>% group_by(Genotype, Condition) %>% summarise(geno_ft = mean(FT)) %>% pivot_wider(names_from=c(Condition), values_from=geno_ft)
exp2 <- joined_ft %>% filter(Exp_year==2023) %>% filter(!is.na(FT)) %>% select(-c(Plants, Generation, Exp_year)) %>% group_by(Genotype, Condition) %>% summarise(geno_ft = mean(FT)) %>% pivot_wider(names_from=c(Condition), values_from=geno_ft)

cor(exp1$mixed, exp1$single, use="complete.obs")
cor(exp2$mixed, exp2$single, use="complete.obs")
cor(exp1$single, exp2$single, use="complete.obs")
cor(exp1$mixed, exp2$mixed, use="complete.obs")




### Plotting
setwd("../results/")
## distribution (all)
png("distribution.png")
ggplot(joined_ft, aes(x=FT)) +
geom_density(alpha=0.5) +
theme_minimal() +
facet_wrap(~Exp_year, nrow=2)
dev.off()


## boxplot of experimental conditions
png("boxplot_condition.png")
ggplot(joined_ft, aes(y=FT, x=Condition, fill=Condition)) +
geom_boxplot() +
theme_minimal() +
facet_wrap(~Exp_year)
dev.off()

## distribution of experimental conditions
png("distribution_condition.png")
ggplot(joined_ft, aes(x=FT, group=Condition, color=Condition, fill=Condition)) +
geom_density(alpha=0.5) +
theme_minimal() +
facet_wrap(~Exp_year, nrow=2)
dev.off()



## boxplot of Generations
png("boxplot_generation.png")
joined_ft %>% ggplot( aes(y=FT, x=as.factor(Generation), group=as.factor(Generation), fill=as.factor(Generation))) +
geom_boxplot() +
theme_minimal() +
facet_wrap(~Exp_year)
dev.off()

## distribution of Generations
png("distribution_generation.png")
ggplot(joined_ft, aes(x=FT, group=as.factor(Generation), color=as.factor(Generation), fill=as.factor(Generation))) +
geom_density(alpha=0.5) +
theme_minimal() +
facet_wrap(~Exp_year, nrow=2)
dev.off()

## distribution of Generations - single (non-mixed) plots only
png("distribution_generation_single_condition.png")
joined_ft %>% filter(Condition == "single") %>% ggplot(aes(x=FT, group=as.factor(Generation), color=as.factor(Generation), fill=as.factor(Generation))) +
geom_density(alpha=0.5) +
theme_minimal() +
facet_wrap(~Exp_year, nrow=2)
dev.off()

## boxplot of Generations - single plots only
png("boxplot_generation_single_condition.png")
joined_ft %>% filter(Condition == "single") %>% ggplot( aes(y=FT, x=as.factor(Generation), group=as.factor(Generation), fill=as.factor(Generation))) +
geom_boxplot() +
theme_minimal() +
facet_wrap(~Exp_year)
dev.off()



## distribution of Generations - mixed plots only
png("distribution_generation_mixed_condition.png")
joined_ft %>% filter(Condition == "mixed") %>% ggplot(aes(x=FT, group=as.factor(Generation), color=as.factor(Generation), fill=as.factor(Generation))) +
geom_density(alpha=0.5) +
theme_minimal() +
facet_wrap(~Exp_year, nrow=2)
dev.off()

## boxplot of Generations - mixed plots only
png("boxplot_generation_mixed_condition.png")
joined_ft %>% filter(Condition == "mixed") %>% ggplot( aes(y=FT, x=as.factor(Generation), group=as.factor(Generation), fill=as.factor(Generation))) +
geom_boxplot() +
theme_minimal() +
facet_wrap(~Exp_year)
dev.off()



## boxplot of Replicates
png("boxplot_Replicates.png")
ggplot(joined_ft, aes(y=FT, x=Replicate, group=Replicate, fill=Replicate)) +
geom_boxplot() +
theme_minimal() +
facet_wrap(~Exp_year)
dev.off()

## distribution of Replicates
png("distribution_Replicates.png")
ggplot(joined_ft, aes(x=FT, group=Replicate, color=Replicate, fill=Replicate)) +
geom_density(alpha=0.5) +
theme_minimal() +
facet_wrap(~Exp_year, nrow=2)
dev.off()
