#!/usr/bin/env Rscript

#SBATCH --ntasks=1
#SBATCH --mem=50G
#SBATCH --time=02:00:00
#SBATCH --job-name="Plot FT"
#SBATCH --output=/rhome/jmarz001/bigdata/Ag-Competition/001_format_plot_floweringtime.stdout
#SBATCH -p short


setwd("/rhome/jmarz001/bigdata/Ag-Competition")
library(tidyverse)

# load phenotyping data
df <- read_delim("Phenotypes.csv", ",")
ft <- read_delim("FT_2022_2023.csv", ",")
layout <- read_delim("FT_2022-2023.xlsx - FIELD_LAYOUT.csv", ",")

# remove comments and count columns
df <- df[, -(17:19)]

# separate duplicate cols into 2 data frames
df1 <- df[, 1:8]
df2 <- df[, 9:16]

# copy correct col names onto second data frame
colnames(df2) <- colnames(df1)
df3 <- bind_rows(df1, df2)

# remove non-numeric date entries
#df3 <- df3[-which(df3$`Flowering Date` == "x"),]
#df3 <- df3[-which(df3$`Flowering Date` == "X"),]
df3$`Flowering Date` <- as.numeric(df3$`Flowering Date`)

# replace colnames with code-friendly versions
colnames(df3) <- c("Genotype", "Plants", "Condition", "Replicate", "BED_2021", "ROW_2021", "FT", "Notes")


## Remove rows without genotype
df3 <- df3 %>% filter(!is.na(Genotype))

#### Add Generation
df3$Genotype <- str_replace(df3$Genotype, "-", "_")
df3$Genotype <- str_replace(df3$Genotype, "-", "_")

df3$Generation <- str_replace(df3$Genotype, "(\\d)_\\d+", "\\1")
df3$Generation <- str_replace(df3$Generation, "(\\d)_\\d+", "\\1")
#str_split_fixed(df3$Genotype, "_", 3)

df3$Generation <- as.numeric(df3$Generation)
df3[which(df3$Generation > 8), which(colnames(df3) == "Generation")] <- 0



#### generation means and vars
df3 %>% group_by(Generation) %>% summarise(mean = mean(FT, na.rm=T), variance = var(FT, na.rm=T), n=n())
### conditions means and vars
df3 %>% group_by(Condition) %>% summarise(mean = mean(FT, na.rm=T), variance = var(FT, na.rm=T), n=n())
### Replicates means and vars
df3 %>% group_by(Replicate) %>% summarise(mean = mean(FT, na.rm=T), variance = var(FT, na.rm=T), n=n())



#### parent genotypes avg flowering date
report <- df3 %>% filter(Generation == 0) %>% filter(Condition=="single") %>% group_by(Genotype) %>% summarise(mean = mean(FT, na.rm=T), variance = var(FT, na.rm=T), n=n())
write_delim(report, "avg_FT_parents.tsv", "\t")



### Plotting
## boxplot of experimental conditions
png("boxplot_condition.png")
ggplot(df3, aes(y=FT, x=Condition, fill=Condition)) +
geom_boxplot() +
theme_minimal()
dev.off()

## distribution of experimental conditions
png("distribution_condition.png")
ggplot(df3, aes(x=FT, group=Condition, color=Condition, fill=Condition)) +
geom_density(alpha=0.5) +
theme_minimal()
dev.off()



## boxplot of Generations
png("boxplot_generation.png")
df3 %>% ggplot( aes(y=FT, x=as.factor(Generation), group=as.factor(Generation), fill=as.factor(Generation))) +
geom_boxplot() +
theme_minimal()
dev.off()

## distribution of Generations
png("distribution_generation.png")
ggplot(df3, aes(x=FT, group=as.factor(Generation), color=as.factor(Generation), fill=as.factor(Generation))) +
geom_density(alpha=0.5) +
theme_minimal()
dev.off()

## distribution of Generations - single (non-mixed) plots only
png("distribution_generation_single_condition.png")
df3 %>% filter(Condition == "single") %>% ggplot(aes(x=FT, group=as.factor(Generation), color=as.factor(Generation), fill=as.factor(Generation))) +
geom_density(alpha=0.5) +
theme_minimal()
dev.off()

## boxplot of Generations - single plots only
png("boxplot_generation_single_condition.png")
df3 %>% filter(Condition == "single") %>% ggplot( aes(y=FT, x=as.factor(Generation), group=as.factor(Generation), fill=as.factor(Generation))) +
geom_boxplot() +
theme_minimal()
dev.off()


## boxplot of Replicates
png("boxplot_Replicates.png")
ggplot(df3, aes(y=FT, x=Replicate, group=Replicate, fill=Replicate)) +
geom_boxplot() +
theme_minimal()
dev.off()

## distribution of Replicates
png("distribution_Replicates.png")
ggplot(df3, aes(x=FT, group=Replicate, color=Replicate, fill=Replicate)) +
geom_density(alpha=0.5) +
theme_minimal()
dev.off()


### SIGNIFICANCE TESTS
# test relationship between flowering date (days between planting and spike emergence) and...

# Replicate
x <-aov(FT ~ Replicate, df3)
summary(x)

# are Replicates strongly correlated?
#rep1 <- df3 %>% select(Genotype, Replicate, FT) %>% filter(Replicate=="rep 1") #%>% filter(!is.na(FT))
#rep2 <- df3 %>% select(Genotype, Replicate, FT) %>% filter(Replicate=="rep 2") #%>% filter(!is.na(FT))

#cor(rep1$FT, rep2$FT)



# experimental group
x <-aov(FT ~ Condition, df3)
summary(x)

# Replicate and experimental group
x <-aov(FT ~ Replicate + Condition, df3)
summary(x)

# generation - parents vs progeny
#df3 %>% mutate(ParentOrProgeny = )

# generation - F0, 18, 28, 58...
x <-aov(FT ~ Generation, df3)
summary(x)

# combinations of above?
x <-aov(FT ~ Condition + Replicate + Generation + Plants, df3)
summary(x)

x <-aov(FT ~ Condition*Generation, df3)
summary(x)

#library(lme4)
#install_packages("lmeTest")
#library(lmeTest)
#x <-aov(lmer(FT ~ Condition + (1|Generation:Genotype) + (1|Replicate), df3))
#x <-aov(lmer(FT ~ Condition + Generation + (1 + Generation|Genotype) + (1|Replicate)), df3) # try to plot this line!
# Fixed: Condition, Generation
# Random: Genotype, Replicate
# relationship between Generation and genotype?
# relationship between Generation and Condition...


#summary(x)
#anova(lmer()) # anova for fixed effects
#summary(aov())
#ranova(lmer()) # for random effects
#t.test(extra ~ group, data = sleep)
