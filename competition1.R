#!/usr/bin/env Rscript

#SBATCH --ntasks=1
#SBATCH --mem=50G
#SBATCH --time=02:00:00
#SBATCH --job-name="Ag-Competition 1"
#SBATCH --output=/rhome/jmarz001/bigdata/Ag-Competition/competition1.stdout
#SBATCH -p short


setwd("/rhome/jmarz001/bigdata/Ag-Competition")
library(tidyverse)

# load phenotyping data
df <- read_delim("Phenotyping Sheets - germination, survival, flowering time, height - CC II Competition Phenotyping - Formatted.csv", ",")
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
colnames(df3) <- c("Genotypes", "number_of_plants", "Condition", "replicate", "2021BED", "2021ROW", "Flowering_Date", "Notes")


## Remove rows without genotype
df3 <- df3 %>% filter(!is.na(Genotypes))

#### Add Generation
df3$Genotypes <- str_replace(df3$Genotypes, "-", "_")
df3$Genotypes <- str_replace(df3$Genotypes, "-", "_")

df3$Generation <- str_replace(df3$Genotypes, "(\\d)_\\d+", "\\1")
df3$Generation <- str_replace(df3$Generation, "(\\d)_\\d+", "\\1")
#str_split_fixed(df3$Genotypes, "_", 3)

df3$Generation <- as.numeric(df3$Generation)
df3[which(df3$Generation > 8), which(colnames(df3) == "Generation")] <- 0



#### generation means and vars
df3 %>% group_by(Generation) %>% summarise(mean = mean(Flowering_Date, na.rm=T), variance = var(Flowering_Date, na.rm=T), n=n())
### conditions means and vars
df3 %>% group_by(Condition) %>% summarise(mean = mean(Flowering_Date, na.rm=T), variance = var(Flowering_Date, na.rm=T), n=n())
### replicates means and vars
df3 %>% group_by(replicate) %>% summarise(mean = mean(Flowering_Date, na.rm=T), variance = var(Flowering_Date, na.rm=T), n=n())



### Plotting
## boxplot of experimental conditions
png("boxplot_condition.png")
ggplot(df3, aes(y=Flowering_Date, x=Condition, fill=Condition)) +
geom_boxplot() +
theme_minimal()
dev.off()

## distribution of experimental conditions
png("distribution_condition.png")
ggplot(df3, aes(x=Flowering_Date, group=Condition, color=Condition, fill=Condition)) +
geom_density(alpha=0.5) +
theme_minimal()
dev.off()



## boxplot of Generations
png("boxplot_generation.png")
df3 %>% ggplot( aes(y=Flowering_Date, x=as.factor(Generation), group=as.factor(Generation), fill=as.factor(Generation))) +
geom_boxplot() +
theme_minimal()
dev.off()

## distribution of Generations
png("distribution_generation.png")
ggplot(df3, aes(x=Flowering_Date, group=as.factor(Generation), color=as.factor(Generation), fill=as.factor(Generation))) +
geom_density(alpha=0.5) +
theme_minimal()
dev.off()

## distribution of Generations - single (non-mixed) plots only
png("distribution_generation_single_condition.png")
df3 %>% filter(Condition == "single") %>% ggplot(aes(x=Flowering_Date, group=as.factor(Generation), color=as.factor(Generation), fill=as.factor(Generation))) +
geom_density(alpha=0.5) +
theme_minimal()
dev.off()

## boxplot of Generations - single plots only
png("boxplot_generation_single_condition.png")
df3 %>% filter(Condition == "single") %>% ggplot( aes(y=Flowering_Date, x=as.factor(Generation), group=as.factor(Generation), fill=as.factor(Generation))) +
geom_boxplot() +
theme_minimal()
dev.off()


## boxplot of replicates
png("boxplot_replicates.png")
ggplot(df3, aes(y=Flowering_Date, x=replicate, group=replicate, fill=replicate)) +
geom_boxplot() +
theme_minimal()
dev.off()

## distribution of replicates
png("distribution_replicates.png")
ggplot(df3, aes(x=Flowering_Date, group=replicate, color=replicate, fill=replicate)) +
geom_density(alpha=0.5) +
theme_minimal()
dev.off()


### SIGNIFICANCE TESTS
# test relationship between flowering date (days between planting and spike emergence) and...

# replicate
x <-aov(Flowering_Date ~ replicate, df3)
summary(x)

# are replicates strongly correlated?
#rep1 <- df3 %>% select(Genotypes, replicate, Flowering_Date) %>% filter(replicate=="rep 1") #%>% filter(!is.na(Flowering_Date))
#rep2 <- df3 %>% select(Genotypes, replicate, Flowering_Date) %>% filter(replicate=="rep 2") #%>% filter(!is.na(Flowering_Date))

#cor(rep1$Flowering_Date, rep2$Flowering_Date)



# experimental group
x <-aov(Flowering_Date ~ Condition, df3)
summary(x)

# replicate and experimental group
x <-aov(Flowering_Date ~ replicate + Condition, df3)
summary(x)

# generation - parents vs progeny
#df3 %>% mutate(ParentOrProgeny = )

# generation - F0, 18, 28, 58...
x <-aov(Flowering_Date ~ Generation, df3)
summary(x)

# combinations of above?
x <-aov(Flowering_Date ~ Condition + replicate + Generation + number_of_plants, df3)
summary(x)

x <-aov(Flowering_Date ~ Condition*Generation, df3)
summary(x)

#library(lme4)
#install_packages("lmeTest")
#library(lmeTest)
#x <-aov(lmer(Flowering_Date ~ Condition + (1|Generation:Genotypes) + (1|replicate), df3))
#x <-aov(lmer(Flowering_Date ~ Condition + Generation + (1 + Generation|Genotypes) + (1|replicate)), df3) # try to plot this line!
# Fixed: Condition, Generation
# Random: Genotype, replicate
# relationship between Generation and genotype?
# relationship between Generation and Condition...


#summary(x)
#anova(lmer()) # anova for fixed effects
#summary(aov())
#ranova(lmer()) # for random effects
#t.test(extra ~ group, data = sleep)
