#!/usr/bin/env Rscript
#SBATCH --mem=50G
#SBATCH --time=02:00:00
#SBATCH --output=/rhome/jmarz001/bigdata/Ag-Competition/scripts/003_change_over_gens.stdout
#SBATCH -p short

setwd("/rhome/jmarz001/bigdata/Ag-Competition/data")
library(tidyverse)

source("../scripts/CUSTOM_FNS.R")

## Remove rows without genotype
df4 <- df4 %>% filter(!is.na(Genotype))
df4 <- df4 %>% select(-PLANT_ID)

df1 <- read_delim("PHENOS_2021_2022.tsv")
df2 <- read_delim("PHENOS_2022_2023.tsv")



## Add a Generation Col
df$Generation <- 0
df[grep("^1$", df$FAM_ID), which(colnames(df)=="Generation")] <- 18
df[grep("^2$", df$FAM_ID), which(colnames(df)=="Generation")] <- 28
df[grep("^3$", df$FAM_ID), which(colnames(df)=="Generation")] <- 50
df[grep("^7$", df$FAM_ID), which(colnames(df)=="Generation")] <- 58

#### generation means and vars
df %>% select(c(Generation, Condition, SEED_WEIGHT_100)) %>% group_by(Generation, Condition) %>% summarise(avg_seed_weight = mean(SEED_WEIGHT_100, na.rm=T), seed_weight_var = var(SEED_WEIGHT_100, na.rm=T)) %>% pivot_wider(names_from=Condition, values_from=c(avg_seed_weight, seed_weight_var))

df$FT <- as.numeric(df$FT)
df %>% select(c(Generation, Condition, FT)) %>% group_by(Generation, Condition) %>% summarise(avg_ft = mean(FT, na.rm=T), ft_var = var(FT, na.rm=T)) %>% pivot_wider(names_from=Condition, values_from=c(avg_ft, ft_var))



#### parent genotypes avg flowering date
report <- df %>% filter(Generation == 0) %>% filter(Condition=="single") %>% group_by(Genotype) %>% summarise(mean = mean(SEED_WEIGHT_100, na.rm=T), variance = var(SEED_WEIGHT_100, na.rm=T), n=n())
write_delim(report, "avg_100seed_weight_parents.tsv", "\t")




#### generation means and vars
df4 %>% group_by(Generation) %>% summarise(mean = mean(FT, na.rm=T), variance = var(FT, na.rm=T), n=n())
### conditions means and vars
df4 %>% group_by(Condition) %>% summarise(mean = mean(FT, na.rm=T), variance = var(FT, na.rm=T), n=n())
### Replicates means and vars
df4 %>% group_by(Replicate) %>% summarise(mean = mean(FT, na.rm=T), variance = var(FT, na.rm=T), n=n())



#### parent genotypes avg flowering date
report <- df4 %>% filter(Generation == 0) %>% filter(Condition=="single") %>% group_by(Genotype) %>% summarise(mean = mean(FT, na.rm=T), variance = var(FT, na.rm=T), n=n())
write_delim(report, "avg_FT_parents_2023.tsv", "\t")






### SIGNIFICANCE TESTS
# test relationship between flowering date (days between planting and spike emergence) and...

# Replicate
x <-aov(FT ~ Replicate, df4)
summary(x)

# are Replicates strongly correlated?
#rep1 <- df4 %>% select(Genotype, Replicate, FT) %>% filter(Replicate=="rep 1") #%>% filter(!is.na(FT))
#rep2 <- df4 %>% select(Genotype, Replicate, FT) %>% filter(Replicate=="rep 2") #%>% filter(!is.na(FT))

#cor(rep1$FT, rep2$FT)



# experimental group
summary(aov(FT ~ Condition, df4))

# Replicate and experimental group
summary(aov(FT ~ Replicate + Condition, df4))

# generation - parents vs progeny
#df4 %>% mutate(ParentOrProgeny = )

# generation - F0, 18, 28, 58...
summary(aov(FT ~ Generation, df4))

# combinations of above?
summary(aov(FT ~ Condition + Replicate + Generation + Plot_Survival, df4))

summary(aov(FT ~ Condition*Generation, df4))

summary(aov(FT ~ Condition + Replicate + Generation + Plot_Survival, df4))

summary(aov(FT ~ Replicate + Generation + Condition + Plot_Survival + Generation*Plot_Survival + Replicate*Plot_Survival + Generation*Plot_Survival + Condition*Generation, df4))



#!/usr/bin/env Rscript

#SBATCH --ntasks=1
#SBATCH --mem=50G
#SBATCH --time=02:00:00
#SBATCH --job-name="2022"
#SBATCH --output=/rhome/jmarz001/bigdata/Ag-Competition/scripts/001_format_2022_phenotypes.stdout
#SBATCH -p short


setwd("/rhome/jmarz001/bigdata/Ag-Competition/data")

# load phenotyping data
ft <- read_delim("FT_2021_2022.tsv")
#Genotype        Plants  Condition       Replicate       BED_2021        ROW_2021        FT
sw <- read_delim("SEED_WEIGHTS_2021_2022.csv")



# remove rows without genotype
#df <- df %>% filter(!is.na(Genotype))
# keep na genotypes for now to keep field bed/row position filled







#### generation means and vars
df3 %>% group_by(Generation) %>% summarise(mean = mean(FT, na.rm=T), variance = var(FT, na.rm=T), n=n())
### conditions means and vars
df3 %>% group_by(Condition) %>% summarise(mean = mean(FT, na.rm=T), variance = var(FT, na.rm=T), n=n())
### Replicates means and vars
df3 %>% group_by(Replicate) %>% summarise(mean = mean(FT, na.rm=T), variance = var(FT, na.rm=T), n=n())



#### parent genotypes avg flowering date
report <- df3 %>% filter(Generation == 0) %>% filter(Condition=="single") %>% group_by(Genotype) %>% summarise(mean = mean(FT, na.rm=T), variance = var(FT, na.rm=T), n=n())
write_delim(report, "avg_FT_parents.tsv", "\t")




### SIGNIFICANCE TESTS
# test relationship between flowering date (days between planting and spike emergence) and...

# Replicate
summary(aov(FT ~ Replicate, df3))

# are Replicates strongly correlated?
#rep1 <- df3 %>% select(Genotype, Replicate, FT) %>% filter(Replicate=="rep 1") #%>% filter(!is.na(FT))
#rep2 <- df3 %>% select(Genotype, Replicate, FT) %>% filter(Replicate=="rep 2") #%>% filter(!is.na(FT))

#cor(rep1$FT, rep2$FT)



# experimental group
summary(aov(FT ~ Condition, df3))

# Replicate and experimental group
summary(aov(FT ~ Replicate + Condition, df3))

# generation - parents vs progeny
#df3 %>% mutate(ParentOrProgeny = )

# generation - F0, 18, 28, 58...
summary(aov(FT ~ Generation, df3))

# combinations of above?
summary(aov(FT ~ Condition + Replicate + Generation + Plants, df3))

summary(aov(FT ~ Condition*Generation, df3))

summary(aov(FT ~ Replicate + Generation + Condition + Plants + Generation*Plants + Replicate*Plants + Generation*Plants + Condition*Generation, df3))












setwd("/rhome/jmarz001/bigdata/Ag-Competition")
library(tidyverse)

# load phenotyping data
df <- read_delim("Phenotypes.csv", ",")
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
