#!/usr/bin/env Rscript

#SBATCH --ntasks=1
#SBATCH --mem=50G
#SBATCH --time=02:00:00
#SBATCH --job-name="2022"
#SBATCH --output=/rhome/jmarz001/bigdata/Ag-Competition/scripts/001_format_2022.stdout
#SBATCH -p short


setwd("/rhome/jmarz001/bigdata/Ag-Competition/data")
library(tidyverse)

# load phenotyping data
df <- read_delim("Phenotyping Sheets - germination, survival, flowering time, height - CC II Competition Phenotyping - Formatted.csv", ",")

# remove comments and count columns
df <- df[, -(17:19)]

# separate duplicate cols into 2 data frames
df1 <- df[, 1:8]
df2 <- df[, 9:16]

# replace colnames with code-friendly versions
colnames(df2) <- c("Genotypes", "number_of_plants", "Condition", "replicate", "2021BED", "2021ROW", "Flowering_Date", "Notes")
colnames(df1) <- c("Genotypes", "number_of_plants", "Condition", "replicate", "2021BED", "2021ROW", "Flowering_Date", "Notes")
df3 <- bind_rows(df1, df2)

# remove non-numeric date entries
#df3 <- df3[-which(df3$`Flowering Date` == "x"),]
#df3 <- df3[-which(df3$`Flowering Date` == "X"),]
df3$Flowering_Date <- as.numeric(df3$Flowering_Date)
#colnames(df3) <- c("Genotypes", "number_of_plants", "Condition", "replicate", "2021BED", "2021ROW", "Flowering_Date", "Notes")
df3 <- df3 %>% select(-Notes)

## Remove rows without genotype
df3 <- df3 %>% filter(!is.na(Genotypes))

# add generation column
df3$Genotypes <- str_replace(df3$Genotypes, "-", "_")
#df3$Generation <- str_replace(df3$Genotypes, "(\\d)_\\d+", "\\1")
#df3$Generation <- str_replace(df3$Generation, "(\\d)_\\d+", "\\1")
#str_split_fixed(df3$Genotypes, "_", 3)
#df3[which(df3$Generation > 8), which(colnames(df3) == "Generation")] <- 0
df3$Generation <- 0
df3[grep("^1_", df3$Genotypes), which(colnames(df3)=="Generation")] <- 18
df3[grep("^2_", df3$Genotypes), which(colnames(df3)=="Generation")] <- 28
df3[grep("^3_", df3$Genotypes), which(colnames(df3)=="Generation")] <- 50
df3[grep("^7_", df3$Genotypes), which(colnames(df3)=="Generation")] <- 58

# simplify replicate col to a number
df3$replicate <- gsub("rep (\\d)", "\\1", df3$replicate)

# remove extra numbers and characters from germination/survival count
df3$number_of_plants <- gsub("(\\d+)-\\d+", "\\1", df3$number_of_plants)
df3$number_of_plants <- gsub(".* (\\d+)", "\\1", df3$number_of_plants)
df3$number_of_plants <- as.numeric(df3$number_of_plants)

write_delim(df3, "FT_2022.tsv", "\t")






#### generation means and vars
df3 %>% group_by(Generation) %>% summarise(mean = mean(Flowering_Date, na.rm=T), variance = var(Flowering_Date, na.rm=T), n=n())
### conditions means and vars
df3 %>% group_by(Condition) %>% summarise(mean = mean(Flowering_Date, na.rm=T), variance = var(Flowering_Date, na.rm=T), n=n())
### replicates means and vars
df3 %>% group_by(replicate) %>% summarise(mean = mean(Flowering_Date, na.rm=T), variance = var(Flowering_Date, na.rm=T), n=n())



#### parent genotypes avg flowering date
report <- df3 %>% filter(Generation == 0) %>% filter(Condition=='single') %>% group_by(Genotypes) %>% summarise(mean = mean(Flowering_Date, na.rm=T), variance = var(Flowering_Date, na.rm=T), n=n())
write_delim(report, "avg_flowering_date_parents.tsv", "\t")




### SIGNIFICANCE TESTS
# test relationship between flowering date (days between planting and spike emergence) and...

# replicate
summary(aov(Flowering_Date ~ replicate, df3))

# are replicates strongly correlated?
#rep1 <- df3 %>% select(Genotypes, replicate, Flowering_Date) %>% filter(replicate=="rep 1") #%>% filter(!is.na(Flowering_Date))
#rep2 <- df3 %>% select(Genotypes, replicate, Flowering_Date) %>% filter(replicate=="rep 2") #%>% filter(!is.na(Flowering_Date))

#cor(rep1$Flowering_Date, rep2$Flowering_Date)



# experimental group
summary(aov(Flowering_Date ~ Condition, df3))

# replicate and experimental group
summary(aov(Flowering_Date ~ replicate + Condition, df3))

# generation - parents vs progeny
#df3 %>% mutate(ParentOrProgeny = )

# generation - F0, 18, 28, 58...
summary(aov(Flowering_Date ~ Generation, df3))

# combinations of above?
summary(aov(Flowering_Date ~ Condition + replicate + Generation + number_of_plants, df3))

summary(aov(Flowering_Date ~ Condition*Generation, df3))

summary(aov(Flowering_Date ~ replicate + Generation + Condition + number_of_plants + Generation*number_of_plants + replicate*number_of_plants + Generation*number_of_plants + Condition*Generation, df3))
