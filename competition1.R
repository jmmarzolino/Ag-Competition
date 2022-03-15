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
ggplot(df3, aes(y=Flowering_Date, x=Generation, group=Generation, fill=Generation)) +
geom_boxplot() +
theme_minimal()
dev.off()

## distribution of Generations
png("distribution_generation.png")
ggplot(df3, aes(x=Flowering_Date, group=Generation, color=Generation, fill=Condition)) +
geom_density(alpha=0.5) +
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

summarise(iris, sd(Sepal.Length, na.rm = T), groups=Species)
