setwd("/rhome/jmarz001/bigdata/Competition")

df <- read_delim("Phenotyping Sheets - germination, survival, flowering time, height - CC II Competition Phenotyping - Formatted.csv", ",")

df <- df[, -(17:19)]

df1 <- df[, 1:8]
df2 <- df[, 9:16]


colnames(df2) <- colnames(df1)
df3 <- bind_rows(df1, df2)


#df3 <- df3[-which(df3$`Flowering Date` == "x"),]
#df3 <- df3[-which(df3$`Flowering Date` == "X"),]
df3$`Flowering Date` <- as.numeric(df3$`Flowering Date`)
colnames(df3) <- c("Genotypes", "number_of_plants", "Condition", "replicate", "2021BED", "2021ROW", "Flowering_Date", "Notes")


t.test(extra ~ group, data = sleep)
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

## boxplot of experimental conditions
png("boxplot.png")
ggplot(df3, aes(y=Flowering_Date, x=Condition, fill=Condition)) +
geom_boxplot() +
theme_minimal()
dev.off()

## distribution
png("distribution.png")
ggplot(df3, aes(x=Flowering_Date, group=Condition, color=Condition, fill=Condition)) +
geom_density(alpha=0.5) +
theme_minimal()
dev.off()


x <-aov(Flowering_Date ~ replicate, df3)
summary(x)

 df3 %>% group_by(replicate) %>% summarise(mean(Flowering_Date, na.rm=T))

summarise(iris, sd(Sepal.Length, na.rm = T), groups=Species)
