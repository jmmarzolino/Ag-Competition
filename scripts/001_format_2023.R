#!/usr/bin/env Rscript

#SBATCH --ntasks=1
#SBATCH --mem=50G
#SBATCH --time=02:00:00
#SBATCH --job-name="2023"
#SBATCH --output=/rhome/jmarz001/bigdata/Ag-Competition/scripts/001_format_2023.stdout
#SBATCH -p short


setwd("/rhome/jmarz001/bigdata/Ag-Competition/data")
library(tidyverse)

# load phenotyping data
competition <- read_delim("Field 2022-2023 Genotype List - Competition.csv", ",")
competition <- competition %>% select(-c(contains("albino"), 'Plot_Germination'))

ft <- read_delim("FT_DAYS_2022-2023.xlsx - FT_DAYS.csv", ",")
layout <- read_delim("FT_DAYS_2022-2023.xlsx - FIELD_LAYOUT.csv", ",")

merge <- full_join(ft, layout)
merge <- merge[order(merge$PLOT_ID),]
# subset merged data to only include CC II Competition part of experiment
df <- merge[grep("CC II", merge$PLANT_ID),]
# remove notes col
df <- df %>% select(-Notes)
# rows 19 and 30 are short, so that's where the replicates split
# separate duplicate cols into 2 data frames
rep1 <- df %>% filter(ROW <= 19)
rep2 <- df %>% filter(ROW > 19)

rep1$REP <- "rep 1"
rep2$REP <- "rep 2"

df3 <- bind_rows(rep1, rep2)

df3$FT_DAYS <- as.numeric(df3$FT_DAYS)
df3$PLANT_ID <- gsub("CC II (\\w+)", "\\1", df3$PLANT_ID)
df3$ROW <- df3$ROW - 8



df4 <- full_join(competition, df3, by=c('Genotypes'='PLANT_ID', 'replicate'='REP', 'Bed_2022'='ROW', 'PLOT_ID'))
## Remove rows without genotype
df4 <- df4 %>% filter(!is.na(Genotypes))
df4 <- df4 %>% select(-PLANT_ID)

# add generation column
df4$Generation <- 0
df4[grep("^1_", df4$Genotypes), which(colnames(df4)=="Generation")] <- 18
df4[grep("^2_", df4$Genotypes), which(colnames(df4)=="Generation")] <- 28
df4[grep("^3_", df4$Genotypes), which(colnames(df4)=="Generation")] <- 50
df4[grep("^7_", df4$Genotypes), which(colnames(df4)=="Generation")] <- 58

write_delim(df4, "FT_2023.tsv", "\t")




########### USEFUL FOR CURRENT FIELD WORK
# how many have a FT value already?
print(paste0("Percent that hasn't flowered: ", round((sum(is.na(df4$FT_DAYS)) / length(df4$FT_DAYS))*100)))

uncollected <- df4[which(is.na(df4$FT_DAYS)),]
print("number of plots for each generation that haven't flowered")
table(uncollected$Generation)

uncollected %>% filter(Condition=="single")

print("genotypes that haven't flowered")
unique(uncollected$Genotypes)

print("parents that haven't flowered")
uncollected %>% filter(Generation==0) %>% select(Genotypes) %>% table
uncollected %>% filter(Condition=="single") %>% filter(Generation==0) %>% select(Genotypes) %>% table
# 11 parents haven't flowered yet...(one almost definitely won't though)

# field positions that haven't flowered
not_flowered <- df[which(is.na(df$FT_DAYS)),]
write_delim(not_flowered, "not_flowered.tsv")

########################## REMOVE FOR LATER ANALYSIS




#### generation means and vars
df4 %>% group_by(Generation) %>% summarise(mean = mean(FT_DAYS, na.rm=T), variance = var(FT_DAYS, na.rm=T), n=n())
### conditions means and vars
df4 %>% group_by(Condition) %>% summarise(mean = mean(FT_DAYS, na.rm=T), variance = var(FT_DAYS, na.rm=T), n=n())
### replicates means and vars
df4 %>% group_by(replicate) %>% summarise(mean = mean(FT_DAYS, na.rm=T), variance = var(FT_DAYS, na.rm=T), n=n())



#### parent genotypes avg flowering date
report <- df4 %>% filter(Generation == 0) %>% filter(Condition=='single') %>% group_by(Genotypes) %>% summarise(mean = mean(FT_DAYS, na.rm=T), variance = var(FT_DAYS, na.rm=T), n=n())
write_delim(report, "avg_FT_DAYS_parents_2023.tsv", "\t")






### SIGNIFICANCE TESTS
# test relationship between flowering date (days between planting and spike emergence) and...

# replicate
x <-aov(FT_DAYS ~ replicate, df4)
summary(x)

# are replicates strongly correlated?
#rep1 <- df4 %>% select(Genotypes, replicate, FT_DAYS) %>% filter(replicate=="rep 1") #%>% filter(!is.na(FT_DAYS))
#rep2 <- df4 %>% select(Genotypes, replicate, FT_DAYS) %>% filter(replicate=="rep 2") #%>% filter(!is.na(FT_DAYS))

#cor(rep1$FT_DAYS, rep2$FT_DAYS)



# experimental group
summary(aov(FT_DAYS ~ Condition, df4))

# replicate and experimental group
summary(aov(FT_DAYS ~ replicate + Condition, df4))

# generation - parents vs progeny
#df4 %>% mutate(ParentOrProgeny = )

# generation - F0, 18, 28, 58...
summary(aov(FT_DAYS ~ Generation, df4))

# combinations of above?
summary(aov(FT_DAYS ~ Condition + replicate + Generation + Plot_Survival, df4))

summary(aov(FT_DAYS ~ Condition*Generation, df4))

summary(aov(FT_DAYS ~ Condition + replicate + Generation + Plot_Survival, df4))

summary(aov(FT_DAYS ~ replicate + Generation + Condition + Plot_Survival + Generation*Plot_Survival + replicate*Plot_Survival + Generation*Plot_Survival + Condition*Generation, df4))
