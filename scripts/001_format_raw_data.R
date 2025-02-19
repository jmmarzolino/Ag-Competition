#!/usr/bin/env Rscript

#SBATCH --mem=30G
#SBATCH --time=02:00:00
#SBATCH --output=/rhome/jmarz001/bigdata/Ag-Competition/scripts/001_format_raw_data.stdout
#SBATCH -p koeniglab


# format raw data sources, join the two years of data, and check for basic errors or missing data
library(tidyverse)
library(ggpubr)

setwd("/bigdata/koeniglab/jmarz001/Ag-Competition/data")
source("/bigdata/koeniglab/jmarz001/Ag-Competition/scripts/CUSTOM_FNS.R")


#### Data from 2021-2022
sw <- read_delim("SEED_WEIGHTS_2021_2022.tsv")
# seed weight 2021-2022 already contains flowering time data

# average extra replicates for genotypes 1_6 and 7_69
tmp <- sw %>% filter(Genotype %in% c("1_6", "7_69"))
sw2 <- sw %>% filter(!(Genotype %in% c("1_6", "7_69")))

tmp[which(tmp$Replicate == "rep 3"), 4] <- "rep 1"
tmp[which(tmp$Replicate == "rep 4"), 4] <- "rep 2"

tmp2 <- tmp %>%
        group_by(Genotype, Condition, Replicate) %>%
        summarise("Plants" = mean(Plants), "TOTAL_MASS" = mean(TOTAL_MASS, na.rm = TRUE), "SEED_WEIGHT_100" = mean(SEED_WEIGHT_100, na.rm = TRUE), "FT" = mean(FT, na.rm=TRUE))

ord <- tmp %>% select(c(Genotype, Condition, Replicate,  BED_2021, ROW_2021)) 
ord2 <- ord[order(ord$Genotype, ord$Condition, ord$Replicate),][c(1,3,5,7,9,11,13,15),]

avgd <- full_join(tmp2, ord2, by = c("Genotype", "Condition", "Replicate"))
sw3 <- full_join(sw2, avgd, by = c("TOTAL_MASS", "SEED_WEIGHT_100", "FT", "Genotype", "Condition", "Replicate", "BED_2021", "ROW_2021", "Plants"))

# add year col
sw3$Exp_year <- 2022



#### Data from 2022-2023
## combine flowering time, seed weights, line IDs and plant experiment info
sw4 <- read_delim("SEED_WEIGHTS_2022_2023.tsv")
ft <- read_delim("FT_2022_2023.tsv")
geno_list <- read_delim("Genotype_List_2023_2023.tsv")

# join data sets
ft <- full_join(ft, geno_list, by = c("PLOT_ID", "ROW", "PLANT_ID"))
# join the seed weights and remove redundant column
sw5 <- full_join(sw4, ft, by = c("PLOT_ID")) %>% select(-"PLANT_ID")

## Tare weights
# Subtract the avg brown bag weight and
# average coin envelope weight from phenotypes
sw4$TOTAL_MASS <- sw4$TOTAL_MASS - 11.24
sw4$SEED_WEIGHT_100 <- sw4$SEED_WEIGHT_100 - 1.61

# average extra replicates for genotypes 1_6 and 7_69
tmp <- sw5 %>% filter(Genotype %in% c("1_6", "7_69"))
sw6 <- sw5 %>% filter(!(Genotype %in% c("1_6", "7_69")))

tmp2 <- tmp %>%
        group_by(Genotype, Condition, Replicate) %>%
        summarise("Plants" = mean(Plants), "TOTAL_MASS" = mean(TOTAL_MASS, na.rm = TRUE), "SEED_WEIGHT_100" = mean(SEED_WEIGHT_100, na.rm = TRUE), "FT" = mean(FT, na.rm=TRUE))

ord <- tmp %>% select(c(Genotype, PLOT_ID, Condition, Replicate, ROW, BED_2022, ROW_2022)) 
ord2 <- ord[order(ord$Genotype, ord$Condition, ord$Replicate),][c(1,3,5,7,9,11,13,15),]


avgd <- full_join(tmp2, ord2, by = c("Genotype", "Condition", "Replicate"))
sw7 <- full_join(sw6, avgd, by = c("PLOT_ID", "TOTAL_MASS", "SEED_WEIGHT_100", "FT", "ROW", "Genotype", "Condition", "Replicate", "BED_2022", "ROW_2022", "Plants"))

#sw4 <- sw4 %>%
 #       group_by(Genotype, Condition, Replicate) %>%
  #      summarise(across(where(is.numeric), \(x) mean(x, na.rm = TRUE)))

# add year col
sw7$Exp_year <- 2023



### Join years of data
sw_join <- full_join(sw3, sw7, by = c("Genotype", "Plants", "Condition", "Replicate", "FT", "TOTAL_MASS", "SEED_WEIGHT_100", "Exp_year", "BED_2021" = "BED_2022", "ROW_2021" = "ROW_2022")) 



### QC checks
# check for duplicated lines
which(table(sw_join$PLOT_ID)>1)
which(table(sw_join$Genotype)>16)
which(table(sw_join$Genotype)>8)
sw_join <- sw_join %>% select(-c("PLOT_ID", "ROW"))
colnames(sw_join) <- c("Genotype", "Plants", "Condition", "Replicate", "BED", "ROW", "FT", "TOTAL_MASS", "SEED_WEIGHT_100", "Exp_year")

# change replicate from character to number
sw_join$Replicate <- as.numeric(gsub("rep (\\d)", "\\1", sw_join$Replicate))


# verify phenotype measurements
# are there any phenotyped lines with plant count of 0?
sw_join[which(sw_join$Plants==0 & !is.na(sw_join$TOTAL_MASS)), ]
sw_join[which(sw_join$Plants==0 & !is.na(sw_join$SEED_WEIGHT_100)), ]
sw_join[which(sw_join$Plants==0 & !is.na(sw_join$FT)), ]



# how many empty rows are there?
# empty genotypes?
sw_join[which(is.na(sw_join$Genotype)),]
# 4 empty rows in 2022 data
# remove those empty rows
sw_join <- sw_join %>% filter(!is.na(Genotype))

# check out plots w 0 plants
(x <- sw_join[which(sw_join$Plants == 0),])
table(x$Genotype)
# all Genotype 2_168 plots have 0 plants
# line 2_168 were all albino, though this doesn't explain the observation in 'mixed' plots
# filter all empty plot rows
sw_join <- sw_join %>% filter(Plants != 0)


# empty conditions, replicates, etc
sw_join[which(is.na(sw_join$Condition)),]
sw_join[which(is.na(sw_join$Replicate)),]
sw_join[which(is.na(sw_join$Exp_year)),]
sw_join[which(is.na(sw_join$BED)),]
sw_join[which(is.na(sw_join$ROW)),]


# any other missing phenotypes?
sw_join[which(is.na(sw_join$FT)),] # missing FT is fine
sw_join[which(is.na(sw_join$TOTAL_MASS)),]
sw_join[which(is.na(sw_join$SEED_WEIGHT_100)),]
sw_join[which(!is.na(sw_join$TOTAL_MASS) & is.na(sw_join$SEED_WEIGHT_100)),] 
# plots with harvested seed but no 100-seed-weight

sw_join[which(is.na(sw_join$TOTAL_MASS) & is.na(sw_join$SEED_WEIGHT_100)),] 
# rows missing both total mass and 100-seed-weight still have FT records
# keeping all records in basic joined data, filter subsequently as needed


# check range of all columns
sw_join %>% reframe(across(where(is.numeric), \(x) range(x, na.rm=TRUE)))

sw_join %>% group_by(Exp_year, Condition) %>% reframe(across(where(is.numeric), \(x) range(x, na.rm=TRUE)))


write_delim(sw_join, "JOINED_PHENOTYPES.tsv", "\t")


## additionally, print the statistics for how many raw values there are, how many values were measured for each trait across years
df <- fread("JOINED_PHENOTYPES.tsv")
df <- df %>% filter(Condition=="single")


df1 <- df %>% filter(Exp_year==2022)
df2 <- df %>% filter(Exp_year==2023)

print('number of rows per year')
nrow(df1)
nrow(df2)

print("amount of missing data for the 3 measured traits")
df %>% group_by(Exp_year) %>% reframe(across(c(FT, TOTAL_MASS, SEED_WEIGHT_100), \(x) sum(is.na(x))))

print("percent of missing data per year")
df1 %>% summarise(across(c(FT, TOTAL_MASS, SEED_WEIGHT_100), \(x) (sum(is.na(x))/508)*100))

df2 %>% summarise(across(c(FT, TOTAL_MASS, SEED_WEIGHT_100), \(x) (sum(is.na(x))/508)*100))