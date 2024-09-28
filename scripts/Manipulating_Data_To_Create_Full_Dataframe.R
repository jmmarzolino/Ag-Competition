library(tidyverse)
library(readr)
library(dplyr)
library(ggpubr)
library(ggplot2)
library(googlesheets4)
library(tidyr)

setwd("/bigdata/koeniglab/jmarz001/Ag-Competition/data")
### Load Data

Seed_weights_2022_2023 <- read_delim("Seed Weights - Field 2023.csv") %>% select(-c(Date, Notes))
FT_2022_2023 <- read_delim("FT_DAYS_2022-2023.xlsx - FT_DAYS.csv")
FT_2023 <- read_delim('FT_2023.tsv')
FT_2022 <- read_delim('FT_2022.tsv')
Genotype_List_2022_2023 <- read_delim("Field 2022-2023 Genotype List - Competition.csv")
Haplo_raw <- read_delim("Competition Lines - Sheet1 - Working - Competition Lines - Sheet1.csv")
Seed_weights_2021_2022 <- read_delim("Seed Weights 2021-2022 - Sheet1.csv") %>% select(c(Genotypes, germinated, Condition, replicate,  `2021BED`, `2021ROW`, Flowering_Date, total_seed_mass_g, subset_seed_count, seed_subset_mass, per_seed_weight_g, `100_seed_weight`))

# join seed weights to FT for 2021-2022
#FT_2022$number_of_plants
Seed_weights_2021_2022$replicate <- as.numeric(gsub("rep (\\d)", "\\1", Seed_weights_2021_2022$replicate))
Seed_weights_2021_2022$Flowering_Date <- as.numeric(Seed_weights_2021_2022$Flowering_Date)

### Creating 2022 Dataframe, Adding in the Exp_Year, and removing "2021BED" and "2021ROW" columns
PHENO2022 <- full_join(FT_2022, Seed_weights_2021_2022, by=c('Genotypes', 'number_of_plants'='germinated', 'Condition', 'replicate', '2021BED', '2021ROW', 'Flowering_Date')) %>% mutate(Exp_year = 2022) %>% select(c("Genotypes", "number_of_plants","Condition","replicate","Flowering_Date","Generation", "total_seed_mass_g", "100_seed_weight","Exp_year"))

### Correcting any typos in the genotype
PHENO2022$Genotypes <- gsub("-", "_", PHENO2022$Genotypes)


Seed_weights_2022_2023$PLOT_ID <- as.numeric(Seed_weights_2022_2023$PLOT_ID)

FT_2023$replicate <- as.numeric(gsub("rep (\\d)", "\\1", FT_2023$replicate))

# Creating Dataframe for the 2022-2023 Year
PHENO2023 <- full_join(FT_2023, Seed_weights_2022_2023, by='PLOT_ID') %>% select(-c('Bed_2022', 'Row_2022', 'ROW')) 
PHENO2023 <- filter(PHENO2023, PLOT_ID <= 1036) %>% mutate(Exp_year = 2023)

# standardize colnames
colnames(PHENO2023) <- c("Genotypes", "Condition","replicate","PLOT_ID","number_of_plants","FT_DAYS", "Generation","total_seed_mass_g", "100_seed_weight",  "Exp_year")
PHENO2023$`100_seed_weight` <- as.numeric(PHENO2023$`100_seed_weight`)

#PHENO2023[which(PHENO2023$total_seed_mass_g <0),]
#range(PHENO2023$total_seed_mass_g, na.rm=T)
#PHENO2023[which(PHENO2023$`100_seed_weight` <0),]
#range(PHENO2023$`100_seed_weight`, na.rm=T)

### Subtract the Average weight of a brown bag and the average weight of an envelope to get the true weights
PHENO2023$total_seed_mass_g <- PHENO2023$total_seed_mass_g - 11.24
PHENO2023$`100_seed_weight` <- PHENO2023$`100_seed_weight` - 1.61


which(PHENO2023$number_of_plants < 0)
which(PHENO2022$total_seed_mass_g < 0)

PHENO_FULL <- full_join(PHENO2023, PHENO2022, by=c('Genotypes', 'number_of_plants', 'Condition', 'replicate', 'FT_DAYS'='Flowering_Date', 'Generation', 'total_seed_mass_g', '100_seed_weight', 'Exp_year')) %>% select(- 'PLOT_ID', 'Generation')
#range(PHENO2023$total_seed_mass_g, na.rm=T)
#range(PHENO2023$`100_seed_weight`, na.rm=T)

# remove a few empty rows
PHENO_FULL <- PHENO_FULL %>% filter(!is.na(Genotypes))
PHENO_FULL <- PHENO_FULL %>% filter(replicate<3)

#
PHENO_FT <- PHENO_FULL %>% select(c( "Genotypes","Condition","replicate","number_of_plants","FT_DAYS","Generation", "Exp_year"))
write_delim(PHENO_FT, "FT_per_year.tsv", "\t")

#PHENO_FULL <- PHENO_FULL %>% group_by(Genotypes, Condition, Exp_year) %>% summarise(across(where(is.numeric), \(x) mean(x, na.rm=T))) %>% select(-c(replicate, FT_DAYS))

# seed count based on seed weight and seed weight per 100 seeds
PHENO_FULL$TOTAL_SEED_COUNT <- round(PHENO_FULL$total_seed_mass_g * (100 / PHENO_FULL$`100_seed_weight`))
# one plot with germination of 0 has a seed weight, so replace the 0 with 1
PHENO_FULL[which(PHENO_FULL$Genotypes=="2_156" & PHENO_FULL$number_of_plants==0), 4] <- 1

# seed produced per individual
PHENO_FULL$FECUNDITY <- PHENO_FULL$TOTAL_SEED_COUNT/ PHENO_FULL$number_of_plants
PHENO_FULL$SURVIVAL <- PHENO_FULL$number_of_plants / 10
PHENO_FULL$ABS_FITNESS <- PHENO_FULL$SURVIVAL * PHENO_FULL$FECUNDITY

PHENO_FULL$REL_FITNESS <- PHENO_FULL$ABS_FITNESS / max(PHENO_FULL$ABS_FITNESS, na.rm=T)


#PHENO_FULL$FIT_SEED_PER_PLANT <- PHENO_FULL$AVG_SEED_PER_PLANT/ mean(PHENO_FULL$AVG_SEED_PER_PLANT, na.rm=T)
#PHENO_FULL$FIT_TOTAL_SEED_COUNT <- PHENO_FULL$TOTAL_SEED_COUNT/ mean(PHENO_FULL$TOTAL_SEED_COUNT, na.rm=T)
#quantile(PHENO_FULL$FIT_SEED_PER_PLANT, na.rm=T)
#quantile(PHENO_FULL$FIT_TOTAL_SEED_COUNT, na.rm=T)
#PHENO_FULL$POP_FIT <- PHENO_FULL$FITNESS

write_delim(PHENO_FULL, "FT_FITNESS.tsv", "\t")

# Writing a Function for the summarise command to ignore NA in the mean calculation
mean_for_summarise <- function(x){
  mean(x, na.rm = TRUE)
}

# Dataframe that contains the avveraged replicates for both years
PHENO_FULL_AVERAGE <- PHENO_FULL %>% group_by(Genotypes, Condition, Exp_year) %>% summarise_at(vars("number_of_plants", "FT_DAYS", "total_seed_mass_g", "100_seed_weight", "TOTAL_SEED_COUNT", "FECUNDITY", "SURVIVAL", "ABS_FITNESS", "REL_FITNESS"), mean_for_summarise)
PHENO_FULL_AVERAGE$Generation <- gsub("^1_.*", 18, PHENO_FULL_AVERAGE$Genotypes)
PHENO_FULL_AVERAGE$Generation <- gsub("^2_.*", 28, PHENO_FULL_AVERAGE$Generation)
PHENO_FULL_AVERAGE$Generation <- gsub("^3_.*", 50, PHENO_FULL_AVERAGE$Generation)
PHENO_FULL_AVERAGEGeneration <- gsub("^7_.*", 58, PHENO_FULL_AVERAGE$Generation)
PHENO_FULL_AVERAGE$Generation <- gsub("^*.*_.*", 0, PHENO_FULL_AVERAGE$Generation)


