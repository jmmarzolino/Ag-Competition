#!/usr/bin/env Rscript

#SBATCH --ntasks=1
#SBATCH --mem=30G
#SBATCH --time=02:00:00
#SBATCH --output=/rhome/jmarz001/bigdata/Ag-Competition/competition1.stdout
#SBATCH -p koeniglab

library(tidyverse)
library(readr)
library(dplyr)
library(ggpubr)
library(ggplot2)

library(tidyr)
source("rhome/jmarz001/bigdata/Ag-Competition/scripts/CUSTOM_FNS.R")
setwd("/bigdata/koeniglab/jmarz001/Ag-Competition/data")
### Load Data

Seed_weights_2022_2023 <- read_delim("SEED_WEIGHTS_2022_2023.csv")
Seed_weights_2022_2023 <- Seed_weights_2022_2023 %>% select(-c(Date, Notes))
FT_2022_2023 <- read_delim("FT_DAYS_2022_2023.csv")
FT_2023 <- read_delim('FT_2023.tsv')
FT_2022 <- read_delim('FT_2021_2022.tsv')

Genotype_List_2022_2023 <- read_delim("Genotype_List_2023_2023.csv")
Haplo_raw <- read_delim("Competition_Lines_Haplotypes.csv")

Seed_weights_2021_2022 <- read_delim("SEED_WEIGHTS_2021_2022.csv")
Seed_weights_2021_2022 <- Seed_weights_2021_2022 %>% select(c(Genotypes, germinated, Condition, replicate,  `2021BED`, `2021ROW`, Flowering_Date, total_seed_mass_g, subset_seed_count, seed_subset_mass, per_seed_weight_g, `100_seed_weight`))


# join seed weights to FT for 2021-2022
#FT_2022$number_of_plants
Seed_weights_2021_2022$replicate <- as.numeric(gsub("rep (\\d)", "\\1", Seed_weights_2021_2022$replicate))
Seed_weights_2021_2022$Flowering_Date <- as.numeric(Seed_weights_2021_2022$Flowering_Date)

PHENO2022 <- full_join(FT_2022, Seed_weights_2021_2022, by=c('Genotypes', 'number_of_plants'='germinated', 'Condition', 'replicate', '2021BED', '2021ROW', 'Flowering_Date'))
# add year col
PHENO2022$Exp_year <- 2022
# filter 2022 phenotype sheet to necessary columns
PHENO2022 <- PHENO2022 %>% select(c("Genotypes", "number_of_plants","Condition","replicate","Flowering_Date","Generation", "total_seed_mass_g", "100_seed_weight","Exp_year"))
#         "2021BED"           "2021ROW"
PHENO2022$Genotypes <- gsub("-", "_", PHENO2022$Genotypes)



Seed_weights_2022_2023$PLOT_ID <- as.numeric(Seed_weights_2022_2023$PLOT_ID)

FT_2023$replicate <- as.numeric(gsub("rep (\\d)", "\\1", FT_2023$replicate))

PHENO2023 <- full_join(FT_2023, Seed_weights_2022_2023, by='PLOT_ID') %>% select(-c('Bed_2022', 'Row_2022', 'ROW'))
PHENO2023 <- filter(PHENO2023, PLOT_ID <= 1036)
PHENO2023$Exp_year <- 2023
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




##########
# 1) There are some impossible negative numbers.
# no negative numbers now
PHENO_FULL %>% reframe(across(where(is.numeric), \(x) range(x, na.rm=T)))

# 2) Typos in the Genotype names for 1_105-1  1_17-2
# no more dashes now
grep("-", PHENO_FULL$Genotypes)

# 3) No data for 2_168
PHENO_FULL[which(PHENO_FULL$Genotypes == "2_168"),]
# I think 2_168 are all albino and w low germination, so they are included in this full data set but it's not a mistake and can be filtered
PHENO_FULL <- PHENO_FULL %>% filter(Genotypes=="2_168")

# 4) More than 8 lines for 7_5    63_4     1_6    7_69 , not always clear what happened.
# 1_6, and 7_69 had extra replicates in 2021-2022 by accident & was continued on purpose the next year. the extra lines are accurate and can be treated as extra replicates and averaged w the others

# 7_5 has one extra line from 2023, not sure why ...

# 63_4 is a merge issue


# 5) Extremely high values of X100_seed_weight
ggplot(PHENO_FULL, aes(`100_seed_weight`)) + geom_histogram()
 PHENO_FULL[which(PHENO_FULL$`100_seed_weight` > 30),]







## calculate measures of fitness

avg_seed_per_plant <- PHENO_FULL %>% group_by(Generation) %>% summarise(gen_avg = mean(FIT_SEED_PER_PLANT, na.rm=T))
f18 <- PHENO_FULL %>% filter(Generation==18)

# FITNESS (single only, by year)
## population relative fitness
### each genotypes' fitness relative to whole pop in field

## generation relative fitness
### genotypes' fitness relative to avg seed of the same generation


## Atalas relative fitness












### Filtering out the Notes column then joining Seed_weights_2022_2023 and FT_2022_2023 to make Conjoined_Data

Seed_weights_2022_2023 <- subset(Seed_weights_2022_2023, select = -Notes)
FT_2022_2023 <- subset(FT_2022_2023, select = -Notes)
Seed_weights_2022_2023$PLOT_ID <- as.numeric(Seed_weights_2022_2023$PLOT_ID)
Conjoined_Data <- full_join(Seed_weights_2022_2023, FT_2022_2023, by = ("PLOT_ID"))
Conjoined_Data <- subset(Conjoined_Data, Conjoined_Data$PLOT_ID <= 1036)
Conjoined_Data$`100 seed weight` <- as.numeric(Conjoined_Data$`100 seed weight`)


### Join Conjoined_Data with Genotype_List_2022_2023, rename the Generation column to make it easier to work with, add Fecundity and Fitness columns

Full_Data <- full_join(Conjoined_Data, Genotype_List_2022_2023, by = ("PLOT_ID"))
Full_Data <- select(Full_Data, !c(Date, ROW, `albino count (not included in germination / survival since they don't survive)`, PLANT_ID, Plot_Survival)) #'
Full_Data <- add_generation(Full_Data)

Full_Data$Fecundity <- Full_Data$`Brown Bag Weight`/(Full_Data$`100 seed weight`/100)
Full_Data$Fitness <- Full_Data$Fecundity * Full_Data$Plot_Germination
Full_Data <- na.omit(Full_Data)
Full_Data <- Full_Data %>% filter(Genotypes != "7_5")

delete_geno <- c("396", "516", "910", "1030", "442", "444", "956", "958")
Full_Data <- Full_Data %>% filter(!(PLOT_ID %in% delete_geno))

#### Creates a data frame that isolates Atlas Genotypes, then averages the replicates

Atlas_tbl <- filter(Full_Data, Full_Data$Genotypes == "48_5")
Atlas_tbl <- na.omit(Atlas_tbl)
Atlas_tbl <- Atlas_tbl %>% mutate(Atlas_Avg_Fecundity = (sum(Atlas_tbl$Fecundity)/3),
                                  Atlas_Avg_Fitness = (sum(Atlas_tbl$Fitness)/3),
                                  Atlas_Avg_Total_Weight = (sum(Atlas_tbl$`Brown Bag Weight`)/3))

### Importing Haplotype Data and cleaning the dataframe

Haplo_raw <- Haplo_raw %>% select(c("Pedigree", "Haplotype", "Generation", "Family"))
Haplo_raw$Family <- unlist(Haplo_raw$Family)
Haplo_raw$Haplotype <- unlist(Haplo_raw$Haplotype)
Haplo_raw <- Haplo_raw %>% mutate(Pedigree = paste0("UCRKL00000", Haplo_raw$Family))
colnames(Haplo_raw)[which(names(Haplo_raw) == "Family")] <- "Genotypes"

### Averaging Replicates

Average_Haplo_rep <- Full_Data %>% group_by(Genotypes, Condition, Generation) %>% summarise(across(.cols = c(`Brown Bag Weight`, `100 seed weight`, Fecundity, Fitness, FT_DAYS), function(x) mean(x))) %>% ungroup()
Average_Haplo_rep$Generation <- as.numeric(Average_Haplo_rep$Generation)
Average_Haplo_rep <- full_join(Haplo_raw, Average_Haplo_rep, by = c("Genotypes", "Generation"))
Average_Haplo_rep <- Average_Haplo_rep %>% filter(`Brown Bag Weight` != "NA")
Average_Haplo_rep <- Average_Haplo_rep %>%  mutate(Atlas_Avg_Fec = 2363.51,
                                                   Atlas_Avg_Fitness = 21347.22,
                                                   Atlas_Avg_Total_Weight = 126.8267)
Average_Haplo_rep$Numbers <- ifelse(Average_Haplo_rep$Condition == "mixed", 1, 0)


### Creating replicate dataframe for correlation graphs

rep1 <- Full_Data %>% filter(replicate == "rep 1")
rep2 <- Full_Data %>% filter(replicate == 'rep 2')
colnames(rep2) <- paste(colnames(rep2), 2, sep = "_")
Replicate_corr_tbl <- full_join(rep1, rep2, by = c("Condition" = "Condition_2", "Generation" = "Generation_2", "Genotypes" = "Genotypes_2"))
Replicate_corr_tbl <- na.omit(Replicate_corr_tbl)

###### Functions to get expected fitness, fecundity, and yield per plant for both conditions
## Numbers col: 1 = mixed, 0 = single

### Fecundity

Average_Haplo_rep$Exp_Fec_Per_Plant <- ifelse(Average_Haplo_rep$Numbers == 1,
       Exp_Fec_Mixed(Average_Haplo_rep$Fecundity),
       Exp_Single(Average_Haplo_rep$Fecundity))

### Fitness

Average_Haplo_rep$Exp_Fit_Per_Plant <- ifelse(Average_Haplo_rep$Numbers == 1,
                                                    Exp_Fit_Mixed(Average_Haplo_rep$Fitness),
                                                    Exp_Single(Average_Haplo_rep$Fitness))

### Total Weight

Exp_TW_mix <- function(x){
  TW_mix <- (x/2) + (Average_Haplo_rep$Atlas_Avg_Total_Weight/2)
  Exp_TW <- TW_mix/10
  return(Exp_TW)
}

Average_Haplo_rep$Exp_TW_Per_Plant <- ifelse(Average_Haplo_rep$Numbers == 1,
                                                   Exp_TW_mix(Average_Haplo_rep$`Brown Bag Weight`),
                                                   Exp_Single(Average_Haplo_rep$`Brown Bag Weight`))

### Adding a column for centered data
Average_Haplo_rep <- Average_Haplo_rep %>% mutate(Centered_Fit = Fitness - mean(Fitness),
                                    Centered_FT = FT_DAYS - mean(FT_DAYS),
                                    Centered_Fec = Fecundity - mean(Fecundity),
                                    Centered_TW = `Brown Bag Weight` - mean(`Brown Bag Weight`, na.rm = TRUE))


### Creating Single Condition Dataframe

Rep_Single <- Average_Haplo_rep %>% filter(Condition == "single")

### Creating Mixed Condition Dataframe

Rep_Mixed <- Average_Haplo_rep %>% filter(Condition == 'mixed')


### Write Delims
write_delim(Full_Data, "Full_Data")
write_delim(Average_Haplo_rep, "Average_Haplo_rep")
write_delim(Rep_Mixed, "Rep_Mixed")
write_delim(Rep_Single, "Rep_Single")
write_delim(Replicate_corr_tbl, "Replicate_corr_tbl")
