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
library(googlesheets4)
library(tidyr)

setwd("/bigdata/koeniglab/jmarz001/Ag-Competition/data")
# Load Data

Seed_weights_2022_2023 <- read_delim("SEED_WEIGHTS_2022_2023.csv") %>% select(-c(Date, Notes))
FT_2022_2023 <- read_delim("FT_DAYS_2022_2023.csv")
FT_2023 <- read_delim('FT_2023.tsv')
FT_2022 <- read_delim('FT_2021_2022.tsv')
Genotype_List_2022_2023 <- read_delim("Genotype_List_2023_2023.csv")
Haplo_raw <- read_delim("Competition_Lines_Haplotypes.csv")
Seed_weights_2021_2022 <- read_delim("SEED_WEIGHTS_2021_2022.csv") %>% select(c(Genotypes, germinated, Condition, replicate,  `2021BED`, `2021ROW`, Flowering_Date, total_seed_mass_g, subset_seed_count, seed_subset_mass, per_seed_weight_g, `100_seed_weight`))
# Don't need to subtract weight of the brown bag or seed envelope because raw data has already taken these into account in the total weight and 100 SW calculations 2022
# ^This applies to the data coming from the google sheet, have to subtract bag weights from updated

# join seed weights to FT for 2021-2022
#FT_2022$number_of_plants
Seed_weights_2021_2022$replicate <- as.numeric(gsub("rep (\\d)", "\\1", Seed_weights_2021_2022$replicate))
Seed_weights_2021_2022$Flowering_Date <- as.numeric(Seed_weights_2021_2022$Flowering_Date)

# Changing the Flowering Dates for Geno 63_4 for both single replicates. They were marked down as flowering after 108 and 112 days but since they were a hooded variety the confidence in exact flowering was low. For this reason these values will be changed to NA

Seed_weights_2021_2022$Flowering_Date <- replace(Seed_weights_2021_2022$Flowering_Date, Seed_weights_2021_2022$Genotypes == "63_4" & Seed_weights_2021_2022$Condition == "single", NA)

# Updating Seed subset mass for 2022 with reweighed values (residual outliers)

Seed_tmp <- Seed_weights_2021_2022 %>% filter(Condition == "single")
Seed_tmp <- Seed_tmp %>% filter(Genotypes == '7_87' | Genotypes == '56_5' | Genotypes == '72_7_1' | Genotypes == "7_95" | Genotypes == "7_219" | Genotypes == "7_207" | Genotypes == "7_198" | Genotypes == "2_156" | Genotypes == "60_4" | Genotypes == "3_72" | Genotypes == "7_177" | Genotypes == "7_158" | Genotypes == "59_5_1" | Genotypes == "7_84" | Genotypes == "1_8")
Seed_tmp$seed_subset_mass <- replace(Seed_tmp$seed_subset_mass, Seed_tmp$Genotypes == '7_158' & Seed_tmp$replicate == 1, 7.8)
Seed_tmp$seed_subset_mass <- replace(Seed_tmp$seed_subset_mass, Seed_tmp$Genotypes == '7_158' & Seed_tmp$replicate == 2, 13.9)
Seed_tmp$seed_subset_mass <- replace(Seed_tmp$seed_subset_mass, Seed_tmp$Genotypes == '59_5_1' & Seed_tmp$replicate == 1, 10.5)
Seed_tmp$seed_subset_mass <- replace(Seed_tmp$seed_subset_mass, Seed_tmp$Genotypes == '59_5_1' & Seed_tmp$replicate == 2, 7.8)
Seed_tmp$seed_subset_mass <- replace(Seed_tmp$seed_subset_mass, Seed_tmp$Genotypes == '7_84' & Seed_tmp$replicate == 1, 12.6)
Seed_tmp$seed_subset_mass <- replace(Seed_tmp$seed_subset_mass, Seed_tmp$Genotypes == '7_84' & Seed_tmp$replicate == 2, 10)
Seed_tmp$seed_subset_mass <- replace(Seed_tmp$seed_subset_mass, Seed_tmp$Genotypes == '1_8' & Seed_tmp$replicate == 1, 10.6)
Seed_tmp$seed_subset_mass <- replace(Seed_tmp$seed_subset_mass, Seed_tmp$Genotypes == '1_8' & Seed_tmp$replicate == 2, 11)
Seed_tmp$seed_subset_mass <- replace(Seed_tmp$seed_subset_mass, Seed_tmp$Genotypes == '7_177' & Seed_tmp$replicate == 1, 19.6)
Seed_tmp$seed_subset_mass <- replace(Seed_tmp$seed_subset_mass, Seed_tmp$Genotypes == '7_177' & Seed_tmp$replicate == 2, 8.7)
Seed_tmp$seed_subset_mass <- replace(Seed_tmp$seed_subset_mass, Seed_tmp$Genotypes == '60_4' & Seed_tmp$replicate == 1, 10.2)
Seed_tmp$seed_subset_mass <- replace(Seed_tmp$seed_subset_mass, Seed_tmp$Genotypes == '60_4' & Seed_tmp$replicate == 2, 7.1)
Seed_tmp$seed_subset_mass <- replace(Seed_tmp$seed_subset_mass, Seed_tmp$Genotypes == '7_207' & Seed_tmp$replicate == 1, 11.7)
Seed_tmp$seed_subset_mass <- replace(Seed_tmp$seed_subset_mass, Seed_tmp$Genotypes == '7_207' & Seed_tmp$replicate == 2, 9.3)
Seed_tmp$seed_subset_mass <- replace(Seed_tmp$seed_subset_mass, Seed_tmp$Genotypes == '7_219' & Seed_tmp$replicate == 1, 12.8)
Seed_tmp$seed_subset_mass <- replace(Seed_tmp$seed_subset_mass, Seed_tmp$Genotypes == '7_219' & Seed_tmp$replicate == 2, 8.8)
Seed_tmp$seed_subset_mass <- replace(Seed_tmp$seed_subset_mass, Seed_tmp$Genotypes == '3_72' & Seed_tmp$replicate == 1, 12.2)
Seed_tmp$seed_subset_mass <- replace(Seed_tmp$seed_subset_mass, Seed_tmp$Genotypes == '3_72' & Seed_tmp$replicate == 2, 8.5)
Seed_tmp$seed_subset_mass <- replace(Seed_tmp$seed_subset_mass, Seed_tmp$Genotypes == '72_7_1' & Seed_tmp$replicate == 1, 6.9)
Seed_tmp$seed_subset_mass <- replace(Seed_tmp$seed_subset_mass, Seed_tmp$Genotypes == '72_7_1' & Seed_tmp$replicate == 2, 8.3)
Seed_tmp$seed_subset_mass <- replace(Seed_tmp$seed_subset_mass, Seed_tmp$Genotypes == '7_87' & Seed_tmp$replicate == 1, 7.4)
Seed_tmp$seed_subset_mass <- replace(Seed_tmp$seed_subset_mass, Seed_tmp$Genotypes == '56_5' & Seed_tmp$replicate == 1, 7.5)
Seed_tmp$seed_subset_mass <- replace(Seed_tmp$seed_subset_mass, Seed_tmp$Genotypes == '56_5' & Seed_tmp$replicate == 2, 8.7)
Seed_tmp$seed_subset_mass <- replace(Seed_tmp$seed_subset_mass, Seed_tmp$Genotypes == '7_95' & Seed_tmp$replicate == 1, 12.4)
Seed_tmp$seed_subset_mass <- replace(Seed_tmp$seed_subset_mass, Seed_tmp$Genotypes == '7_95' & Seed_tmp$replicate == 2, 9.8)
Seed_tmp$seed_subset_mass <- replace(Seed_tmp$seed_subset_mass, Seed_tmp$Genotypes == '7_198' & Seed_tmp$replicate == 1, 13.6)
Seed_tmp$seed_subset_mass <- replace(Seed_tmp$seed_subset_mass, Seed_tmp$Genotypes == '2_156' & Seed_tmp$replicate == 1, 11.3)

# Subtracting envelope weight from these values to get true weight

Seed_tmp$seed_subset_mass <- Seed_tmp$seed_subset_mass - 1.61

# Updating the subset seed count for 2_156 rep 1. Previously, 2058 but after recount, 162

Seed_tmp$subset_seed_count <- replace(Seed_tmp$subset_seed_count, Seed_tmp$Genotypes == '2_156' & Seed_tmp$replicate == 1, 162)

# Updating 100 SW by dividing subset mass by subset seed number, then multiplying by 100

Seed_tmp <- Seed_tmp %>% mutate(`100_seed_weight` = (seed_subset_mass/subset_seed_count) * 100)

# Integrating updated 2022 single 100 SW into 2022 data

single_tmp <- subset(Seed_weights_2021_2022, Condition == "single")
vc <- c('7_87', '7_158', '72_7_1', '2_156', '7_219', '7_95', '7_198', '7_207', '60_4', '59_5_1', '3_72', '1_8', '56_5', '7_177', '7_84')

`%!in%` <- Negate(`%in%`)

single_tmp <- single_tmp %>% filter(Genotypes %!in% vc)
single_tmp <- rbind(single_tmp, Seed_tmp)

# Don't need to update TW for Mixed 2022 bc reweighed values are close to previous values
# Updating mixed 2022 100 SW values (residual errors, so reweighed both replicates)

mix_tmp <- subset(Seed_weights_2021_2022, Condition == 'mixed')
mix_tmp$`100_seed_weight` <- replace(mix_tmp$`100_seed_weight`, mix_tmp$Genotypes == '3_146' & mix_tmp$replicate == 1, 6.7- 1.61)
mix_tmp$`100_seed_weight` <- replace(mix_tmp$`100_seed_weight`, mix_tmp$Genotypes == '3_146' & mix_tmp$replicate == 2, 6.7- 1.61)
mix_tmp$`100_seed_weight` <- replace(mix_tmp$`100_seed_weight`, mix_tmp$Genotypes == '7_96' & mix_tmp$replicate == 1, 6.6- 1.61)
mix_tmp$`100_seed_weight` <- replace(mix_tmp$`100_seed_weight`, mix_tmp$Genotypes == '7_96' & mix_tmp$replicate == 2, 7.6- 1.61)
mix_tmp$`100_seed_weight` <- replace(mix_tmp$`100_seed_weight`, mix_tmp$Genotypes == '1_17_2' & mix_tmp$replicate == 1, 4.8- 1.61)
mix_tmp$`100_seed_weight` <- replace(mix_tmp$`100_seed_weight`, mix_tmp$Genotypes == '1_17_2' & mix_tmp$replicate == 2, 7.1- 1.61)
mix_tmp$`100_seed_weight` <- replace(mix_tmp$`100_seed_weight`, mix_tmp$Genotypes == '62_5' & mix_tmp$replicate == 1, 6.1- 1.61)
mix_tmp$`100_seed_weight` <- replace(mix_tmp$`100_seed_weight`, mix_tmp$Genotypes == '62_5' & mix_tmp$replicate == 2, 5.2- 1.61)
mix_tmp$`100_seed_weight` <- replace(mix_tmp$`100_seed_weight`, mix_tmp$Genotypes == '60_4' & mix_tmp$replicate == 1, 6.9 - 1.61)
mix_tmp$`100_seed_weight` <- replace(mix_tmp$`100_seed_weight`, mix_tmp$Genotypes == '60_4' & mix_tmp$replicate == 2, 5.5 - 1.61)
mix_tmp$`100_seed_weight` <- replace(mix_tmp$`100_seed_weight`, mix_tmp$Genotypes == '2_219' & mix_tmp$replicate == 1, 6.5 - 1.61)
mix_tmp$`100_seed_weight` <- replace(mix_tmp$`100_seed_weight`, mix_tmp$Genotypes == '2_219' & mix_tmp$replicate == 2, 6 - 1.61)
mix_tmp$`100_seed_weight` <- replace(mix_tmp$`100_seed_weight`, mix_tmp$Genotypes == '1_187' & mix_tmp$replicate == 1, 6.7 - 1.61)
mix_tmp$`100_seed_weight` <- replace(mix_tmp$`100_seed_weight`, mix_tmp$Genotypes == '1_187' & mix_tmp$replicate == 2, 4.9 - 1.61)
mix_tmp$`100_seed_weight` <- replace(mix_tmp$`100_seed_weight`, mix_tmp$Genotypes == '1_179_1' & mix_tmp$replicate == 1, 7.1 - 1.61)
mix_tmp$`100_seed_weight` <- replace(mix_tmp$`100_seed_weight`, mix_tmp$Genotypes == '1_179_1' & mix_tmp$replicate == 2, 6.5 - 1.61)
mix_tmp$`100_seed_weight` <- replace(mix_tmp$`100_seed_weight`, mix_tmp$Genotypes == '45_3_2' & mix_tmp$replicate == 1, 6.6 - 1.61)
mix_tmp$`100_seed_weight` <- replace(mix_tmp$`100_seed_weight`, mix_tmp$Genotypes == '45_3_2' & mix_tmp$replicate == 2, 7.4 - 1.61)
mix_tmp$`100_seed_weight` <- replace(mix_tmp$`100_seed_weight`, mix_tmp$Genotypes == '1_8' & mix_tmp$replicate == 1, 6.4 - 1.61)
mix_tmp$`100_seed_weight` <- replace(mix_tmp$`100_seed_weight`, mix_tmp$Genotypes == '1_8' & mix_tmp$replicate == 2, 5.5 - 1.61)
mix_tmp$`100_seed_weight` <- replace(mix_tmp$`100_seed_weight`, mix_tmp$Genotypes == '2_58' & mix_tmp$replicate == 1, 7.2- 1.61)
mix_tmp$`100_seed_weight` <- replace(mix_tmp$`100_seed_weight`, mix_tmp$Genotypes == '2_58' & mix_tmp$replicate == 2, 5.6- 1.61)
mix_tmp$`100_seed_weight` <- replace(mix_tmp$`100_seed_weight`, mix_tmp$Genotypes == '2_62' & mix_tmp$replicate == 1, 6.6- 1.61)
mix_tmp$`100_seed_weight` <- replace(mix_tmp$`100_seed_weight`, mix_tmp$Genotypes == '2_62' & mix_tmp$replicate == 2, 5.8- 1.61)
mix_tmp$`100_seed_weight` <- replace(mix_tmp$`100_seed_weight`, mix_tmp$Genotypes == '2_56' & mix_tmp$replicate == 1, 7.2- 1.61)
mix_tmp$`100_seed_weight` <- replace(mix_tmp$`100_seed_weight`, mix_tmp$Genotypes == '2_56' & mix_tmp$replicate == 2, 7.2 - 1.61)
mix_tmp$`100_seed_weight` <- replace(mix_tmp$`100_seed_weight`, mix_tmp$Genotypes == '7_69' & mix_tmp$`2021BED` == 41, 6.5 - 1.61)
mix_tmp$`100_seed_weight` <- replace(mix_tmp$`100_seed_weight`, mix_tmp$Genotypes == '7_69' & mix_tmp$`2021BED` == 35, 6.5- 1.61)
mix_tmp$`100_seed_weight` <- replace(mix_tmp$`100_seed_weight`, mix_tmp$Genotypes == '7_69' & mix_tmp$`2021BED` == 49, 6.2- 1.61)
mix_tmp$`100_seed_weight` <- replace(mix_tmp$`100_seed_weight`, mix_tmp$Genotypes == '7_69' & mix_tmp$`2021BED` == 55, 5.1- 1.61)

# Updating 2022 dataframe by joining subsetted data

Seed_weights_2021_2022 <- rbind(single_tmp, mix_tmp)

# Updating 100 SW because of outliers (mixed)

Seed_weights_2021_2022$`100_seed_weight` <- replace(Seed_weights_2021_2022$`100_seed_weight`, Seed_weights_2021_2022$Genotypes == '2_35' & Seed_weights_2021_2022$Condition == "mixed" & Seed_weights_2021_2022$replicate == 1, 5.1 - 1.61)
Seed_weights_2021_2022$`100_seed_weight` <- replace(Seed_weights_2021_2022$`100_seed_weight`, Seed_weights_2021_2022$Genotypes == '1_75' & Seed_weights_2021_2022$Condition == "mixed" & Seed_weights_2021_2022$replicate == 2, 5.2- 1.61)
Seed_weights_2021_2022$`100_seed_weight` <- replace(Seed_weights_2021_2022$`100_seed_weight`, Seed_weights_2021_2022$Genotypes == '7_165' & Seed_weights_2021_2022$Condition == "mixed" & Seed_weights_2021_2022$replicate == 1, 5.9- 1.61)
Seed_weights_2021_2022$`100_seed_weight` <- replace(Seed_weights_2021_2022$`100_seed_weight`, Seed_weights_2021_2022$Genotypes == '3_4' & Seed_weights_2021_2022$Condition == "mixed" & Seed_weights_2021_2022$replicate == 1, 5.9- 1.61)
Seed_weights_2021_2022$`100_seed_weight` <- replace(Seed_weights_2021_2022$`100_seed_weight`, Seed_weights_2021_2022$Genotypes == '2_148' & Seed_weights_2021_2022$Condition == "mixed" & Seed_weights_2021_2022$replicate == 1, 5.5- 1.61)
Seed_weights_2021_2022$`100_seed_weight` <- replace(Seed_weights_2021_2022$`100_seed_weight`, Seed_weights_2021_2022$Genotypes == '60_4' & Seed_weights_2021_2022$Condition == "mixed" & Seed_weights_2021_2022$replicate == 1, 6.8- 1.61)
Seed_weights_2021_2022$`100_seed_weight` <- replace(Seed_weights_2021_2022$`100_seed_weight`, Seed_weights_2021_2022$Genotypes == '7_58' & Seed_weights_2021_2022$Condition == "mixed" & Seed_weights_2021_2022$replicate == 1, 6.3- 1.61)
Seed_weights_2021_2022$`100_seed_weight` <- replace(Seed_weights_2021_2022$`100_seed_weight`, Seed_weights_2021_2022$Genotypes == '7_211' & Seed_weights_2021_2022$Condition == "mixed" & Seed_weights_2021_2022$replicate == 2, 4.1- 1.61)
Seed_weights_2021_2022$`100_seed_weight` <- replace(Seed_weights_2021_2022$`100_seed_weight`, Seed_weights_2021_2022$Genotypes == '3_146' & Seed_weights_2021_2022$Condition == "mixed" & Seed_weights_2021_2022$replicate == 2, 5.8- 1.61)

# Fixing Typos and averaging the extra replicates for genotypes 1_6 and 7_69

FT_2022$Genotypes <- gsub("-", "_", FT_2022$Genotypes)
FT_2022 <- FT_2022 %>% select(!c("2021BED", "2021ROW"))
tmp <- FT_2022 %>% filter(Genotypes == "1_6" | Genotypes == "7_69") %>% group_by(Genotypes, Condition, replicate, Generation) %>% summarise(across(.cols = c("Flowering_Date", "number_of_plants"), .fns = mean, na.rm = T)) %>% ungroup()
FT_2022 <- FT_2022 %>% filter(Genotypes != "1_6" & Genotypes != "7_69")
FT_2022 <- rbind(tmp, FT_2022)

# Averaging genotypes 1_6 and 7_69 because we have 8 replicates of each, so this script averages the replicates into 2 single reps and 2 mixed reps.

Seed_weights_2021_2022 <- Seed_weights_2021_2022 %>% select(!c("2021BED", "2021ROW"))
Seed_weights_2021_2022$replicate <- gsub(3, 1, Seed_weights_2021_2022$replicate)
Seed_weights_2021_2022$replicate <- gsub(4, 2, Seed_weights_2021_2022$replicate)
Seed_weights_2021_2022$replicate <- as.numeric(Seed_weights_2021_2022$replicate)

tmp <- Seed_weights_2021_2022 %>% filter(Genotypes == "1_6" | Genotypes == "7_69") %>% group_by(Genotypes, Condition, replicate) %>% summarise(across(.cols = where(is.numeric), mean, na.rm= T)) %>% ungroup()
Seed_weights_2021_2022 <- Seed_weights_2021_2022 %>% filter(Genotypes != "1_6" & Genotypes != "7_69")
Seed_weights_2021_2022 <- rbind(tmp, Seed_weights_2021_2022)

# Creating 2022 Dataframe, Adding in the Exp_Year, and removing "2021BED" and "2021ROW" columns
PHENO2022 <- full_join(FT_2022, Seed_weights_2021_2022) %>% mutate(Exp_year = 2022) %>% select(!c("germinated"))

# one plot with germination of 0 has a seed weight, so replace the 0 with 1
PHENO2022[which(PHENO2022$Genotypes=="2_156" & PHENO2022$number_of_plants==0), 6] <- 1

# Calculating some of the phenotypes for 2022

PHENO2022$FECUNDITY <- PHENO2022$total_seed_mass_g / PHENO2022$number_of_plants
PHENO2022$SURVIVAL <- PHENO2022$number_of_plants / 10

### We have 5 genotypes where we accidentally planted 12 seeds instead of 10. For those individuals, it makes sense to adjust their survival rate relative to the 12 seeds planted

# Isolating those individuals and adjusting their survival rates

hmp <- PHENO2022 %>% filter(number_of_plants > 10)
hmp$SURVIVAL <- hmp$number_of_plants / 12

# Adding back into original dataframe

PHENO2022 <- PHENO2022 %>% filter(number_of_plants <= 10)
PHENO2022 <- rbind(PHENO2022, hmp)

PHENO2022$ABS_FITNESS <- PHENO2022$SURVIVAL * PHENO2022$FECUNDITY
PHENO2022$REL_FITNESS <- PHENO2022$ABS_FITNESS / max(PHENO2022$ABS_FITNESS, na.rm=T)

# Adding Centered data for the phenos of 2022

PHENO2022$FECUNDITY <- as.vector(scale(PHENO2022$FECUNDITY, center = TRUE, scale =T))
PHENO2022$ABS_FITNESS <- as.vector(scale(PHENO2022$ABS_FITNESS, center = TRUE, scale = T))

# Removing Columns that aren't necessarily needed and these columns are also not included in the 2023 data set

PHENO2022 <- PHENO2022 %>% select(!c("per_seed_weight_g", "subset_seed_count", "seed_subset_mass"))

# Preparing Raw data to be joined into the 2023 dataframe

Seed_weights_2022_2023$PLOT_ID <- as.numeric(Seed_weights_2022_2023$PLOT_ID)
Seed_weights_2022_2023 <- Seed_weights_2022_2023 %>% filter(PLOT_ID <= 1036)

# Removing one of the duplicated PLOT_ID 839 entries and replacing the existing values with updated TW and 100 SW

Seed_weights_2022_2023 <- Seed_weights_2022_2023 %>% rowid_to_column()
Seed_weights_2022_2023 <- Seed_weights_2022_2023 %>% filter(PLOT_ID != 839 | rowid != 804)
Seed_weights_2022_2023 <- Seed_weights_2022_2023 %>% select(-c(rowid))
Seed_weights_2022_2023$`Brown Bag Weight` <- replace(Seed_weights_2022_2023$`Brown Bag Weight`, Seed_weights_2022_2023$PLOT_ID == 839, 72.4)
Seed_weights_2022_2023$`100 seed weight` <- replace(Seed_weights_2022_2023$`100 seed weight`, Seed_weights_2022_2023$PLOT_ID == 839, 5.3)

# Creating a function to easily update 100 SW values for 2023 season

Hundred_SW_up <- function(plot_id, new_val){
  replace(Seed_weights_2022_2023$`100 seed weight`, Seed_weights_2022_2023$PLOT_ID == plot_id, new_val)
}

# Creating a function to easily update TW for 2023 values

TW_function <- function(id, new_val_tw){
  replace(Seed_weights_2022_2023$`Brown Bag Weight`, Seed_weights_2022_2023$PLOT_ID == id, new_val_tw)
}


# Updating 100 SW for the 2023 season with newly weighed values (Single values)

Seed_weights_2022_2023$`100 seed weight` <- Hundred_SW_up(186, 4.2)
Seed_weights_2022_2023$`100 seed weight` <- Hundred_SW_up(788, 4.7)
Seed_weights_2022_2023$`100 seed weight` <- Hundred_SW_up(813, 8.5)

#Updating 100 SW for the 2023 season with newly weighed values (Mixed values) - Correcting general outliers

Seed_weights_2022_2023$`100 seed weight` <- Hundred_SW_up(142, 5.7)
Seed_weights_2022_2023$`100 seed weight` <- Hundred_SW_up(654, 5.7)
Seed_weights_2022_2023$`100 seed weight` <- Hundred_SW_up(537, 4.6)
Seed_weights_2022_2023$`100 seed weight` <- Hundred_SW_up(838, 5.5)
Seed_weights_2022_2023$`100 seed weight` <- Hundred_SW_up(814, 7.6)
Seed_weights_2022_2023$`100 seed weight` <- Hundred_SW_up(927, 7.6)
Seed_weights_2022_2023$`100 seed weight` <- Hundred_SW_up(531, 6.7)
Seed_weights_2022_2023$`100 seed weight` <- Hundred_SW_up(989, 4.8)
Seed_weights_2022_2023$`100 seed weight` <- Hundred_SW_up(87, 5.9)
Seed_weights_2022_2023$`100 seed weight` <- Hundred_SW_up(215, 5.3)

# Updating 100 SW 2023 values(residual outliers, so reweighed both replicates)

Seed_weights_2022_2023$`100 seed weight` <- Hundred_SW_up(18, 7.1)
Seed_weights_2022_2023$`100 seed weight` <- Hundred_SW_up(23, 6.8)
Seed_weights_2022_2023$`100 seed weight` <- Hundred_SW_up(32, 7.2)
Seed_weights_2022_2023$`100 seed weight` <- Hundred_SW_up(39, 7.2)
Seed_weights_2022_2023$`100 seed weight` <- Hundred_SW_up(40, 6.6)
Seed_weights_2022_2023$`100 seed weight` <- Hundred_SW_up(58, 6.5)
Seed_weights_2022_2023$`100 seed weight` <- Hundred_SW_up(68, 6.2)
Seed_weights_2022_2023$`100 seed weight` <- Hundred_SW_up(82, 7.1)
Seed_weights_2022_2023$`100 seed weight` <- Hundred_SW_up(140, 6.4)
Seed_weights_2022_2023$`100 seed weight` <- Hundred_SW_up(204, 6.3)
Seed_weights_2022_2023$`100 seed weight` <- Hundred_SW_up(227, 6.6)
Seed_weights_2022_2023$`100 seed weight` <- Hundred_SW_up(274, 5.9)
Seed_weights_2022_2023$`100 seed weight` <- Hundred_SW_up(278, 5.9)
Seed_weights_2022_2023$`100 seed weight` <- Hundred_SW_up(324, 6.4)
Seed_weights_2022_2023$`100 seed weight` <- Hundred_SW_up(325, 6.4)
Seed_weights_2022_2023$`100 seed weight` <- Hundred_SW_up(326, 6.4)
Seed_weights_2022_2023$`100 seed weight` <- Hundred_SW_up(334, 6.6)
Seed_weights_2022_2023$`100 seed weight` <- Hundred_SW_up(341, 5.5)
Seed_weights_2022_2023$`100 seed weight` <- Hundred_SW_up(399, 6.6)
Seed_weights_2022_2023$`100 seed weight` <- Hundred_SW_up(424, 6.7)
Seed_weights_2022_2023$`100 seed weight` <- Hundred_SW_up(475, 5.4)
Seed_weights_2022_2023$`100 seed weight` <- Hundred_SW_up(498, 6)
Seed_weights_2022_2023$`100 seed weight` <- Hundred_SW_up(532, 7.8)
Seed_weights_2022_2023$`100 seed weight` <- Hundred_SW_up(546, 5.9)
Seed_weights_2022_2023$`100 seed weight` <- Hundred_SW_up(553, 5.2)
Seed_weights_2022_2023$`100 seed weight` <- Hundred_SW_up(554, 5.7)
Seed_weights_2022_2023$`100 seed weight` <- Hundred_SW_up(572, 6.6)
Seed_weights_2022_2023$`100 seed weight` <- Hundred_SW_up(582, 4.9)
Seed_weights_2022_2023$`100 seed weight` <- Hundred_SW_up(596, 6.6)
Seed_weights_2022_2023$`100 seed weight` <- Hundred_SW_up(654, 5.5)
Seed_weights_2022_2023$`100 seed weight` <- Hundred_SW_up(718, 6.1)
Seed_weights_2022_2023$`100 seed weight` <- Hundred_SW_up(741, 8)
Seed_weights_2022_2023$`100 seed weight` <- Hundred_SW_up(788, 4.7)
Seed_weights_2022_2023$`100 seed weight` <- Hundred_SW_up(792, 5.2)
Seed_weights_2022_2023$`100 seed weight` <- Hundred_SW_up(834, 5.3)
Seed_weights_2022_2023$`100 seed weight` <- Hundred_SW_up(838, 5.1)
Seed_weights_2022_2023$`100 seed weight` <- Hundred_SW_up(839, 5.9)
Seed_weights_2022_2023$`100 seed weight` <- Hundred_SW_up(848, 8)
Seed_weights_2022_2023$`100 seed weight` <- Hundred_SW_up(855, 7.4)
Seed_weights_2022_2023$`100 seed weight` <- Hundred_SW_up(913, 7.5)
Seed_weights_2022_2023$`100 seed weight` <- Hundred_SW_up(938, 5.9)
Seed_weights_2022_2023$`100 seed weight` <- Hundred_SW_up(989, 4.9)
Seed_weights_2022_2023$`100 seed weight` <- Hundred_SW_up(1012, 7.3)

# Updating TW for 2023 (residual outliers, so reweighed both replicates)

Seed_weights_2022_2023$`Brown Bag Weight` <- TW_function(64, 120.6)
Seed_weights_2022_2023$`Brown Bag Weight` <- TW_function(76, 87.5)
Seed_weights_2022_2023$`Brown Bag Weight` <- TW_function(85, 131.7)
Seed_weights_2022_2023$`Brown Bag Weight` <- TW_function(100, 82.1)
Seed_weights_2022_2023$`Brown Bag Weight` <- TW_function(119, 157.9)
Seed_weights_2022_2023$`Brown Bag Weight` <- TW_function(127, 185.8)
Seed_weights_2022_2023$`Brown Bag Weight` <- TW_function(178, 107.4)
Seed_weights_2022_2023$`Brown Bag Weight` <- TW_function(267, 105.7)
Seed_weights_2022_2023$`Brown Bag Weight` <- TW_function(269, 65.8)
Seed_weights_2022_2023$`Brown Bag Weight` <- TW_function(274, 175)
Seed_weights_2022_2023$`Brown Bag Weight` <- TW_function(346, 143.3)
Seed_weights_2022_2023$`Brown Bag Weight` <- TW_function(383, 144.2)
Seed_weights_2022_2023$`Brown Bag Weight` <- TW_function(386, 61.1)
Seed_weights_2022_2023$`Brown Bag Weight` <- TW_function(391, 45.9)
Seed_weights_2022_2023$`Brown Bag Weight` <- TW_function(468, 140.9)
Seed_weights_2022_2023$`Brown Bag Weight` <- TW_function(496, 100.8)
Seed_weights_2022_2023$`Brown Bag Weight` <- TW_function(537, 24)
Seed_weights_2022_2023$`Brown Bag Weight` <- TW_function(578, 237.8)
Seed_weights_2022_2023$`Brown Bag Weight` <- TW_function(590, 194.2)
Seed_weights_2022_2023$`Brown Bag Weight` <- TW_function(599, 226.3)
Seed_weights_2022_2023$`Brown Bag Weight` <- TW_function(614, 219.1)
Seed_weights_2022_2023$`Brown Bag Weight` <- TW_function(633, 224.1)
Seed_weights_2022_2023$`Brown Bag Weight` <- TW_function(641, 222.1)
Seed_weights_2022_2023$`Brown Bag Weight` <- TW_function(692, 214.4)
Seed_weights_2022_2023$`Brown Bag Weight` <- TW_function(781, 216.2)
Seed_weights_2022_2023$`Brown Bag Weight` <- TW_function(783, 203.8)
Seed_weights_2022_2023$`Brown Bag Weight` <- TW_function(788, 52)
Seed_weights_2022_2023$`Brown Bag Weight` <- TW_function(860, 46.6)
Seed_weights_2022_2023$`Brown Bag Weight` <- TW_function(897, 224.7)
Seed_weights_2022_2023$`Brown Bag Weight` <- TW_function(900, 210)
Seed_weights_2022_2023$`Brown Bag Weight` <- TW_function(905, 197.8)
Seed_weights_2022_2023$`Brown Bag Weight` <- TW_function(982, 224.3)
Seed_weights_2022_2023$`Brown Bag Weight` <- TW_function(1010, 238.9)

# Reformatting replicate column
FT_2023$replicate <- as.numeric(gsub("rep (\\d)", "\\1", FT_2023$replicate))


# Creating Dataframe for the 2022-2023 Year
PHENO2023 <- full_join(FT_2023, Seed_weights_2022_2023) %>% select(!c("Bed_2022", "Row_2022", "ROW"))
PHENO2023 <- PHENO2023 %>% mutate(Exp_year = 2023)


# standardize colnames
colnames(PHENO2023) <- c("Genotypes", "Condition","replicate","PLOT_ID","number_of_plants","FT_DAYS", "Generation","total_seed_mass_g", "100_seed_weight",  "Exp_year")
PHENO2023$`100_seed_weight` <- as.numeric(PHENO2023$`100_seed_weight`)
PHENO2023 <- PHENO2023 %>% select(!c("PLOT_ID"))

# Averaging the extra replicates to create 2 single replicates and 2 mixed replicates

tmp <- PHENO2023 %>% filter(Genotypes == "1_6" | Genotypes == "7_69") %>% group_by(Genotypes, Generation, Condition, replicate) %>% summarise(across(.cols = where(is.numeric), .fns = mean, na.rm = T)) %>% ungroup()
PHENO2023 <- PHENO2023 %>% filter(Genotypes != "1_6" & Genotypes != "7_69")
PHENO2023 <- rbind(tmp, PHENO2023)

#PHENO2023[which(PHENO2023$total_seed_mass_g <0),]
#range(PHENO2023$total_seed_mass_g, na.rm=T)
#PHENO2023[which(PHENO2023$`100_seed_weight` <0),]
#range(PHENO2023$`100_seed_weight`, na.rm=T)

### Subtract the Average weight of a brown bag and the average weight of an envelope to get the true weights
PHENO2023$total_seed_mass_g <- PHENO2023$total_seed_mass_g - 11.24
PHENO2023$`100_seed_weight` <- PHENO2023$`100_seed_weight` - 1.61

# Calculating some of the phenos for 2023

PHENO2023$FECUNDITY <- PHENO2023$total_seed_mass_g / PHENO2023$number_of_plants


PHENO2023$SURVIVAL <- PHENO2023$number_of_plants / 10

# We have some genotypes where we planted 12 seeds instead of 10. Adjusting survival rate based on the 12 seeds

NT <- PHENO2023 %>% filter(number_of_plants > 10)
NT$SURVIVAL <- NT$number_of_plants/ 12

PHENO2023 <- PHENO2023 %>% filter(number_of_plants <= 10)
PHENO2023 <- rbind(PHENO2023, NT)


PHENO2023$ABS_FITNESS <- PHENO2023$SURVIVAL * PHENO2023$FECUNDITY
PHENO2023$REL_FITNESS <- PHENO2023$ABS_FITNESS / max(PHENO2023$ABS_FITNESS, na.rm=T)

PHENO2023$FECUNDITY <- as.vector(scale(PHENO2023$FECUNDITY, center = T, scale = T))
PHENO2023$ABS_FITNESS <- as.vector(scale(PHENO2023$ABS_FITNESS, center = T, scale = T))

# Creating the dataframe to contain replicates for both years and all conditions

PHENO_FULL <- full_join(PHENO2023, PHENO2022, by = c("FT_DAYS" = 'Flowering_Date', "Genotypes", "Generation", "Condition", "replicate", "number_of_plants", "total_seed_mass_g", "100_seed_weight", "Exp_year", "FECUNDITY", "SURVIVAL", "ABS_FITNESS", "REL_FITNESS"))
#range(PHENO2023$total_seed_mass_g, na.rm=T)
#range(PHENO2023$`100_seed_weight`, na.rm=T)

# remove a few empty rows
PHENO_FULL <- PHENO_FULL %>% filter(!is.na(Genotypes))
PHENO_FULL <- PHENO_FULL %>% filter(Genotypes != "2_168")

# remove fecundity outliers (centered fecundity > 3 | centered fecundity < -3)

PHENO_FULL$FECUNDITY[which(PHENO_FULL$FECUNDITY >= 3 | PHENO_FULL$FECUNDITY <= -3)] <- NA

# Adding Haplotype data to their corresponding genotypes

Haplo_raw <- Haplo_raw %>% select(c("Haplotype", "Family", "Generation"))
PHENO_FULL <- full_join(PHENO_FULL, Haplo_raw, by = c("Generation", "Genotypes" = "Family")) %>% filter(number_of_plants != "NA")

# Outputting a dataframe to be converted to a google sheet

PHENO_FT <- PHENO_FULL %>% select(c( "Genotypes","Condition","replicate","number_of_plants","FT_DAYS","Generation", "Exp_year"))
write_delim(PHENO_FT, "FT_per_year.tsv", "\t")

#PHENO_FULL <- PHENO_FULL %>% group_by(Genotypes, Condition, Exp_year) %>% summarise(across(where(is.numeric), \(x) mean(x, na.rm=T))) %>% select(-c(replicate, FT_DAYS))

# seed count based on seed weight and seed weight per 100 seeds
PHENO_FULL$TOTAL_SEED_COUNT <- round(PHENO_FULL$total_seed_mass_g * (100 / PHENO_FULL$`100_seed_weight`))

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

# Dataframe that contains the averaged replicates for both years
PHENO_FULL_AVERAGE <- PHENO_FULL %>% group_by(Genotypes, Condition, Exp_year, Haplotype) %>% summarise_at(vars("number_of_plants", "FT_DAYS", "total_seed_mass_g", "100_seed_weight", "TOTAL_SEED_COUNT", "FECUNDITY", "SURVIVAL", "ABS_FITNESS", "REL_FITNESS"), mean_for_summarise)
PHENO_FULL_AVERAGE$Generation <- gsub("^1_.*", 18, PHENO_FULL_AVERAGE$Genotypes)
PHENO_FULL_AVERAGE$Generation <- gsub("^2_.*", 28, PHENO_FULL_AVERAGE$Generation)
PHENO_FULL_AVERAGE$Generation <- gsub("^3_.*", 50, PHENO_FULL_AVERAGE$Generation)
PHENO_FULL_AVERAGE$Generation <- gsub("^7_.*", 58, PHENO_FULL_AVERAGE$Generation)
PHENO_FULL_AVERAGE$Generation <- gsub("^*.*_.*", 0, PHENO_FULL_AVERAGE$Generation)

write_delim(PHENO_FULL_AVERAGE, "FT_FITNESS_AVERAGE.tsv", "\t")
