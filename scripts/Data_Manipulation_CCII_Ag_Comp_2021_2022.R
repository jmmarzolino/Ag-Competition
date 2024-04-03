library(tidyverse)
library(readr)
library(dplyr)
library(ggpubr)
library(ggplot2)
library(googlesheets4)
library(tidyr)

### Data from 2021-2022 Folder in Ag_Comp Drive
Seed_weight_2021_2022_raw <- read_sheet('https://docs.google.com/spreadsheets/d/1EDoPrAeOsl0JQ_d7ghsxAR_y238zNHFGeSMEOAfMBFM/edit#gid=592417706')
Haplotype_data_raw <- read_sheet('https://docs.google.com/spreadsheets/d/13CHW_ZFK7BDMoJ2vgQkm1QhlCm068_m7s1lOCF-lVSc/edit#gid=1521625104')

### Cleaning Raw Data for Seed Weights, adding Fitness column
Seed_weight_2021_2022_raw <- Seed_weight_2021_2022_raw %>% select(!c("SUB_LINE_ID", "2021BED", "2021ROW", "FAM_ID","LINE_ID","2021BEDROW", "subset_seed_count", "seed_subset_mass")) %>% filter(total_seed_mass_g != "NA")
Seed_weight_2021_2022_raw$Flowering_Date <- unlist(Seed_weight_2021_2022_raw$Flowering_Date) %>% as.numeric(Seed_weight_2021_2022_raw$Flowering_Date)
Seed_weight_2021_2022_raw$germinated <- unlist(Seed_weight_2021_2022_raw$germinated) %>% as.numeric(Seed_weight_2021_2022_raw$germinated)
Seed_weight_2021_2022_raw$total_seed_mass_g <- unlist(Seed_weight_2021_2022_raw$total_seed_mass_g) %>% as.numeric(Seed_weight_2021_2022_raw$total_seed_mass_g)
Seed_weight_2021_2022_raw$per_seed_weight_g <- unlist(Seed_weight_2021_2022_raw$per_seed_weight_g) %>% as.numeric(Seed_weight_2021_2022_raw$per_seed_weight_g)
Seed_weight_2021_2022_raw$`100_seed_weight` <- unlist(Seed_weight_2021_2022_raw$`100_seed_weight`) %>% as.numeric(Seed_weight_2021_2022_raw$`100_seed_weight`)
Seed_weight_2021_2022_raw$total_seed_estimate <- unlist(Seed_weight_2021_2022_raw$total_seed_estimate) %>% as.numeric(Seed_weight_2021_2022_raw$total_seed_estimate)
colnames(Seed_weight_2021_2022_raw)[which(names(Seed_weight_2021_2022_raw) == "total_seed_estimate")] <- "Fecundity"
Seed_weight_2021_2022_raw$Fitness <- Seed_weight_2021_2022_raw$Fecundity * Seed_weight_2021_2022_raw$germinated

### Creating Generation Column
Seed_weight_2021_2022_raw$Generation <- gsub("^1_.*", 18, Seed_weight_2021_2022_raw$Genotypes)
Seed_weight_2021_2022_raw$Generation <- gsub("^2_.*", 28, Seed_weight_2021_2022_raw$Generation)
Seed_weight_2021_2022_raw$Generation <- gsub("^3_.*", 50, Seed_weight_2021_2022_raw$Generation)
Seed_weight_2021_2022_raw$Generation <- gsub("^7_.*", 58, Seed_weight_2021_2022_raw$Generation)
Seed_weight_2021_2022_raw$Generation <- gsub("^*.*_.*", 0, Seed_weight_2021_2022_raw$Generation)

### Average the replicates
Averaged_Full_2021_2022 <- Seed_weight_2021_2022_raw %>% group_by(Genotypes, Condition, Generation) %>% summarise(across(.cols = c(total_seed_mass_g, '100_seed_weight', Fecundity, Fitness, Flowering_Date), function(x) mean(x))) %>% ungroup()
Averaged_Full_2021_2022$Generation <- as.numeric(Averaged_Full_2021_2022$Generation)

### Cleaning the Haplotype Dataframe and adding it onto Averaged_Full_2021_2022
Haplotype_df <- Haplotype_data_raw %>% select(c("Family", "Haplotype", "Generation"))
Haplotype_df$Haplotype <- unlist(Haplotype_df$Haplotype) %>% as.character(Haplotype_df$Haplotype)
Haplotype_df$Family <- unlist(Haplotype_df$Family) %>% as.character(Haplotype_df$Family)
Averaged_Full_2021_2022 <- full_join(Averaged_Full_2021_2022, Haplotype_df, by = c("Genotypes" = "Family", "Generation"))
Averaged_Full_2021_2022 <- Averaged_Full_2021_2022 %>% filter(total_seed_mass_g != "NA")

### Adding Averaged Atlas values into the table

Averaged_Full_2021_2022 <- Averaged_Full_2021_2022 %>% mutate(Atlas_Avg_Fec = 2363.51,
                                                              Atlas_Avg_Fit = 21347.22,
                                                              Atlas_Avg_TW = 126.8267)

### Creating Single and Mixed Dataframes
Single_2021_2022 <- Averaged_Full_2021_2022 %>% filter(Condition == 'single')
Mixed_2021_2022 <- Averaged_Full_2021_2022 %>% filter(Condition == 'mixed')

### Write Delims
write_delim(Averaged_Full_2021_2022, "Averaged_Full_2021_2022")
write_delim(Single_2021_2022, "Single_2021_2022")
write_delim(Mixed_2021_2022, "Mixed_2021_2022")
