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

### Adding Averaged Atlas values into the table and adding columns for centered data

Atlas_tbl_2021_2022 <- Averaged_Full_2021_2022 %>% filter(Genotypes == "48_5") %>% 
  mutate(Atlas_Avg_Fec = mean(Fecundity),
         Atlas_Avg_Fit = mean(Fitness),
         Atlas_Avg_TW = mean(total_seed_mass_g))

Averaged_Full_2021_2022 <- Averaged_Full_2021_2022 %>% mutate(Atlas_Avg_Fec = 1413.799,
                                                              Atlas_Avg_Fit = 12880.63,
                                                              Atlas_Avg_TW = 64.55571)

Averaged_Full_2021_2022 <- Averaged_Full_2021_2022 %>% mutate(Centered_Fit = Fitness - mean(Fitness, na.rm = TRUE),
                                                  Centered_FT = Flowering_Date - mean(Flowering_Date, na.rm = TRUE),
                                                  Centered_Fec = Fecundity - mean(Fecundity, na.rm = TRUE),
                                                  Centered_TW = total_seed_mass_g - mean(total_seed_mass_g))

### Calculating Contribution of each seed to phenotype

### Fecundity
Exp_Single <- function(x){
  result_single <- x/10
  return(result_single) 
} 

Exp_Fec_Mixed <- function(x){
  TW_mix <- (x/2) + (Averaged_Full_2021_2022$Atlas_Avg_Fec/2) 
  Exp_Fec_mix <- TW_mix/10
  return(Exp_Fec_mix)
}

Averaged_Full_2021_2022$Exp_Fec_Per_Plant <- ifelse(Averaged_Full_2021_2022$Condition == "mixed",
                                              Exp_Fec_Mixed(Averaged_Full_2021_2022$Fecundity),
                                              Exp_Single(Averaged_Full_2021_2022$Fecundity))

### Fitness

Exp_Fit_Mixed <- function(x){
  fit_mix <- (x/2) + (Averaged_Full_2021_2022$Atlas_Avg_Fit/2)
  Exp_Fit_mix <- fit_mix/10
  return(Exp_Fit_mix)
}

Averaged_Full_2021_2022$Exp_Fit_Per_Plant <- ifelse(Averaged_Full_2021_2022$Condition == 'mixed',
                                              Exp_Fit_Mixed(Averaged_Full_2021_2022$Fitness),
                                              Exp_Single(Averaged_Full_2021_2022$Fitness))

### Total Weight

Exp_TW_mix <- function(x){
  TW_mix <- (x/2) + (Averaged_Full_2021_2022$Atlas_Avg_TW/2)
  Exp_TW <- TW_mix/10
  return(Exp_TW)
}

Averaged_Full_2021_2022$Exp_TW_Per_Plant <- ifelse(Averaged_Full_2021_2022$Condition == 'mixed',
                                             Exp_TW_mix(Averaged_Full_2021_2022$total_seed_mass_g),
                                             Exp_Single(Averaged_Full_2021_2022$total_seed_mass_g))


### Creating Single and Mixed Dataframes
Single_2021_2022 <- Averaged_Full_2021_2022 %>% filter(Condition == 'single')
Mixed_2021_2022 <- Averaged_Full_2021_2022 %>% filter(Condition == 'mixed')

### (Single) Making a dataframe suitable for plotting correlation. 7_69 and 1_6 have 4 replicates, this code averages them down to the regular two replicates

tmp <- Seed_weight_2021_2022_raw %>% filter(Condition == "single")
d <- c("1_6", "7_69")
e <- c("rep 3", "rep 4")
tmp <- tmp %>% filter(Genotypes %in% d) %>% filter(replicate %in% e)
tmp1 <- tmp %>% group_by(Genotypes, Generation) %>% summarise(across(where(is.numeric), mean))
tmp1 <- tmp1 %>% mutate(replicate = "rep 2",
                        Condition = "single")

rmp <- Seed_weight_2021_2022_raw %>% filter(Condition == "single")
d <- c("1_6", "7_69")
e <- c("rep 1", "rep 2")
rmp <- rmp %>% filter(Genotypes %in% d) %>% filter(replicate %in% e)
rmp1 <- rmp %>% group_by(Genotypes, Generation) %>% summarise(across(where(is.numeric), mean))
rmp1 <- rmp1 %>% mutate(replicate = "rep 1",
                        Condition = "single")

`%!in%` = Negate(`%in%`)
ymp <- Seed_weight_2021_2022_raw %>% filter(Condition == "single") %>% filter(Genotypes %!in% d)
Rep_2021_2022_Single <- rbind(rmp1, tmp1, ymp)

### (Single) Reformatting dataframe

R1 <- Rep_2021_2022_Single %>% filter(replicate == "rep 1")
R2 <- Rep_2021_2022_Single %>% filter(replicate == "rep 2")
colnames(R2) <- paste(colnames(R2), sep = "_", 2)
Rep_2021_2022_Single <- full_join(R1, R2, by = c("Condition" = "Condition_2", "Generation" = "Generation_2", "Genotypes" = "Genotypes_2"))
Rep_2021_2022_Single <- na.omit(Rep_2021_2022_Single)

### (Mixed) Making a dataframe suitable for plotting correlation. 7_69 and 1_6 have 4 replicates, this code averages them down to the regular two replicates

tmp <- Seed_weight_2021_2022_raw %>% filter(Condition == "mixed")
d <- c("1_6", "7_69")
e <- c("rep 3", "rep 4")
tmp <- tmp %>% filter(Genotypes %in% d) %>% filter(replicate %in% e)
tmp1 <- tmp %>% group_by(Genotypes, Generation) %>% summarise(across(where(is.numeric), mean))
tmp1 <- tmp1 %>% mutate(replicate = "rep 2",
                        Condition = "mixed")

rmp <- Seed_weight_2021_2022_raw %>% filter(Condition == "mixed")
d <- c("1_6", "7_69")
e <- c("rep 1", "rep 2")
rmp <- rmp %>% filter(Genotypes %in% d) %>% filter(replicate %in% e)
rmp1 <- rmp %>% group_by(Genotypes, Generation) %>% summarise(across(where(is.numeric), mean))
rmp1 <- rmp1 %>% mutate(replicate = "rep 1",
                        Condition = "mixed")

`%!in%` = Negate(`%in%`)
ymp <- Seed_weight_2021_2022_raw %>% filter(Condition == "mixed") %>% filter(Genotypes %!in% d)
Rep_2021_2022_Mixed <- rbind(rmp1, tmp1, ymp)

### (Mixed) Reformatting Dataframe

R1 <- Rep_2021_2022_Mixed %>% filter(replicate == "rep 1")
R2 <- Rep_2021_2022_Mixed %>% filter(replicate == "rep 2")
colnames(R2) <- paste(colnames(R2), sep = "_", 2)
Rep_2021_2022_Mixed <- full_join(R1, R2, by = c("Condition" = "Condition_2", "Generation" = "Generation_2", "Genotypes" = "Genotypes_2"))
Rep_2021_2022_Mixed <- na.omit(Rep_2021_2022_Mixed)

### Write Delims
write_delim(Averaged_Full_2021_2022, "Averaged_Full_2021_2022.tsv")
write_delim(Single_2021_2022, "Single_2021_2022.tsv")
write_delim(Mixed_2021_2022, "Mixed_2021_2022.tsv")
write_delim(Seed_weight_2021_2022_raw, "Seed_weight_2021_2022_raw")
write_delim(Rep_2021_2022_Single, "Rep_2021_2022_Single")
write_delim(Rep_2021_2022_Mixed, "Rep_2021_2022_Mixed")







FT_FITNESS <- read_sheet('https://docs.google.com/spreadsheets/d/15-7DX0YVGhldTwaW6nkKnNhryFmxEwo2ZHiZwAtBu58/edit#gid=1001803440')

### Create function to unlist and convert vectors to numeric
unlist_numeric <- function(x){
  unlist(x) %>% 
    as.numeric(x)
}

### Replace negative TW with NA, unlist all vectors and convert them to numeric, fill in missing generation values
FT_FITNESS <- FT_FITNESS %>% relocate(Exp_year, .after = replicate)
FT_FITNESS[5:14] <- lapply(FT_FITNESS[5:14], unlist_numeric)
FT_FITNESS[which(FT_FITNESS$total_seed_mass_g < 0), (7:ncol(FT_FITNESS))] <- NA

FT_FITNESS$Generation <- gsub("^1_.*", 18, FT_FITNESS$Genotypes)
FT_FITNESS$Generation <- gsub("^2_.*", 28, FT_FITNESS$Generation)
FT_FITNESS$Generation <- gsub("^3_.*", 50, FT_FITNESS$Generation)
FT_FITNESS$Generation <- gsub("^7_.*", 58, FT_FITNESS$Generation)
FT_FITNESS$Generation <- gsub("^*.*_.*", 0, FT_FITNESS$Generation)

outlier_cutoff = quantile(test$FECUNDITY,0.75, na.rm = TRUE) + (1.5 * IQR(test$FECUNDITY, na.rm = TRUE))
ggplot(FT_FITNESS, aes(x = FECUNDITY)) +
  geom_histogram(bins = 50) +
  geom_vline(xintercept = outlier_cutoff, color = 'red')






