library(tidyverse)
library(readr)
library(dplyr)
library(ggpubr)
library(ggplot2)
library(googlesheets4)
library(tidyr)

### Load Data

Seed_weights_2022_2023 <- read_sheet('https://docs.google.com/spreadsheets/d/1mxSidDcodD7-Iju9JJZ_YoejhjEqt6JDad3LMGOh61s/edit#gid=1749035245')
FT_2022_2023 <- read_sheet('https://docs.google.com/spreadsheets/d/1Rb1oN4yeqcQDFKfezf3uGLEOlp-s8x1iAIpAO6cqT4E/edit#gid=682364943')
Genotype_List_2022_2023 <- read_sheet('https://docs.google.com/spreadsheets/d/1pKOlthCtyF-T8bbU_96xfoUDcQbbVSnP3jAk6iP3egc/edit#gid=326137907')
Haplo_raw <- read_sheet("https://docs.google.com/spreadsheets/d/13CHW_ZFK7BDMoJ2vgQkm1QhlCm068_m7s1lOCF-lVSc/edit#gid=1521625104")

### Filtering out the Notes column then joining Seed_weights_2022_2023 and FT_2022_2023 to make Conjoined_Data

Seed_weights_2022_2023 <- subset(Seed_weights_2022_2023, select = -Notes)
FT_2022_2023 <- subset(FT_2022_2023, select = -Notes)
Conjoined_Data <- full_join(Seed_weights_2022_2023, FT_2022_2023, by = ("PLOT_ID"))
Conjoined_Data$PLOT_ID <- as.numeric(Conjoined_Data$PLOT_ID)
Conjoined_Data <- subset(Conjoined_Data, Conjoined_Data$PLOT_ID <= 1036)

### Subtract the Average weight of a brown bag and the average weight of an envelope to get the true weights

Conjoined_Data$`Brown Bag Weight` <- Conjoined_Data$`Brown Bag Weight` - 11.24
Conjoined_Data$`100 seed weight` <- Conjoined_Data$`100 seed weight` - 1.61

### Join Conjoined_Data with Genotype_List_2022_2023, rename the Generation column to make it easier to work with, add Fecundity and Fitness columns

Full_Data <- full_join(Conjoined_Data, Genotype_List_2022_2023, by = ("PLOT_ID"))
Full_Data <- select(Full_Data, !c(Date, ROW, `albino count (not included in germination / survival since they don't survive)`, PLANT_ID, Plot_Survival))
Full_Data$Generation <- gsub("^1_.*", 18, Full_Data$Genotypes)
Full_Data$Generation <- gsub("^2_.*", 28, Full_Data$Generation)
Full_Data$Generation <- gsub("^3_.*", 50, Full_Data$Generation)
Full_Data$Generation <- gsub("^7_.*", 58, Full_Data$Generation)
Full_Data$Generation <- gsub("^*.*_.*", 0, Full_Data$Generation)
Full_Data$Fecundity <- Full_Data$`Brown Bag Weight`/(Full_Data$`100 seed weight`/100)
Full_Data$Fitness <- Full_Data$Fecundity * Full_Data$Plot_Germination
Full_Data <- na.omit(Full_Data)
Full_Data <- Full_Data %>% filter(Genotypes != "7_5")

`%!in%` = Negate(`%in%`)
delete_geno <- c("396", "516", "910", "1030", "442", "444", "956", "958")
Full_Data <- Full_Data %>% filter(PLOT_ID %!in% delete_geno)

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

### Creating Single Condition Dataframe

Rep_Single <- Average_Haplo_rep %>% filter(Condition == "single")

### Creating Mixed Condition Dataframe

Rep_Mixed <- Average_Haplo_rep %>% filter(Condition == 'mixed')

### Creating replicate dataframe for correlation graphs

rep1 <- Full_Data %>% filter(replicate == "rep 1")
rep2 <- Full_Data %>% filter(replicate == 'rep 2')
colnames(rep2) <- paste(colnames(rep2), 2, sep = "_")
Replicate_corr_tbl <- full_join(rep1, rep2, by = c("Condition" = "Condition_2", "Generation" = "Generation_2", "Genotypes" = "Genotypes_2"))
Replicate_corr_tbl <- na.omit(Replicate_corr_tbl)

###### Functions to get expected fitness, fecundity, and yield per plant for both conditions
## Numbers col: 1 = mixed, 0 = single

### Fecundity

Exp_Single <- function(x){
  result_single <- x/10
  return(result_single)
}

Exp_Fec_Mixed <- function(x){
  TW_mix <- (x/2) + (Average_Haplo_rep$Atlas_Avg_Fec/2)
  Exp_Fec_mix <- TW_mix/10
  return(Exp_Fec_mix)
}

Average_Haplo_rep$Exp_Fec_Per_Plant <- ifelse(Average_Haplo_rep$Numbers == 1,
       Exp_Fec_Mixed(Average_Haplo_rep$Fecundity),
       Exp_Single(Average_Haplo_rep$Fecundity))

### Fitness

Exp_Fit_Mixed <- function(x){
  fit_mix <- (x/2) + (Average_Haplo_rep$Atlas_Avg_Fitness/2)
  Exp_Fit_mix <- fit_mix/10
  return(Exp_Fit_mix)
}

Average_Haplo_rep$Exp_Fit_Per_Plant <- ifelse(Average_Haplo_rep$Numbers == 1,
                                                    Exp_Fit_Mixed(Average_Haplo_rep$Fitness),
                                                    Exp_Single(Average_Haplo_rep$Fitness))

### Yield

Exp_TW_mix <- function(x){
  TW_mix <- (x/2) + (Average_Haplo_rep$Atlas_Avg_Total_Weight/2)
  Exp_TW <- TW_mix/10
  return(Exp_TW)
}

Average_Haplo_rep$Exp_TW_Per_Plant <- ifelse(Average_Haplo_rep$Numbers == 1,
                                                   Exp_TW_mix(Average_Haplo_rep$`Brown Bag Weight`),
                                                   Exp_Single(Average_Haplo_rep$`Brown Bag Weight`))

### Write Delims
write_delim(Full_Data, "Full_Data")
write_delim(Average_Haplo_rep, "Average_Haplo_rep")
write_delim(Rep_Mixed, "Rep_Mixed")
write_delim(Rep_Single, "Rep_Single")
write_delim(Replicate_corr_tbl, "Replicate_corr_tbl")
