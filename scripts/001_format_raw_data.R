#!/usr/bin/env Rscript

#SBATCH --ntasks=1
#SBATCH --mem=30G
#SBATCH --time=02:00:00
#SBATCH --output=/rhome/jmarz001/bigdata/Ag-Competition/001_format_raw_data.stdout
#SBATCH -p koeniglab

library(tidyverse)
library(ggpubr)

setwd("/bigdata/koeniglab/jmarz001/Ag-Competition/data")
source("/bigdata/koeniglab/jmarz001/Ag-Competition/scripts/CUSTOM_FNS.R")


#### Data from 2021-2022
sw <- read_delim("SEED_WEIGHTS_2021_2022.tsv")
# seed weight 2021-2022 already contains flowering time data

# Averaging genotypes 1_6 and 7_69 because we have 8 Replicates of each, so this script averages the Replicates into 2 single reps and 2 mixed reps.

# two genotypes have 4 Replicates instead of 2 (genotypes 7_69 and 1_6)
# average down to 2 Replicates
tmp <- sw %>% filter(Genotype %in% c("1_6", "7_69"))
tmp[which(tmp$Replicate=="rep 3"),4] <- "rep 1"
tmp[which(tmp$Replicate=="rep 4"),4] <- "rep 2"
tmp2 <- tmp %>% group_by(Genotype, Condition, Replicate) %>% summarise('FT'=mean(FT, na.rm=T), 'TOTAL_MASS'=mean(TOTAL_MASS, na.rm=T), 'SEED_WEIGHT_100'=mean(SEED_WEIGHT_100, na.rm=T))

sw2 <- sw %>% filter(!(Genotype %in% c("1_6", "7_69")))
sw2 <- full_join(sw2, tmp2)

# add year col
sw2$Exp_year <- 2022









#### Data from 2022-2023
## combine flowering time, seed weights, line IDs and plant experiment info
sw3 <- read_delim("SEED_WEIGHTS_2022_2023.tsv")
ft <- read_delim("FT_2022_2023.tsv")
geno_list <- read_delim("Genotype_List_2023_2023.tsv")


# join data sets
ft <- full_join(ft, geno_list, by=c('PLOT_ID', 'ROW', 'PLANT_ID'))
# join the seed weights and remove redundant column
sw4 <- full_join(sw3, ft, by=c('PLOT_ID')) %>% select(-'PLANT_ID')







### Subtract the Average weight of a brown bag and the average weight of an envelope to get the true weights
sw$TOTAL_MASS <- sw$TOTAL_MASS - 11.24
sw$SEED_WEIGHT_100 <- sw$SEED_WEIGHT_100 - 1.61




# add year col
sw4$Exp_year <- 2023


sw4 <- full_join(sw2, sw4, by = c('Genotype', 'Plants', 'Condition', 'Replicate', 'FT', 'TOTAL_MASS', 'SEED_WEIGHT_100', 'Exp_year', 'BED_2021'='BED_2022', 'ROW_2021'='ROW_2022')

which(sw4$Genotype != sw4$PLANT_ID)



 %>% select(- 'PLOT_ID', 'Generation')









# one plot with germination of 0 has a seed weight, so replace the 0 with 1
sw4[which(sw$Plants==0 & !is.na(sw$TOTAL_MASS)), ]


 <- 1



# Calculating some of the phenotypes for 2022

sw$FEC <- sw$TOTAL_MASS / sw$Plants
sw$SURVIVAL <- sw$Plants / 10

### We have 5 genotypes where we accidentally planted 12 seeds instead of 10. For those individuals, it makes sense to adjust their survival rate relative to the 12 seeds planted

# Isolating those individuals and adjusting their survival rates

hmp <- sw %>% filter(Plants > 10)
hmp$SURVIVAL <- hmp$Plants / 12

# Adding back into original dataframe

sw <- sw %>% filter(Plants <= 10)
sw <- rbind(sw, hmp)

sw$ABS_FITNESS <- sw$SURVIVAL * sw$FEC
sw$REL_FITNESS <- sw$ABS_FITNESS / max(sw$ABS_FITNESS, na.rm=T)

# Adding Centered data for the phenos of 2022

sw$FEC <- as.vector(scale(sw$FEC, center = TRUE, scale =T))
sw$ABS_FITNESS <- as.vector(scale(sw$ABS_FITNESS, center = TRUE, scale = T))

# Removing one of the duplicated PLOT_ID 839 entries and replacing the existing values with updated TW and 100 SW

Seed_weights_2022_2023 <- Seed_weights_2022_2023 %>% rowid_to_column()
Seed_weights_2022_2023 <- Seed_weights_2022_2023 %>% filter(PLOT_ID != 839 | rowid != 804)
Seed_weights_2022_2023 <- Seed_weights_2022_2023 %>% select(-c(rowid))
Seed_weights_2022_2023$`Brown Bag Weight` <- replace(Seed_weights_2022_2023$`Brown Bag Weight`, Seed_weights_2022_2023$PLOT_ID == 839, 72.4)
Seed_weights_2022_2023$`100 seed weight` <- replace(Seed_weights_2022_2023$`100 seed weight`, Seed_weights_2022_2023$PLOT_ID == 839, 5.3)


### Subtract the Average weight of a brown bag and the average weight of an envelope to get the true weights
sw$TOTAL_MASS <- sw$TOTAL_MASS - 11.24
sw$SEED_WEIGHT_100 <- sw$SEED_WEIGHT_100 - 1.61

# Subtracting envelope weight from these values to get true weight
Seed_tmp$seed_subset_mass <- Seed_tmp$seed_subset_mass - 1.61

# Don't need to subtract weight of the brown bag or seed envelope because raw data has already taken these into account in the total weight and 100 SW calculations 2022
# ^This applies to the data coming from the google sheet, have to subtract bag weights from updated



which(sw$Plants < 0)
which(sw$TOTAL_MASS < 0)
# remove a few empty rows
PHENO_FULL <- PHENO_FULL %>% filter(!is.na(Genotype))
PHENO_FULL <- PHENO_FULL %>% filter(Replicate<3)




hist(sw$FT)

summary(sw$FT)

## identify outlier values
# total mass
summary(sw$TOTAL_MASS)
high <- median(sw$TOTAL_MASS, na.rm=T) + 2*IQR(sw$TOTAL_MASS, na.rm=T)
low <- median(sw$TOTAL_MASS, na.rm=T) - 2*IQR(sw$TOTAL_MASS, na.rm=T)

ggplot(sw, aes(TOTAL_MASS)) + geom_histogram() + geom_vline(aes(xintercept=high)) + geom_vline(aes(xintercept=low)) + geom_vline(aes(xintercept=median(sw$TOTAL_MASS, na.rm=T)), linetype="dashed") + theme_bw()


sw <- sw %>% filter(!is.na(TOTAL_MASS))

# estimate FECUNDITY from total and 100 seed mass
sw$FEC
sw$FITNESS <- sw$FEC * (sw$Plants/10)

### Average the Replicates
sw_avg <- sw %>% group_by(Genotype, Condition, Replicate) %>% summarise(across(.cols = c(TOTAL_MASS, SEED_WEIGHT_100, FECUNDITY, FITNESS, FT), function(x) mean(x))) %>% ungroup()

### negative TW values?
FT_FITNESS[which(FT_FITNESS$TOTAL_MASS < 0), (7:ncol(FT_FITNESS))]


outlier_cutoff = quantile(test$FEC,0.75, na.rm = TRUE) + (1.5 * IQR(test$FEC, na.rm = TRUE))
ggplot(FT_FITNESS, aes(x = FECUNDITY)) +
  geom_histogram(bins = 50) +
  geom_vline(xintercept = outlier_cutoff, color = 'red')











range(sw$TOTAL_MASS, na.rm=T)
range(sw$SEED_WEIGHT_100, na.rm=T)

# remove a few empty rows
PHENO_FULL <- PHENO_FULL %>% filter(!is.na(Genotype))




# seed count based on seed weight and seed weight per 100 seeds
PHENO_FULL$TOTAL_SEED_COUNT <- round(PHENO_FULL$TOTAL_MASS * (100 / PHENO_FULL$SEED_WEIGHT_100))
# one plot with germination of 0 has a seed weight, so replace the 0 with 1
PHENO_FULL[which(PHENO_FULL$Genotype=="2_156" & PHENO_FULL$Plants==0), 4] <- 1

# seed produced per individual
PHENO_FULL$FEC <- PHENO_FULL$TOTAL_SEED_COUNT/ PHENO_FULL$Plants
PHENO_FULL$SURVIVAL <- PHENO_FULL$Plants / 10
PHENO_FULL$ABS_FITNESS <- PHENO_FULL$SURVIVAL * PHENO_FULL$FEC

PHENO_FULL$REL_FITNESS <- PHENO_FULL$ABS_FITNESS / max(PHENO_FULL$ABS_FITNESS, na.rm=T)

write_delim(PHENO_FULL, "FT_FITNESS.tsv", "\t")

# seed produced per individual
PHENO_FULL$FEC <- PHENO_FULL$TOTAL_SEED_COUNT/ PHENO_FULL$Plants
PHENO_FULL$SURVIVAL <- PHENO_FULL$Plants / 10
PHENO_FULL$ABS_FITNESS <- PHENO_FULL$SURVIVAL * PHENO_FULL$FEC

PHENO_FULL$REL_FITNESS <- PHENO_FULL$ABS_FITNESS / max(PHENO_FULL$ABS_FITNESS, na.rm=T)

# Dataframe that contains the avveraged Replicates for both years
PHENO_FULL_AVERAGE <- PHENO_FULL %>% group_by(Genotype, Condition, Exp_year) %>% summarise_at(vars("Plants", "FT", "TOTAL_MASS", "SEED_WEIGHT_100", "TOTAL_SEED_COUNT", "FECUNDITY", "SURVIVAL", "ABS_FITNESS", "REL_FITNESS"), mean_for_summarise)


  ##########
  # 1) There are some impossible negative numbers.
  # no negative numbers now
  PHENO_FULL %>% reframe(across(where(is.numeric), \(x) range(x, na.rm=T)))

  # 3) No data for 2_168
  PHENO_FULL[which(PHENO_FULL$Genotype == "2_168"),]
  # I think 2_168 are all albino and w low germination, so they are included in this full data set but it's not a mistake and can be filtered
  PHENO_FULL <- PHENO_FULL %>% filter(Genotype=="2_168")

  # 4) More than 8 lines for 7_5    63_4     1_6    7_69 , not always clear what happened.
  # 1_6, and 7_69 had extra Replicates in 2021-2022 by accident & was continued on purpose the next year. the extra lines are accurate and can be treated as extra Replicates and averaged w the others

  # 7_5 has one extra line from 2023, not sure why ...

  # 63_4 is a merge issue


  # 5) Extremely high values of SEED_WEIGHT_100
  ggplot(PHENO_FULL, aes(SEED_WEIGHT_100)) + geom_histogram()
   PHENO_FULL[which(PHENO_FULL$SEED_WEIGHT_100 > 30),]







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

  Full_Data$FEC <- Full_Data$`Brown Bag Weight`/(Full_Data$`100 seed weight`/100)
  Full_Data$Fitness <- Full_Data$FEC * Full_Data$Plot_Germination
  Full_Data <- na.omit(Full_Data)
  Full_Data <- Full_Data %>% filter(Genotype != "7_5")

  delete_geno <- c("396", "516", "910", "1030", "442", "444", "956", "958")
  Full_Data <- Full_Data %>% filter(!(PLOT_ID %in% delete_geno))

  #### Creates a data frame that isolates Atlas Genotype, then averages the Replicates

  Atlas_tbl <- filter(Full_Data, Full_Data$Genotype == "48_5")
  Atlas_tbl <- na.omit(Atlas_tbl)
  Atlas_tbl <- Atlas_tbl %>% mutate(Atlas_Avg_Fecundity = (sum(Atlas_tbl$FEC)/3),
                                    Atlas_Avg_Fitness = (sum(Atlas_tbl$Fitness)/3),
                                    Atlas_Avg_Total_Weight = (sum(Atlas_tbl$`Brown Bag Weight`)/3))

  ### Importing Haplotype Data and cleaning the dataframe

  Haplo_raw <- Haplo_raw %>% select(c("Pedigree", "Haplotype", "Generation", "Family"))
  Haplo_raw$Family <- unlist(Haplo_raw$Family)
  Haplo_raw$Haplotype <- unlist(Haplo_raw$Haplotype)
  Haplo_raw <- Haplo_raw %>% mutate(Pedigree = paste0("UCRKL00000", Haplo_raw$Family))
  colnames(Haplo_raw)[which(names(Haplo_raw) == "Family")] <- "Genotype"

  ### Averaging Replicates

  Average_Haplo_rep <- Full_Data %>% group_by(Genotype, Condition, Generation) %>% summarise(across(.cols = c(`Brown Bag Weight`, `100 seed weight`, Fecundity, Fitness, FT), function(x) mean(x))) %>% ungroup()
  Average_Haplo_rep$Generation <- as.numeric(Average_Haplo_rep$Generation)
  Average_Haplo_rep <- full_join(Haplo_raw, Average_Haplo_rep, by = c("Genotype", "Generation"))
  Average_Haplo_rep <- Average_Haplo_rep %>% filter(`Brown Bag Weight` != "NA")
  Average_Haplo_rep <- Average_Haplo_rep %>%  mutate(Atlas_Avg_Fec = 2363.51,
                                                     Atlas_Avg_Fitness = 21347.22,
                                                     Atlas_Avg_Total_Weight = 126.8267)
  Average_Haplo_rep$Numbers <- ifelse(Average_Haplo_rep$Condition == "mixed", 1, 0)


  ### Creating Replicate dataframe for correlation graphs

  rep1 <- Full_Data %>% filter(Replicate == "rep 1")
  rep2 <- Full_Data %>% filter(Replicate == 'rep 2')
  colnames(rep2) <- paste(colnames(rep2), 2, sep = "_")
  Replicate_corr_tbl <- full_join(rep1, rep2, by = c("Condition" = "Condition_2", "Generation" = "Generation_2", "Genotype" = "Genotype_2"))
  Replicate_corr_tbl <- na.omit(Replicate_corr_tbl)

  ###### Functions to get expected fitness, fecundity, and yield per plant for both conditions
  ## Numbers col: 1 = mixed, 0 = single

  ### Fecundity

  Average_Haplo_rep$Exp_Fec_Per_Plant <- ifelse(Average_Haplo_rep$Numbers == 1,
         Exp_Fec_Mixed(Average_Haplo_rep$FEC),
         Exp_Single(Average_Haplo_rep$FEC))

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
                                      Centered_FT = FT - mean(FT),
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

  source("/bigdata/koeniglab/jmarz001/Ag-Competition/scripts/CUSTOM_FNS.R")
  setwd("/bigdata/koeniglab/jmarz001/Ag-Competition/data")

  ### Load Data

  ## 2022-2023
  Seed_weights_2022_2023 <- read_delim("SEED_WEIGHTS_2022_2023.tsv")
  FT_2022_2023 <- read_delim("FT_2022_2023.tsv")
  FT_2022 <- read_delim('FT_2021_2022.tsv')
  Genotype_List_2022_2023 <- read_delim("Genotype_List_2023_2023.tsv")
  Haplo_raw <- read_delim("Competition_Lines_Haplotypes.csv")

  ## 2021-2022
  Seed_weights_2021_2022 <- read_delim("SEED_WEIGHTS_2021_2022.tsv")
  Seed_weights_2021_2022 <- Seed_weights_2021_2022 %>% select(c(Genotype, germinated, Condition, Replicate,  `BED_2021`, `ROW_2021`, FT, TOTAL_MASS, subset_seed_count, seed_subset_mass, per_seed_weight_g, SEED_WEIGHT_100))


  # join seed weights to FT for 2021-2022
  #FT_2022$Plants
  Seed_weights_2021_2022$Replicate <- as.numeric(gsub("rep (\\d)", "\\1", Seed_weights_2021_2022$Replicate))
  Seed_weights_2021_2022$FT <- as.numeric(Seed_weights_2021_2022$FT)

  PHENO2022 <- full_join(FT_2022, Seed_weights_2021_2022, by=c('Genotype', 'Plants'='germinated', 'Condition', 'Replicate', 'BED_2021', 'ROW_2021', 'FT'))
  # add year col
  PHENO2022$Exp_year <- 2022
  # filter 2022 phenotype sheet to necessary columns
  PHENO2022 <- PHENO2022 %>% select(c("Genotype", "Plants","Condition","Replicate","FT","Generation", "TOTAL_MASS", "SEED_WEIGHT_100","Exp_year"))
  #         "BED_2021"           "ROW_2021"
  PHENO2022$Genotype <- gsub("-", "_", PHENO2022$Genotype)



  Seed_weights_2022_2023$PLOT_ID <- as.numeric(Seed_weights_2022_2023$PLOT_ID)

  FT_2023$Replicate <- as.numeric(gsub("rep (\\d)", "\\1", FT_2023$Replicate))

  PHENO2023 <- full_join(FT_2023, Seed_weights_2022_2023, by='PLOT_ID') %>% select(-c('BED_2022', 'ROW_2022', 'ROW'))
  PHENO2023 <- filter(PHENO2023, PLOT_ID <= 1036)
  PHENO2023$Exp_year <- 2023
  # standardize colnames
  colnames(PHENO2023) <- c("Genotype", "Condition","Replicate","PLOT_ID","Plants","FT", "Generation","TOTAL_MASS", "SEED_WEIGHT_100",  "Exp_year")
  PHENO2023$SEED_WEIGHT_100 <- as.numeric(PHENO2023$SEED_WEIGHT_100)

  #PHENO2023[which(PHENO2023$TOTAL_MASS <0),]
  #range(PHENO2023$TOTAL_MASS, na.rm=T)
  #PHENO2023[which(PHENO2023$SEED_WEIGHT_100 <0),]
  #range(PHENO2023$SEED_WEIGHT_100, na.rm=T)
  ### Subtract the Average weight of a brown bag and the average weight of an envelope to get the true weights
  PHENO2023$TOTAL_MASS <- PHENO2023$TOTAL_MASS - 11.24
  PHENO2023$SEED_WEIGHT_100 <- PHENO2023$SEED_WEIGHT_100 - 1.61



  PHENO_FULL <- full_join(PHENO2023, PHENO2022, by=c('Genotype', 'Plants', 'Condition', 'Replicate', 'FT'='FT', 'Generation', 'TOTAL_MASS', 'SEED_WEIGHT_100', 'Exp_year')) %>% select(- 'PLOT_ID', 'Generation')
  #range(PHENO2023$TOTAL_MASS, na.rm=T)
  #range(PHENO2023$SEED_WEIGHT_100, na.rm=T)

  # remove a few empty rows
  PHENO_FULL <- PHENO_FULL %>% filter(!is.na(Genotype))
  PHENO_FULL <- PHENO_FULL %>% filter(Replicate<3)




  #
  PHENO_FT <- PHENO_FULL %>% select(c( "Genotype","Condition","Replicate","Plants","FT","Generation", "Exp_year"))
  write_delim(PHENO_FT, "FT_per_year.tsv", "\t")



  #PHENO_FULL <- PHENO_FULL %>% group_by(Genotype, Condition, Exp_year) %>% summarise(across(where(is.numeric), \(x) mean(x, na.rm=T))) %>% select(-c(Replicate, FT))

  # seed count based on seed weight and seed weight per 100 seeds
  PHENO_FULL$TOTAL_SEED_COUNT <- round(PHENO_FULL$TOTAL_MASS * (100 / PHENO_FULL$SEED_WEIGHT_100))
  # one plot with germination of 0 has a seed weight, so replace the 0 with 1
  PHENO_FULL[which(PHENO_FULL$Genotype=="2_156" & PHENO_FULL$Plants==0), 4] <- 1

  # seed produced per individual
  PHENO_FULL$FEC <- PHENO_FULL$TOTAL_SEED_COUNT/ PHENO_FULL$Plants
  PHENO_FULL$SURVIVAL <- PHENO_FULL$Plants / 10
  PHENO_FULL$ABS_FITNESS <- PHENO_FULL$SURVIVAL * PHENO_FULL$FEC

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
  grep("-", PHENO_FULL$Genotype)

  # 3) No data for 2_168
  PHENO_FULL[which(PHENO_FULL$Genotype == "2_168"),]
  # I think 2_168 are all albino and w low germination, so they are included in this full data set but it's not a mistake and can be filtered
  PHENO_FULL <- PHENO_FULL %>% filter(Genotype=="2_168")

  # 4) More than 8 lines for 7_5    63_4     1_6    7_69 , not always clear what happened.
  # 1_6, and 7_69 had extra Replicates in 2021-2022 by accident & was continued on purpose the next year. the extra lines are accurate and can be treated as extra Replicates and averaged w the others

  # 7_5 has one extra line from 2023, not sure why ...

  # 63_4 is a merge issue


  # 5) Extremely high values of SEED_WEIGHT_100
  ggplot(PHENO_FULL, aes(SEED_WEIGHT_100)) + geom_histogram()
   PHENO_FULL[which(PHENO_FULL$SEED_WEIGHT_100 > 30),]







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



  Full_Data$FEC <- Full_Data$`Brown Bag Weight`/(Full_Data$`100 seed weight`/100)
  Full_Data$Fitness <- Full_Data$FEC * Full_Data$Plot_Germination
  Full_Data <- na.omit(Full_Data)
  Full_Data <- Full_Data %>% filter(Genotype != "7_5")

  delete_geno <- c("396", "516", "910", "1030", "442", "444", "956", "958")
  Full_Data <- Full_Data %>% filter(!(PLOT_ID %in% delete_geno))

  #### Creates a data frame that isolates Atlas Genotype, then averages the Replicates

  Atlas_tbl <- filter(Full_Data, Full_Data$Genotype == "48_5")
  Atlas_tbl <- na.omit(Atlas_tbl)
  Atlas_tbl <- Atlas_tbl %>% mutate(Atlas_Avg_Fecundity = (sum(Atlas_tbl$FEC)/3),
                                    Atlas_Avg_Fitness = (sum(Atlas_tbl$Fitness)/3),
                                    Atlas_Avg_Total_Weight = (sum(Atlas_tbl$`Brown Bag Weight`)/3))

  ### Importing Haplotype Data and cleaning the dataframe

  Haplo_raw <- Haplo_raw %>% select(c("Pedigree", "Haplotype", "Generation", "Family"))
  Haplo_raw$Family <- unlist(Haplo_raw$Family)
  Haplo_raw$Haplotype <- unlist(Haplo_raw$Haplotype)
  Haplo_raw <- Haplo_raw %>% mutate(Pedigree = paste0("UCRKL00000", Haplo_raw$Family))
  colnames(Haplo_raw)[which(names(Haplo_raw) == "Family")] <- "Genotype"

  ### Averaging Replicates

  Average_Haplo_rep <- Full_Data %>% group_by(Genotype, Condition, Generation) %>% summarise(across(.cols = c(`Brown Bag Weight`, `100 seed weight`, Fecundity, Fitness, FT), function(x) mean(x))) %>% ungroup()
  Average_Haplo_rep$Generation <- as.numeric(Average_Haplo_rep$Generation)
  Average_Haplo_rep <- full_join(Haplo_raw, Average_Haplo_rep, by = c("Genotype", "Generation"))
  Average_Haplo_rep <- Average_Haplo_rep %>% filter(`Brown Bag Weight` != "NA")
  Average_Haplo_rep <- Average_Haplo_rep %>%  mutate(Atlas_Avg_Fec = 2363.51,
                                                     Atlas_Avg_Fitness = 21347.22,
                                                     Atlas_Avg_Total_Weight = 126.8267)
  Average_Haplo_rep$Numbers <- ifelse(Average_Haplo_rep$Condition == "mixed", 1, 0)


  ### Creating Replicate dataframe for correlation graphs

  rep1 <- Full_Data %>% filter(Replicate == "rep 1")
  rep2 <- Full_Data %>% filter(Replicate == 'rep 2')
  colnames(rep2) <- paste(colnames(rep2), 2, sep = "_")
  Replicate_corr_tbl <- full_join(rep1, rep2, by = c("Condition" = "Condition_2", "Generation" = "Generation_2", "Genotype" = "Genotype_2"))
  Replicate_corr_tbl <- na.omit(Replicate_corr_tbl)

  ###### Functions to get expected fitness, fecundity, and yield per plant for both conditions
  ## Numbers col: 1 = mixed, 0 = single

  ### Fecundity

  Average_Haplo_rep$Exp_Fec_Per_Plant <- ifelse(Average_Haplo_rep$Numbers == 1,
         Exp_Fec_Mixed(Average_Haplo_rep$FEC),
         Exp_Single(Average_Haplo_rep$FEC))

  ### Fitness

  Average_Haplo_rep$Exp_Fit_Per_Plant <- ifelse(Average_Haplo_rep$Numbers == 1,
                                                      Exp_Fit_Mixed(Average_Haplo_rep$Fitness),
                                                      Exp_Single(Average_Haplo_rep$Fitness))

  ### Total Weight


  Average_Haplo_rep$Exp_TW_Per_Plant <- ifelse(Average_Haplo_rep$Numbers == 1,
                                                     Exp_TW_mix(Average_Haplo_rep$`Brown Bag Weight`),
                                                     Exp_Single(Average_Haplo_rep$`Brown Bag Weight`))

  ### Adding a column for centered data
  Average_Haplo_rep <- Average_Haplo_rep %>% mutate(Centered_Fit = Fitness - mean(Fitness),
                                      Centered_FT = FT - mean(FT),
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
