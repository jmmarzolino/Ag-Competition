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

# average extra replicates for genotypes 1_6 and 7_69
tmp <- sw %>% filter(Genotype %in% c("1_6", "7_69"))
sw2 <- sw %>% filter(!(Genotype %in% c("1_6", "7_69")))

tmp[which(tmp$Replicate == "rep 3"), 4] <- "rep 1"
tmp[which(tmp$Replicate == "rep 4"), 4] <- "rep 2"
#tmp2 <- tmp %>%
 #       group_by(Genotype, Condition, Replicate) %>%
  #      summarise(across(where(is.numeric), \(x) mean(x, na.rm = TRUE)))



tmp2 <- tmp %>%
        group_by(Genotype, Condition, Replicate) %>%
        summarise("Plants" = mean(Plants), "TOTAL_MASS" = mean(TOTAL_MASS, na.rm = TRUE), "SEED_WEIGHT_100" = mean(SEED_WEIGHT_100, na.rm = TRUE), "FT" = mean(FT, na.rm=TRUE))

ord <- tmp %>% select(c(Genotype, Condition, Replicate,  BED_2021, ROW_2021)) 
ord2 <- ord[order(ord$Genotype, ord$Condition, ord$Replicate),][c(1,3,5,7,9,11,13,15),]

sw2 <- full_join(sw2, tmp2)
avgd <- full_join(tmp2, ord2, by = c("Genotype", "Condition", "Replicate"))
sw6 <- full_join(sw5, avgd, by = c("PLOT_ID", "TOTAL_MASS", "SEED_WEIGHT_100", "FT", "ROW", "Genotype", "Condition", "Replicate", "BED_2022", "ROW_2022", "Plants"))

# add year col
sw2$Exp_year <- 2022



#### Data from 2022-2023
## combine flowering time, seed weights, line IDs and plant experiment info
sw3 <- read_delim("SEED_WEIGHTS_2022_2023.tsv")
ft <- read_delim("FT_2022_2023.tsv")
geno_list <- read_delim("Genotype_List_2023_2023.tsv")

# join data sets
ft <- full_join(ft, geno_list, by = c("PLOT_ID", "ROW", "PLANT_ID"))
# join the seed weights and remove redundant column
sw4 <- full_join(sw3, ft, by = c("PLOT_ID")) %>% select(-"PLANT_ID")

## Tare weights
# Subtract the avg brown bag weight and
# average coin envelope weight from phenotypes
sw4$TOTAL_MASS <- sw4$TOTAL_MASS - 11.24
sw4$SEED_WEIGHT_100 <- sw4$SEED_WEIGHT_100 - 1.61

# average extra replicates for genotypes 1_6 and 7_69
tmp <- sw4 %>% filter(Genotype %in% c("1_6", "7_69"))
sw5 <- sw4 %>% filter(!(Genotype %in% c("1_6", "7_69")))

tmp2 <- tmp %>%
        group_by(Genotype, Condition, Replicate) %>%
        summarise("Plants" = mean(Plants), "TOTAL_MASS" = mean(TOTAL_MASS, na.rm = TRUE), "SEED_WEIGHT_100" = mean(SEED_WEIGHT_100, na.rm = TRUE), "FT" = mean(FT, na.rm=TRUE))

ord <- tmp %>% select(c(Genotype, PLOT_ID, Condition, Replicate, ROW, BED_2022, ROW_2022)) 
ord2 <- ord[order(ord$Genotype, ord$Condition, ord$Replicate),][c(1,3,5,7,9,11,13,15),]


avgd <- full_join(tmp2, ord2, by = c("Genotype", "Condition", "Replicate"))
sw6 <- full_join(sw5, avgd, by = c("PLOT_ID", "TOTAL_MASS", "SEED_WEIGHT_100", "FT", "ROW", "Genotype", "Condition", "Replicate", "BED_2022", "ROW_2022", "Plants"))

#sw4 <- sw4 %>%
 #       group_by(Genotype, Condition, Replicate) %>%
  #      summarise(across(where(is.numeric), \(x) mean(x, na.rm = TRUE)))

# add year col
sw6$Exp_year <- 2023



### Join years of data
sw5 <- full_join(sw2, sw4, by = c("Genotype", "Plants", "Condition", "Replicate", "FT", "TOTAL_MASS", "SEED_WEIGHT_100", "Exp_year", "BED_2021" = "BED_2022", "ROW_2021" = "ROW_2022")) 



### QC checks
# check for duplicated lines
which(table(sw5$PLOT_ID)>1)
which(table(sw5$Genotype)>16)
sw5 <- sw5 %>% select(-c("PLOT_ID", "ROW"))
colnames(sw5) <- c("Genotype", "Plants", "Condition", "Replicate", "BED", "ROW", "FT", "TOTAL_MASS", "SEED_WEIGHT_100", "Exp_year")

# verify phenotype measurements
# are there any phenotyped lines with plant count of 0?
sw5[which(sw5$Plants==0 & !is.na(sw5$TOTAL_MASS)), ]
sw5[which(sw5$Plants==0 & !is.na(sw5$SEED_WEIGHT_100)), ]
sw5[which(sw5$Plants==0 & !is.na(sw5$FT)), ]



# how many empty rows are there?
# empty genotypes?
sw5[which(is.na(sw5$Genotype)),]
# 4 empty rows in 2022 data
# remove those empty rows
sw5 <- sw5 %>% filter(!is.na(Genotype))

# check out plots w 0 plants
(x <- sw5[which(sw5$Plants == 0),])
table(x$Genotype)
# all Genotype 2_168 plots have 0 plants
# line 2_168 were all albino, though this doesn't explain the observation in 'mixed' plots
# filter all empty plot rows
sw5 <- sw5 %>% filter(Plants != 0)


# empty conditions, replicates, etc
sw5[which(is.na(sw5$Condition)),]
sw5[which(is.na(sw5$Replicate)),]
sw5[which(is.na(sw5$Exp_year)),]
sw5[which(is.na(sw5$BED_2021)),]
sw5[which(is.na(sw5$ROW_2021)),]


# any other missing phenotypes?
sw5[which(is.na(sw5$FT)),] # missing FT is fine
sw5[which(is.na(sw5$TOTAL_MASS)),]
sw5[which(is.na(sw5$SEED_WEIGHT_100)),]
sw5[which(!is.na(sw5$TOTAL_MASS) & is.na(sw5$SEED_WEIGHT_100)),] 
# plots with harvested seed but no 100-seed-weight

sw5[which(is.na(sw5$TOTAL_MASS) & is.na(sw5$SEED_WEIGHT_100)),] 
# rows missing both total mass and 100-seed-weight still have FT records
# keeping all records in basic joined data, filter subsequently as needed


# check range of all columns
sw5 %>% reframe(across(where(is.numeric), \(x) range(x, na.rm=TRUE)))

sw5 %>% group_by(Exp_year, Condition) %>% reframe(across(where(is.numeric), \(x) range(x, na.rm=TRUE)))

  # 4) More than 8 lines for 7_5    63_4     1_6    7_69 , not always clear what happened.


  # 7_5 has one extra line from 2023, not sure why ...

  # 63_4 is a merge issue


  # 5) Extremely high values of SEED_WEIGHT_100
  ggplot(PHENO_FULL, aes(SEED_WEIGHT_100)) + geom_histogram()
   PHENO_FULL[which(PHENO_FULL$SEED_WEIGHT_100 > 30),]



  FT_2023$Replicate <- as.numeric(gsub("rep (\\d)", "\\1", FT_2023$Replicate))


  #PHENO_FULL$FIT_SEED_PER_PLANT <- PHENO_FULL$AVG_SEED_PER_PLANT/ mean(PHENO_FULL$AVG_SEED_PER_PLANT, na.rm=TRUE)
  #PHENO_FULL$FIT_TOTAL_SEED_COUNT <- PHENO_FULL$TOTAL_SEED_COUNT/ mean(PHENO_FULL$TOTAL_SEED_COUNT, na.rm=TRUE)
  #quantile(PHENO_FULL$FIT_SEED_PER_PLANT, na.rm=TRUE)
  #quantile(PHENO_FULL$FIT_TOTAL_SEED_COUNT, na.rm=TRUE)
  #PHENO_FULL$POP_FIT <- PHENO_FULL$FITNESS
  #PHENO2023[which(PHENO2023$TOTAL_MASS <0),]
  #range(PHENO2023$TOTAL_MASS, na.rm=TRUE)
  #PHENO2023[which(PHENO2023$SEED_WEIGHT_100 <0),]
  #range(PHENO2023$SEED_WEIGHT_100, na.rm=TRUE)


hist(sw$FT)

summary(sw$FT)

## identify outlier values
outlier_cutoff = quantile(test$FEC,0.75, na.rm = TRUE) + (1.5 * IQR(test$FEC, na.rm = TRUE))

# total mass
summary(sw$TOTAL_MASS)
high <- median(sw$TOTAL_MASS, na.rm=TRUE) + 2*IQR(sw$TOTAL_MASS, na.rm=TRUE)
low <- median(sw$TOTAL_MASS, na.rm=TRUE) - 2*IQR(sw$TOTAL_MASS, na.rm=TRUE)

ggplot(sw, aes(TOTAL_MASS)) + geom_histogram() + geom_vline(aes(xintercept=high)) + geom_vline(aes(xintercept=low)) + geom_vline(aes(xintercept=median(sw$TOTAL_MASS, na.rm=TRUE)), linetype="dashed") + theme_bw()


ggplot(FT_FITNESS, aes(x = FECUNDITY)) +
  geom_histogram(bins = 50) +
  geom_vline(xintercept = outlier_cutoff, color = "#F31919")








# Calculating derived phenotypes

sw$FEC <- sw$TOTAL_MASS / sw$Plants
sw$SURVIVAL <- sw$Plants / 10

# seed count based on seed weight and seed weight per 100 seeds
PHENO_FULL$TOTAL_SEED_COUNT <- round(PHENO_FULL$TOTAL_MASS * (100 / PHENO_FULL$SEED_WEIGHT_100))

# estimate FECUNDITY from total and 100 seed mass
sw$FEC
sw$FITNESS <- sw$FEC * (sw$Plants/10)

### We have 5 genotypes where we accidentally planted 12 seeds instead of 10. For those individuals, it makes sense to adjust their survival rate relative to the 12 seeds planted

# Isolating those individuals and adjusting their survival rates

hmp <- sw %>% filter(Plants > 10)
hmp$SURVIVAL <- hmp$Plants / 12

# Adding back into original dataframe

sw <- sw %>% filter(Plants <= 10)
sw <- rbind(sw, hmp)

sw$ABS_FITNESS <- sw$SURVIVAL * sw$FEC
sw$REL_FITNESS <- sw$ABS_FITNESS / max(sw$ABS_FITNESS, na.rm=TRUE)

# Adding Centered data for the phenos of 2022

sw$FEC <- as.vector(scale(sw$FEC, center = TRUE, scale =TRUE))
sw$ABS_FITNESS <- as.vector(scale(sw$ABS_FITNESS, center = TRUE, scale = TRUE))








### Average the Replicates
sw_avg <- sw %>% group_by(Genotype, Condition, Replicate) %>% summarise(across(.cols = c(TOTAL_MASS, SEED_WEIGHT_100, FECUNDITY, FITNESS, FT), function(x) mean(x))) %>% ungroup()




# seed produced per individual
PHENO_FULL$FEC <- PHENO_FULL$TOTAL_SEED_COUNT/ PHENO_FULL$Plants
PHENO_FULL$SURVIVAL <- PHENO_FULL$Plants / 10
PHENO_FULL$ABS_FITNESS <- PHENO_FULL$SURVIVAL * PHENO_FULL$FEC

PHENO_FULL$REL_FITNESS <- PHENO_FULL$ABS_FITNESS / max(PHENO_FULL$ABS_FITNESS, na.rm=TRUE)









## calculate measures of fitness

avg_seed_per_plant <- PHENO_FULL %>% group_by(Generation) %>% summarise(gen_avg = mean(FIT_SEED_PER_PLANT, na.rm=TRUE))
f18 <- PHENO_FULL %>% filter(Generation==18)

# FITNESS (single only, by year)
## population relative fitness
### each genotypes" fitness relative to whole pop in field

## generation relative fitness
### genotypes" fitness relative to avg seed of the same generation


## Atalas relative fitness

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

  Average_Haplo_rep <- full_join(Haplo_raw, Average_Haplo_rep, by = c("Genotype", "Generation"))
  Average_Haplo_rep <- Average_Haplo_rep %>% filter(`Brown Bag Weight` != "NA")
  Average_Haplo_rep <- Average_Haplo_rep %>%  mutate(Atlas_Avg_Fec = 2363.51,
                    Atlas_Avg_Fitness = 21347.22,
                    Atlas_Avg_Total_Weight = 126.8267)
  Average_Haplo_rep$Numbers <- ifelse(Average_Haplo_rep$Condition == "mixed", 1, 0)


  ### Creating Replicate dataframe for correlation graphs

  rep1 <- Full_Data %>% filter(Replicate == "rep 1")
  rep2 <- Full_Data %>% filter(Replicate == "rep 2")
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





  # seed count based on seed weight and seed weight per 100 seeds
  PHENO_FULL$TOTAL_SEED_COUNT <- round(PHENO_FULL$TOTAL_MASS * (100 / PHENO_FULL$SEED_WEIGHT_100))


  # seed produced per individual
  PHENO_FULL$FEC <- PHENO_FULL$TOTAL_SEED_COUNT/ PHENO_FULL$Plants
  PHENO_FULL$SURVIVAL <- PHENO_FULL$Plants / 10
  PHENO_FULL$ABS_FITNESS <- PHENO_FULL$SURVIVAL * PHENO_FULL$FEC

  PHENO_FULL$REL_FITNESS <- PHENO_FULL$ABS_FITNESS / max(PHENO_FULL$ABS_FITNESS, na.rm=TRUE)













  ## calculate measures of fitness

  avg_seed_per_plant <- PHENO_FULL %>% group_by(Generation) %>% summarise(gen_avg = mean(FIT_SEED_PER_PLANT, na.rm=TRUE))
  f18 <- PHENO_FULL %>% filter(Generation==18)

  # FITNESS (single only, by year)
  ## population relative fitness
  ### each genotypes" fitness relative to whole pop in field

  ## generation relative fitness
  ### genotypes" fitness relative to avg seed of the same generation

  write_delim(PHENO_FULL, "FT_FITNESS.tsv", "\t")
  write_delim(Full_Data, "Full_Data")
  write_delim(Average_Haplo_rep, "Average_Haplo_rep")
