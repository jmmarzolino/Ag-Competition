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

### Split Single_Condition into two data frames in order to average the replicates of the single Condition

Single_Condition <- subset(Full_Data, Full_Data$Condition == "single")
Single_rep_1 <- subset(Single_Condition, replicate == "rep 1")
Single_rep_2 <- subset(Single_Condition, replicate == "rep 2")
colnames(Single_rep_2)[1:14] <- paste(colnames(Single_rep_2)[c(1:14)], '_2', sep = '_')
Rep_Single <- inner_join(Single_rep_1, Single_rep_2, by = c("Genotypes" = "Genotypes__2","Condition" = "Condition__2"))
Rep_Single <- mutate(Rep_Single, Avg_Fecundity = (Rep_Single$Fecundity + Rep_Single$Fecundity__2)/2,
                     Avg_Fit = (Rep_Single$Fitness + Rep_Single$Fitness__2)/2,
                     Avg_FT = (Rep_Single$FT_DAYS + Rep_Single$FT_DAYS__2)/2,
                     Avg_Total_Weight = (Rep_Single$`Brown Bag Weight` + Rep_Single$`Brown Bag Weight__2`)/2,
                     Expect_Seed_Mass_per_Plant = (Rep_Single$`Brown Bag Weight` + Rep_Single$`Brown Bag Weight__2`)/10)


#### Creates a data frame that isolates Atlas Genotypes, then averages the replicates

Atlas_tbl <- filter(Full_Data, Full_Data$Genotypes == "48_5")
Atlas_tbl <- na.omit(Atlas_tbl)
Atlas_tbl <- Atlas_tbl %>% mutate(Atlas_Avg_Fecundity = (sum(Atlas_tbl$Fecundity)/3),
                                  Atlas_Avg_Fitness = (sum(Atlas_tbl$Fitness)/3),
                                  Atlas_Avg_Total_Weight = (sum(Atlas_tbl$`Brown Bag Weight`)/3))

### Average the mixed replicate data and add average total weight, fecundity, fitness, and expected fecundity

Mixed_Condition <- subset(Full_Data, Full_Data$Condition == "mixed")
Mixed_rep_1 <- subset(Mixed_Condition, replicate == "rep 1")
Mixed_rep_2 <- subset(Mixed_Condition, replicate == "rep 2")
colnames(Mixed_rep_2)[1:14] <- paste(colnames(Mixed_rep_2)[c(1:14)], '_2', sep = '_')
Rep_Mixed <- inner_join(Mixed_rep_1, Mixed_rep_2, by = c("Genotypes" = "Genotypes__2","Condition" = "Condition__2"))
Rep_Mixed <- Rep_Mixed %>% mutate(Avg_Total_weight = (Rep_Mixed$`Brown Bag Weight` + Rep_Mixed$`Brown Bag Weight__2`)/2,
                                  Avg_Fecundity = Rep_Mixed$Fecundity + Rep_Mixed$Fecundity__2,
                                  Atlas_Avg_Fec = 2363.51,
                                  Atlas_Avg_Fitness = 21347.22)
Rep_Mixed <- Rep_Mixed %>% mutate(Expected_Fecundity = (Avg_Fecundity - Atlas_Avg_Fec)/2 + (Atlas_Avg_Fec)/2)

### Averaging replicates to correlate the replicate data, adding average fitness, fecundity, 100 seed weight, flowering time

rep1 <- subset(Full_Data, replicate == "rep 1")
rep2 <- subset(Full_Data, replicate == "rep 2")
colnames(rep2)[1:13] <- paste(colnames(rep2)[c(1:13)], '_2', sep = '_')
Side_By_Side_Replicates <- inner_join(rep1, rep2, by = c("Genotypes" = "Genotypes__2","Condition" = "Condition__2"))
Side_By_Side_Replicates <- Side_By_Side_Replicates %>% mutate(Avg_Fec = (Side_By_Side_Replicates$Fecundity + Side_By_Side_Replicates$Fecundity__2)/2,
                                                              Avg_Total_Weight = (Side_By_Side_Replicates$`Brown Bag Weight`+Side_By_Side_Replicates$`Brown Bag Weight__2`)/2,
                                                              Avg_Fit = (Side_By_Side_Replicates$Fitness + Side_By_Side_Replicates$Fitness__2)/2,
                                                              Avg_100_SW = (Side_By_Side_Replicates$`100 seed weight` + Side_By_Side_Replicates$`100 seed weight__2`)/2,
                                                              Avg_FT = (Side_By_Side_Replicates$FT_DAYS + Side_By_Side_Replicates$FT_DAYS__2)/2,
                                                              Atlas_Avg_Fec = 2363.51,
                                                              Atlas_Avg_Fitness = 21347.22,
                                                              Atlas_Avg_Total_Weight = 126.8267)
Side_By_Side_Replicates$Numbers <- ifelse(Side_By_Side_Replicates$Condition == "mixed", 1, 0)

### Functions to get expected per plant fecundity for both conditions
### Numbers col: 1 = mixed, 0 = single

Exp_Single <- function(x){
  result_single <- x/10
  return(result_single) 
} 

Exp_Fec_Mixed <- function(x){
  TW_mix <- (x/2) + (Side_By_Side_Replicates$Atlas_Avg_Fec/2) 
  Exp_Fec_mix <- TW_mix/10
  return(Exp_Fec_mix)
}


Side_By_Side_Replicates$Exp_Fec_Per_Plant <- ifelse(Side_By_Side_Replicates$Numbers == 1,
       NA,
       Exp_Single(Side_By_Side_Replicates$Avg_Fec))






### Functions to get expected per plant fitness for both conditions

Exp_Fit_Mixed <- function(x){
  fit_mix <- x + (Side_By_Side_Replicates$Atlas_Avg_Fitness/2)
  Exp_Fit_mix <- fit_mix/10
  return(Exp_Fit_mix)
}

Side_By_Side_Replicates$Exp_Fit_Per_Plant <- ifelse(Side_By_Side_Replicates$Numbers == 1,
                                                    Exp_Fit_Mixed(Side_By_Side_Replicates$Avg_Fit),
                                                    Exp_Single(Side_By_Side_Replicates$Avg_Fit))

### Functions to get expected total weight per plant for both conditions

Exp_TW_mix <- function(x){
  TW_mix <- x + (Side_By_Side_Replicates$Atlas_Avg_Total_Weight/2)
  Exp_TW <- TW_mix/10
  return(Exp_TW)
}

Side_By_Side_Replicates$Exp_TW_Per_Plant <- ifelse(Side_By_Side_Replicates$Numbers == 1,
                                                   Exp_TW_mix(Side_By_Side_Replicates$Avg_Total_Weight),
                                                   Exp_Single(Side_By_Side_Replicates$Avg_Total_Weight))

### Write Delims
write_delim(Full_Data, "Full_Data")
write_delim(Side_By_Side_Replicates, "Side_By_Side_Replicates")
write_delim(Rep_Mixed, "Rep_Mixed")
write_delim(Rep_Single, "Rep_Single")

