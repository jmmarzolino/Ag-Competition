## TEMP script
library(tidyverse)
library(data.table)
library(corrplot)

setwd("/rhome/jmarz001/bigdata/Ag-Competition/data")
source("../scripts/CUSTOM_FNS.R")

df <- fread("DERIVED_PHENOTYPES.tsv")
df <- df %>% 
        filter(Condition == "single") %>% 
        select(-Condition) %>%
        group_by(Genotype) %>%
        summarise(across(where(is.numeric), mean)) %>%
        ungroup() %>%
        select(-c(Replicate, Exp_year, SEED_COUNT))

blup <- fread("trait_BLUPs.tsv")
test <- full_join(df, blups)

cor(test$FT, test$FT_blup)
cor(test$SEED_WEIGHT_100, test$SEED_WEIGHT_100_blup)
cor(test$TOTAL_MASS, test$TOTAL_MASS_blup)
cor(test$Germination, test$Germination_blup)
cor(test$FECUNDITY, test$FECUNDITY_blup)
cor(test$MASS_PER_PLANT, test$MASS_PER_PLANT_blup)

## all traits have high correlation from averged genotype trait values w genotype's blup values

# lowest value is 0.94 for 100-seed weight
# so changing avg values for blups should not substantially change things