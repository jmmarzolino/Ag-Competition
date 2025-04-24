#!/usr/bin/env Rscript
#SBATCH --mem=30G
#SBATCH --time=02:00:00
#SBATCH --output=/rhome/jmarz001/bigdata/Ag-Competition/scripts/004_combine_exp_data.stdout
#SBATCH -p koeniglab

library(tidyverse)
library(ggpubr)
library(data.table)
library(corrplot)

setwd("/bigdata/koeniglab/jmarz001/Ag-Competition/data")
source("../scripts/CUSTOM_FNS.R")

#####    Greenhouse Exp
# CC II greenhouse experiment data downloaded from paper supplement
df <- fread("Prog_pheno.txt")
prt <- fread("Parents_pheno.txt")

## explore data
# how many lines are 2 vs. 6 row
table(round(df$rows))
# Fam_1 = generation number
df %>% group_by(Fam_1) %>% summarise(across(c(days_to_heading, days_to_heading_2017), \(x) mean(x, na.rm=T)))

sum(is.na(df$days_to_heading))
sum(is.na(df$days_to_heading_2017))
# 2017 is missing >60% data
# I'm going to focus on just days_to_heading until I've considered methods: averaging values, modeling, etc.

# select only necessary columns
df <- df %>% 
    select(c(V1, Samples, X100_seed_mass, seed_mass, seed_estimate, days_to_heading))
# add ID number to parent data
prt$Samples <- as.character(45:72)
# select necessary columns
prt <- prt %>% 
    select(c(V1, Samples, X100_seed_mass, seed_mass, seed_estimate, days_to_heading))

# join parent and progeny data from greenhouse experiment
gh <- full_join(df, prt) %>% select(-V1)




#####    Field Exp
# CC II field experiment data
ag <- fread("trait_BLUPs.tsv") 

# format genotype IDs to 'family' groups, #_##, removing any third numbers
ag$Genotype <- gsub("^(\\d+_\\d+).*", "\\1", ag$Genotype)
# for parent 'generation' numbers (>7),
# format genotype to only first number
ag$gen <- (as.numeric(gsub("^(\\d+)_\\d+.*", "\\1", ag$Genotype)))
ag[which(ag$gen > 10), 1] <- ag[which(ag$gen > 10), 8]


####    Join data
# join data by genotypes in both data sets
combo <- inner_join(ag, gh, by=c("Genotype"="Samples")) %>% select(-gen)
write_delim(combo, "combined_CCII.tsv")


####    Trait correlations
# save table of trait correlation values
x_cor <- cor(combo[,2:ncol(combo)], use="na.or.complete", method="spearman")
dimnames(x_cor)[[1]] <- tidy_text_substitution(dimnames(x_cor)[[1]])
dimnames(x_cor)[[2]] <- tidy_text_substitution(dimnames(x_cor)[[2]])
x <- data.table(x_cor)
x$TRAIT <- colnames(x)
write_delim(tibble(x), "over_exps_correlation_table.tsv", "\t")



#####    Plotting
#plot(combo[,2:10])

# trait pairs for plotting
#"X100_seed_mass"="SEED_WEIGHT_100", 
#"seed_mass"="TOTAL_MASS"
#"seed_estimate"="FECUNDITY", 
#"days_to_heading"="FT"

# plot trait relationships 
png("over_exps_trait_corrs.png")
corrplot(x_cor, method="color", type="upper", order="original", title="", mar=c(0,0,4,0), addCoef.col = "black", tl.col = "black")
dev.off()
