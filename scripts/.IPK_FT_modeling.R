#!/usr/bin/env Rscript

#SBATCH --ntasks=1
#SBATCH --mem=10G
#SBATCH --time=02:00:00
#SBATCH --job-name="FT correlation"
#SBATCH --output=/rhome/jmarz001/bigdata/IPK_Analysis/scripts/000_FT_trait_corr.stdout
#SBATCH -p short

setwd("/bigdata/koeniglab/jmarz001/IPK_Analysis/")
library(tidyverse)
library(corrplot)
library(lme4)

#########   MODEL FLOWERING TIME TO ONE VALUE PER GENOTYPE
FT <- read_delim("/rhome/jmarz001/bigdata/Ag-Competition/data/FT_BOTH_YEARS.tsv")
FT <- FT %>% filter(Condition=="single") %>% select(-c(Plants, Condition))
FT_avg <- FT %>% group_by(Exp_year, Genotype, Generation) %>% summarise(avg_FT = mean(FT))
FT_avg_2yr <- FT %>% group_by(Genotype, Generation) %>% summarise(avg_FT = mean(FT))
write_delim(FT_avg_2yr, "FT_AVG_2YRS.tsv")


FT$dummy_yr <- 1
FT[which(FT$Exp_year == 2023), 6] <- 2
# FT[which(FT$Genotype == "1_100"),]

model <- lmer(FT ~ dummy_yr + (1| Genotype), data=FT)
#lmer(FT ~ dummy_yr + (1| Genotype), data=FT) %>% coef

# divide up and store genotypes modelled FT
model_coefs <- coef(model)
model_coefs <- model_coefs$Genotype
model_coefs <- rownames_to_column(model_coefs)
colnames(model_coefs) <- c("Genotype", "Intercept", "Year_Coef")
model_coefs$FT <- model_coefs$Intercept + model_coefs$Year_Coef
model_coefs <- model_coefs %>% select(c(Genotype, FT)) %>% tibble
write_delim(model_coefs, "data/modelled_FT.tsv")
write_delim(model_coefs, "../Ag-Competition/data/modelled_FT.tsv")


#######   NOW CORRELATE THE FT W OTHER TRAITS
data <- read_delim("data/joined_derived_traits_data_GenotypeFiltered.tsv")
# check that genotypes were collapsed/summarised
which(table(data$Genotype) > 1)

# match genotype format (\\d_\\d+_*\\d*)
# remove UCRKL and leading zeros, sub - w _
data$Genotype <- gsub("UCRKL0000", "", data$Genotype)
data$Genotype <- gsub("0", "", data$Genotype)
data$Genotype <- gsub("-", "_", data$Genotype)

# put both dfs in order of genotype
data <- data[order(data$Genotype), ]
model_coefs <- model_coefs[order(model_coefs$Genotype), ]

# try to standardize some more genotyeps (parents)
#data[grep("^\\d{2}_\\d.*", data$Genotype),1]
data$Genotype <- gsub("^(\\d{2})_\\d.*", "\\1", data$Genotype)
model_coefs$Genotype <- gsub("^(\\d{2})_\\d.*", "\\1", model_coefs$Genotype)

# capture a few more genotypes by removing third numbers (out of #_#_#)
data$Genotype <- gsub("^(\\d+_\\d+)_\\d+$", "\\1", data$Genotype)
model_coefs$Genotype <- gsub("^(\\d+_\\d+)_\\d+$", "\\1", model_coefs$Genotype)


## line up data genotypes w FT genotypes

#full_join(data[,1:2], model_coefs)
#left_join(model_coefs, data[,1:4])
#setdiff(model_coefs$Genotype, data$Genotype)
#setdiff(data$Genotype, model_coefs$Genotype)
df <- inner_join(data, model_coefs, by="Genotype")
apply(df[,3:ncol(df)], 2, sd, na.rm=T)
cor(df$FT, df[, 3:ncol(df)])
cor(df$FT, df[, 3:ncol(df)]) %>% quantile(na.rm=T)



by_gen <- df %>% group_by(Generation) %>% summarise(across(where(is.numeric), mean))
#apply(by_gen[,3:ncol(by_gen)], 2, sd, na.rm=T)
cor(by_gen$FT, by_gen[, 3:ncol(by_gen)])
cor(by_gen$FT, by_gen[, 3:ncol(by_gen)]) %>% quantile(na.rm=T)

by_gen_cor <- cor(by_gen$FT, by_gen[, 3:ncol(by_gen)])
by_gen_cor[which(abs(by_gen_cor)>0.6)]



x <- tibble("FT_correlation" = c(by_gen_cor))
x$trait <- dimnames(by_gen_cor)[2][[1]]

x %>% filter(abs(FT_correlation)>0.6)

# combine FT correlation values with table of sig trait change
sig <- read_delim("results/002A_change_over_time_derived_traits.tsv")


tmp <- full_join(x, sig, by="trait")
write_delim(tmp, "results/002A_change_over_time_derived_traits_FT_corr.tsv", "\t")
