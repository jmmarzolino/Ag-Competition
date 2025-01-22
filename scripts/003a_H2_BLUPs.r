#!/usr/bin/env Rscript
#SBATCH --mem=50G
#SBATCH --time=02:00:00
#SBATCH --output=/rhome/jmarz001/bigdata/Ag-Competition/scripts/003a_H2_BLUPs.stdout
#SBATCH -p koeniglab

library(tidyverse)
library(lme4)
library(inti)

setwd("/rhome/jmarz001/bigdata/Ag-Competition/data")
source("../scripts/CUSTOM_FNS.R")
## extract blups & H2 for each trait & model

# load & format data
pheno <- read_delim("DERIVED_PHENOTYPES.tsv")

pheno <- pheno %>% 
    filter(Condition == "single") %>% 
    select(c(Genotype, Generation, Replicate, Exp_year, FT, TOTAL_MASS, GERMINATION, SEED_WEIGHT_100, FECUNDITY, FITNESS)) 

pheno$Exp_year <- as.factor(pheno$Exp_year)


# init pdf to store plot outputs
pdf("heritability_model_fits.pdf", width=20, height=10)

################### FLOWERING TIME
h2_ft <- 
H2cal(data = pheno
          , trait = "FT"
          , gen.name = "Genotype"
          , rep.n = 2
          , year.n = 2
          , year.name = "Exp_year"
          , fixed.model = "0 + Exp_year + Genotype + (1|Genotype:Exp_year)"
          , random.model = "1 + Exp_year + (1|Genotype) + (1|Genotype:Exp_year)"
          , plot_diag = TRUE
          , outliers.rm = TRUE
          , summary = TRUE
          )


################### TOTAL MASS
h2_totmass <- 
H2cal(data = pheno
          , trait = "TOTAL_MASS"
          , gen.name = "Genotype"
          , rep.n = 2
          , year.n = 2
          , year.name = "Exp_year"
          , fixed.model = "0 + GERMINATION + (1|Replicate) + Genotype + (1|Genotype:Exp_year)"
          , random.model = "1 + GERMINATION + (1|Replicate)  + (1|Genotype) + (1|Genotype:Exp_year)"
          , plot_diag = TRUE
          , outliers.rm = TRUE
          , summary = TRUE
          )


h2_sw100 <- 
H2cal(data = pheno
          , trait = "SEED_WEIGHT_100"
          , gen.name = "Genotype"
          , rep.n = 2
          , year.n = 2
          , year.name = "Exp_year"
          , fixed.model = "0 + Genotype + (1|Genotype:Exp_year)"
          , random.model = "1 + (1|Genotype) + (1|Genotype:Exp_year)"
          , plot_diag = TRUE
          , outliers.rm = TRUE
          , summary = TRUE
          )


## not fitable
h2_germ <- 
H2cal(data = pheno
        , trait = "GERMINATION"
        , gen.name = "Genotype"
        , rep.n = 2
        , year.n = 2
        , year.name = "Exp_year"
        , fixed.model = "0 + (1|Replicate) + Genotype + (1|Genotype:Exp_year)"
        , random.model = "1 + (1|Replicate) + (1|Genotype) + (1|Genotype:Exp_year)"
        , plot_diag = TRUE
        , outliers.rm = TRUE
        )


h2_fec <- 
H2cal(data = pheno
          , trait = "FECUNDITY"
          , gen.name = "Genotype"
          , rep.n = 2
          , year.n = 2
          , year.name = "Exp_year"
          , fixed.model = "0 + GERMINATION + FT + Genotype + (1|Genotype:Exp_year)"
          , random.model = "1 + GERMINATION + FT + (1|Genotype) + (1|Genotype:Exp_year)"
          , plot_diag = TRUE
          , outliers.rm = TRUE
          )


h2_fit <- 
H2cal(data = pheno
        , trait = "FITNESS"
        , gen.name = "Genotype"
        , rep.n = 2
        , year.n = 2
        , year.name = "Exp_year"
        , fixed.model = "0 + FT + Genotype + (1|Genotype:Exp_year)"
        , random.model = "1 + FT + (1|Genotype) + (1|Genotype:Exp_year)"
        , plot_diag = TRUE
        , outliers.rm = TRUE
        )

dev.off()


H2_table <- tibble("trait" = colnames(pheno)[5:ncol(pheno)], 
                    "H2_s" = c(h2_ft$tabsmr$H2.s, h2_totmass$tabsmr$H2.s, h2_sw100$tabsmr$H2.s, h2_germ$tabsmr$H2.s, h2_fec$tabsmr$H2.s, h2_fit$tabsmr$H2.s) )
write_delim(H2_table, "trait_heritability.tsv", "\t")


# combine BLUP dataframes
blup_output <- full_join(h2_ft$blups, h2_totmass$blups)
blup_output <- full_join(blup_output, h2_sw100$blups)
blup_output <- full_join(blup_output, h2_germ$blups)
blup_output <- full_join(blup_output, h2_fec$blups)
blup_output <- full_join(blup_output, h2_fit$blups)

write_delim(blup_output, "trait_BLUPs.tsv", "\t")

## Extract intercept + genotype data to re-scale data to original units
### trait average + genotypes deviance from average
#inter <- summary(h2_ft$model)$coefficients[1, 1]
#h2_ft$blups$FT <- h2_ft$blups$FT + inter
## OK SO ACTUALLY
# H2cal DOES RE-SCALE blups to the intercept
#mean(h2_ft$blups$FT - inter)
# the average of the blups minus the intercept is near 0,
# indicating those values are the indv. genotypes perturbations from the intercept
# so the blup values extracted from h2_ft$blups$FT ALREADY HAVE THE INTERCEPT ADDED

