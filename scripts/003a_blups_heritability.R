#!/usr/bin/env Rscript
#SBATCH --mem=50G
#SBATCH --time=02:00:00
#SBATCH --output=/rhome/jmarz001/bigdata/Ag-Competition/scripts/003a_blups_heritability.stdout
#SBATCH -p koeniglab

library(tidyverse)
library(lme4)
library(inti)

setwd("/rhome/jmarz001/bigdata/Ag-Competition/data")
source("../scripts/CUSTOM_FNS.R")

## extract blups & H2 for each trait & model

# Loading Data
pheno <- read_delim("DERIVED_PHENOTYPES.tsv")
#pheno <- pheno %>% 
#    select(c(Genotype, Replicate, Exp_year, FT, TOTAL_MASS, SEED_WEIGHT_100, FECUNDITY)) #%>% 
    #mutate(across(-c(Genotype, Generation, Replicate, Plants), ~(scale(.) %>% as.vector)))

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
          , fixed.model = "0 + (1|Exp_year) + Genotype + (1|Genotype:Exp_year)"
          , random.model = "1 + (1|Exp_year) + (1|Genotype) + (1|Genotype:Exp_year)"
          , plot_diag = TRUE
          , outliers.rm = TRUE
          , summary = TRUE
          )


h2_totmass <- 
H2cal(data = pheno
          , trait = "TOTAL_MASS"
          , gen.name = "Genotype"
          , rep.n = 2
          , year.n = 2
          , year.name = "Exp_year"
          , fixed.model = "0 + (1|Exp_year) + Plants + (1|Replicate) + Genotype + (1|Genotype:Exp_year)"
          , random.model = "1 + (1|Exp_year) + Plants + (1|Replicate)  + (1|Genotype) + (1|Genotype:Exp_year)"
          , plot_diag = TRUE
          , outliers.rm = TRUE
          , summary = TRUE
          )

h2_massper <- 
H2cal(data = pheno
          , trait = "MASS_PER_PLANT"
          , gen.name = "Genotype"
          , rep.n = 2
          , year.n = 2
          , year.name = "Exp_year"
          , fixed.model = "0 + (1|Exp_year) + Plants + (1|Replicate) + Genotype + (1|Genotype:Exp_year)"
          , random.model = "1 + (1|Exp_year) + Plants + (1|Replicate)  + (1|Genotype) + (1|Genotype:Exp_year)"
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
          , fixed.model = "0 + TOTAL_MASS + (1|Exp_year) + Genotype + FT + (1|Genotype:Exp_year)"
          , random.model = "1 + TOTAL_MASS + (1|Exp_year) + (1|Genotype) + FT + (1|Genotype:Exp_year)"
          , plot_diag = TRUE
          , outliers.rm = TRUE
          , summary = TRUE
          )


## not fitable
h2_plants <- 
H2cal(data = pheno
        , trait = "Plants"
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
          , fixed.model = "0 + Plants + FT + Genotype + (1|Genotype:Exp_year)"
          , random.model = "1 + Plants + FT + (1|Genotype) + (1|Genotype:Exp_year)"
          , plot_diag = TRUE
          , outliers.rm = TRUE
          )


dev.off()


H2_table <- tibble("trait" = c("FT", "TOTAL_MASS", "SEED_WEIGHT_100", "MASS_PER_PLANT", "Plants", "FECUNDITY"), 
                    "H2_s" = c(h2_ft$tabsmr$H2.s, h2_totmass$tabsmr$H2.s, h2_sw100$tabsmr$H2.s, h2_massper$tabsmr$H2.s,  h2_plants$tabsmr$H2.s, 
                    h2_fec$tabsmr$H2.s) )
write_delim(H2_table, "trait_heritability.tsv", "\t")


# combine BLUP dataframes
blup_output <- full_join(h2_ft$blups, h2_totmass$blups)
blup_output <- full_join(blup_output, h2_sw100$blups)
blup_output <- full_join(blup_output, h2_massper$blups)
blup_output <- full_join(blup_output, h2_plants$blups)
blup_output <- full_join(blup_output, h2_fec$blups)

write_delim(blup_output, "trait_BLUPs.tsv", "\t")

