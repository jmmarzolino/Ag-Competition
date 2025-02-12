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

# Loading Data
pheno <- read_delim("DERIVED_PHENOTYPES.tsv")

#pheno <- pheno %>% 
#    select(c(Genotype, Replicate, Exp_year, FT, TOTAL_MASS, SEED_WEIGHT_100, FECUNDITY)) #%>% 
    #mutate(across(-c(Genotype, Generation, Replicate, Plants), ~(scale(.) %>% as.vector)))


## extract blups & H2 for each trait & model

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


#lmer(FT ~ 1 + (1|Genotype) + (1|Genotype:Exp_year), pheno) %>% plot

#mod <-  lmer(FT ~ 1 + Exp_year + (1|Genotype), pheno)
#pheno$resids <- resid(mod)
#pheno$fitted <- fitted(mod)
#ggplot(pheno, aes(fitted, resids, color=Plants)) + geom_point() + geom_smooth() + geom_hline(yintercept=0)


## Extract intercept + genotype data to re-scale data to original units
### trait average + genotypes deviance from average
#inter <- summary(h2_ft$model)$coefficients[1, 1]
#h2_ft$blups$FT <- h2_ft$blups$FT



################### TOTAL MASS
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


#lmer(TOTAL_MASS ~ 1 + Exp_year + Plants + Replicate + (1|Genotype) + (1|Genotype:Exp_year), pheno) 

#lmer(TOTAL_MASS ~ 1 + Exp_year + FITNESS + Generation + Replicate + (1|Genotype) + (1|Genotype:Exp_year), pheno)


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


#inter <- summary(h2_totmass$model)$coefficients[1, 1]
#h2_totmass$blups$TOTAL_MASS <- h2_totmass$blups$TOTAL_MASS
# interesting note in mass per plant - modeling with 'Plants' shows negative relationship btwn mass per plant and number of plants


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

#inter <- summary(h2_sw100$model)$coefficients[1, 1]
#h2_sw100$blups$SEED_WEIGHT_100 <- h2_sw100$blups$SEED_WEIGHT_100 + inter




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

#inter <- summary(h2_plants$model)$coefficients[1, 1]
#h2_plants$blups$Plants <- h2_plants$blups$Plants + inter



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

#mod <- lmer(FECUNDITY ~ 1 + Plants + FT + (1|Genotype) + (1|Genotype:Exp_year), pheno)
#pheno$resids <- resid(mod)
#pheno$fits <- fitted(mod)

#ggplot(pheno, aes(fitted, resids)) + geom_point() + geom_smooth() + geom_hline(yintercept=0)

#inter <- summary(h2_fec$model)$coefficients[1, 1]
#h2_fec$blups$FECUNDITY <- h2_fec$blups$FECUNDITY + inter

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


###### H2cal arguments ---
#gen.name: Name of the genotypes.
#rep.n: Number of replications in the experiment.
#year.n: Number of years (default = 1). See details.

#year.name: Name of the years (default = NULL). See details.

#fixed.model: The fixed effects in the model (BLUEs). See examples.
# emmeans: Use emmeans for calculate the BLUEs (default = FALSE).

#random.model: The random effects in the model (BLUPs). See examples.
# summary: Print summary from random model (default = FALSE).

#plot_diag: Show diagnostic plots for fixed and random effects (default = FALSE). Options: "base", "ggplot". .

#outliers.rm: Remove outliers (default = FALSE)

#For MET experiments you should ‘env.n’ and ‘env.name’ and/or ‘year.n’ and ‘year.name’ according your experiment.
