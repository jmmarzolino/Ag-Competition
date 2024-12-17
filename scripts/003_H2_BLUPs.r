#!/usr/bin/env Rscript
#SBATCH --mem=50G
#SBATCH --time=02:00:00
#SBATCH --output=/rhome/jmarz001/bigdata/Ag-Competition/scripts/003_BLUPs.stdout
#SBATCH -p koeniglab

library(tidyverse)
library(lme4)
library(inti)

setwd("/rhome/jmarz001/bigdata/Ag-Competition/data")
source("../scripts/CUSTOM_FNS.R")

# Loading Data
pheno <- read_delim("DERIVED_PHENOTYPES.tsv")

pheno_year_scale <- pheno %>% 
    filter(Condition == "single") %>% 
    select(c(Genotype, Generation, Replicate, Exp_year, FT, TOTAL_MASS, GERMINATION, SEED_WEIGHT_100, FECUNDITY, FITNESS)) %>% 
    mutate(across(-c(Genotype, Generation, Replicate, GERMINATION), ~(scale(.) %>% as.vector))) 


## extract blups & H2 for each trait & model


# init pdf to store plot outputs
pdf("heritability_model_fits.pdf", width=20, height=10)

################### FLOWERING TIME
h2_ft <- 
H2cal(data = pheno_year_scale
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


#lmer(FT ~ 1 + (1|Genotype) + (1|Genotype:Exp_year), pheno_year_scale) %>% plot

#mod <-  lmer(FT ~ 1 + Exp_year + (1|Genotype), pheno_year_scale)
#pheno_year_scale$resids <- resid(mod)
#pheno_year_scale$fitted <- fitted(mod)
#ggplot(pheno_year_scale, aes(fitted, resids, color=GERMINATION)) + geom_point() + geom_smooth() + geom_hline(yintercept=0)



################### TOTAL MASS
h2_totmass <- 
H2cal(data = pheno_year_scale
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


#lmer(TOTAL_MASS ~ 1 + Exp_year + GERMINATION + Replicate + (1|Genotype) + (1|Genotype:Exp_year), pheno_year_scale) 

#lmer(TOTAL_MASS ~ 1 + Exp_year + FITNESS + Generation + Replicate + (1|Genotype) + (1|Genotype:Exp_year), pheno_year_scale)





h2_sw100 <- 
H2cal(data = pheno_year_scale
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
H2cal(data = pheno_year_scale
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
H2cal(data = pheno_year_scale
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

#mod <- lmer(FECUNDITY ~ 1 + GERMINATION + FT + (1|Genotype) + (1|Genotype:Exp_year), pheno_year_scale)
#pheno_year_scale$resids <- resid(mod)
#pheno_year_scale$fits <- fitted(mod)

#ggplot(pheno_year_scale, aes(fitted, resids)) + geom_point() + geom_smooth() + geom_hline(yintercept=0)






h2_fit <- 
H2cal(data = pheno_year_scale
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


H2_table <- tibble("trait" = colnames(pheno_year_scale)[5:ncol(pheno_year_scale)], 
                    #"H2" = c(h2_ft$tabsmr$H2.c, h2_totmass$tabsmr$H2.c, h2_sw100$tabsmr$H2.c, h2_germ$tabsmr$H2.c, h2_fec$tabsmr$H2.c, h2_fit$tabsmr$H2.c),
                    "H2_s" = c(h2_ft$tabsmr$H2.s, h2_totmass$tabsmr$H2.s, h2_sw100$tabsmr$H2.s, h2_germ$tabsmr$H2.s, h2_fec$tabsmr$H2.s, h2_fit$tabsmr$H2.s)
                    )
write_delim(H2_table, "trait_heritability.tsv", "\t")



# Creating BLUP dataframe
blup_output <- full_join(h2_ft$blups, h2_totmass$blups)
blup_output <- full_join(blup_output, h2_sw100$blups)
blup_output <- full_join(blup_output, h2_germ$blups)
blup_output <- full_join(blup_output, h2_fec$blups)
blup_output <- full_join(blup_output, h2_fit$blups)

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
