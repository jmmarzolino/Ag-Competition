#!/usr/bin/env Rscript

library(tidyverse)
library(lme4)
library(nlme)
library(inti)

setwd("/rhome/jmarz001/bigdata/Ag-Competition/data")

# Loading Data
pheno <- read_delim("DERIVED_PHENOTYPES.tsv")
# set experiment year as factor
pheno$Exp_year <- as.factor(pheno$Exp_year)


# model with inter::H2cal
h2_totmass <- 
H2cal(data = pheno
          , trait = "TOTAL_MASS"
          , gen.name = "Genotype"
          , rep.n = 2
          , year.n = 2
          , year.name = "Exp_year"
          , fixed.model = "1 + Genotype + Exp_year + (1|Genotype:Exp_year)"
          , random.model = "1 + (1|Genotype) + Exp_year + (1|Genotype:Exp_year)"
          , plot_diag = TRUE
          , outliers.rm = TRUE
          , summary = TRUE
          )


# model with lme4::lmer
lmer(TOTAL_MASS ~ 1 + (1 |Genotype) + Exp_year + (1|Genotype:Exp_year), pheno)
# plot fitted vs residual values
lmer(TOTAL_MASS ~ 1 + (1 |Genotype) + Exp_year + (1|Genotype:Exp_year), pheno) %>% plot

## total-mass model residuals colored by experimental year 
mod <-  lmer(TOTAL_MASS ~ 1 + (1 |Genotype) + Exp_year + (1|Genotype:Exp_year), pheno) 
# extract coresponding col of experiment year
yr <- pheno[as.numeric(names(fitted(mod))), 2]
# combine fitted, residuals, and year & plot
mod_tb <- tibble("year"=as.factor(yr$Exp_year), "residuals"=resid(mod), "fitted"=fitted(mod))
ggplot(mod_tb, aes(fitted, residuals, color=year)) + geom_point() + geom_smooth() + geom_hline(yintercept=0)



# model with inter::lme
# lme handles NA values differently, so filter NA values from phenotype
pheno <- pheno[which(!is.na(pheno$TOTAL_MASS)), ]
# pheno ~ 1+ fixed effects, (pheno) ~ 1 | random effects
lme(TOTAL_MASS ~ 1 + Exp_year, random = ~ 1 | Genotype, data=pheno)

lme(TOTAL_MASS ~ 1 + Exp_year, random = ~ 1 | Genotype, data=pheno) %>% plot
