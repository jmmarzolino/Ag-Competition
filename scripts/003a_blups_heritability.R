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
    #mutate(across(-c(Genotype, Generation, Replicate, Germination), ~(scale(.) %>% as.vector)))
pheno$Exp_year <- as.factor(pheno$Exp_year)

# init pdf to store plot outputs
pdf("heritability_model_fits.pdf", width=20, height=10)


h2_ft <- 
H2cal(data = pheno
          , trait = "FT"
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


h2_massper <- 
H2cal(data = pheno
          , trait = "MASS_PER_PLANT"
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


h2_sw100 <- 
H2cal(data = pheno
          , trait = "SEED_WEIGHT_100"
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


## not well fitable
h2_germ <- 
H2cal(data = pheno
        , trait = "Germination"
        , gen.name = "Genotype"
        , rep.n = 2
        , year.n = 2
        , year.name = "Exp_year"
        , fixed.model = "1 + Genotype + Exp_year + (1|Genotype:Exp_year)"
        , random.model = "1 + (1|Genotype) + Exp_year + (1|Genotype:Exp_year)"
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
          , fixed.model = "1 + Genotype + Exp_year + (1|Genotype:Exp_year)"
          , random.model = "1 + (1|Genotype) + Exp_year + (1|Genotype:Exp_year)"
          , plot_diag = TRUE
          , outliers.rm = TRUE
          )

dev.off()


H2_table <- tibble("trait" = c("FT", "TOTAL_MASS", "SEED_WEIGHT_100", "MASS_PER_PLANT", "Germination", "FECUNDITY"), 
                    "H2_s" = c(h2_ft$tabsmr$H2.s, h2_totmass$tabsmr$H2.s, h2_sw100$tabsmr$H2.s, h2_massper$tabsmr$H2.s,  h2_germ$tabsmr$H2.s, 
                    h2_fec$tabsmr$H2.s) )
write_delim(H2_table, "trait_heritability.tsv", "\t")


# combine BLUP dataframes
blup_output <- full_join(h2_ft$blups, h2_totmass$blups)
blup_output <- full_join(blup_output, h2_sw100$blups)
blup_output <- full_join(blup_output, h2_massper$blups)
blup_output <- full_join(blup_output, h2_germ$blups)
blup_output <- full_join(blup_output, h2_fec$blups)

write_delim(blup_output, "trait_BLUPs.tsv", "\t")


#lmer(SEED_WEIGHT_100 ~ 1 + (1|Genotype) + Exp_year + (1|Genotype:Exp_year), pheno) %>% plot
#mod <-  lmer(SEED_WEIGHT_100 ~ 1 + (1|Genotype) + Exp_year+ (1|Genotype:Exp_year), pheno)
# yr <- pheno[as.numeric(names(fitted(mod))), 2]
#mod_tb <- tibble("year"=as.factor(yr$Exp_year), "residuals"=resid(mod), "fitted"=fitted(mod))
#ggplot(mod_tb, aes(fitted, residuals, color=year)) + geom_point() + geom_smooth() + geom_hline(yintercept=0)

# write out edited genotype list without parental genotypes for filtering vcf
blup_output <- add_generation(blup_output)
blup_output$Genotype <- gsub("^(\\d+_\\d+)_\\d", "\\1", blup_output$Genotype)
geno_list <- blup_output[which(blup_output$Generation != 0), 1][[1]]
geno_tab <- tibble("X1"=geno_list, "X2"=geno_list)
write_delim(geno_tab, "../results/gwas/AgComp_genotypes.tsv", col_names=F)