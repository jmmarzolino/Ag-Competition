#!/usr/bin/env Rscript
#SBATCH --mem=50G
#SBATCH --time=02:00:00
#SBATCH --output=/rhome/jmarz001/bigdata/Ag-Competition/scripts/004_heritability.stdout
#SBATCH -p koeniglab

library(tidyverse)
library(lme4)
#library(methods)
library(ggpubr)

setwd("/rhome/jmarz001/bigdata/Ag-Competition/data")

# Loading Data
pheno <- read_delim("DERIVED_PHENOTYPES.tsv")
# select columns and scale phenotype data
pheno <- pheno %>% filter(Condition == "single") %>% select(c(Genotype, Generation, Condition, Replicate, Exp_year, FT, TOTAL_MASS, SURVIVAL, SEED_WEIGHT_100, FECUNDITY, FITNESS)) %>% mutate(across(-c(Genotype, Generation, Condition, Replicate, Exp_year), ~(scale(.) %>% as.vector))) 



# check model fits for different traits 
for(i in 6:ncol(pheno)) {
  print(i)
  print(colnames(pheno[,i]))
  print(AIC(lmer(as.numeric(unlist(pheno[,i])) ~ Exp_year + Replicate + Generation + (1|Genotype), data=pheno)))
  print(AIC(lmer(as.numeric(unlist(pheno[,i])) ~ Exp_year + Replicate + (1|Genotype), data=pheno)))
  print(AIC(lmer(as.numeric(unlist(pheno[,i])) ~ Exp_year + (1|Exp_year:Replicate) + (1|Genotype), data=pheno)))
  print(AIC(lmer(as.numeric(unlist(pheno[,i])) ~ Exp_year + Generation + (1|Genotype), data=pheno)))
  print(AIC(lmer(as.numeric(unlist(pheno[,i])) ~ Exp_year + (1|Generation:Genotype), data=pheno)))
  print(AIC(lmer(as.numeric(unlist(pheno[,i])) ~ Exp_year + (1+Generation|Genotype), data=pheno)))
  print(AIC(lmer(as.numeric(unlist(pheno[,i])) ~ Exp_year + (1|Genotype), data=pheno)))
}
# Exp_year + (1|Genotype) is a winner across the board




# Breeders function for calculating heritability

Breeders_funct <- function(x) {
  # set up empty df for storing results
  var_tab <- tibble(trait= character(), genetic_var=double(), total_var=double())

      # calculate for each trait
      for(i in c(6:ncol(x))){
          df <- x %>%
              select(c(all_of(i), Genotype, Replicate, Generation, Exp_year))
            # calculate mixed model
            genotype_mod <- lmer(df[,1][[1]] ~ Exp_year + (1| Genotype), data = df)
            ## summary of the model
            geno_mod_summ <- summary(genotype_mod)

            ## genotypic variance
            genotypic_variance <- (data.frame(geno_mod_summ$varcor)$sdcor[1]^2)
            ## total variance
            total_variance <- genotypic_variance + (data.frame(geno_mod_summ$varcor)$sdcor[2]^2)
            new_row <- tibble(trait= colnames(x[,i]), genetic_var=genotypic_variance, total_var=total_variance)
            var_tab <- bind_rows(var_tab, new_row)
          }

          ## ratio of Vg and Vp
          var_tab$H2 <- var_tab$genetic_var / var_tab$total_var
          return(var_tab)
  }

heritability <- Breeders_funct(pheno)
# note - survival isn't well described by the model, probably because of how it was filtered to only plots with high enough survival
# model fit is singular
write_delim(heritability, "trait_heritability.tsv", "\t")




# Extracting Broad Sense Heritability per generation
h0 <- pheno %>% filter(Generation == 0) %>% Breeders_funct()
h18 <- pheno %>% filter(Generation == 18) %>% Breeders_funct()
h28 <- pheno %>% filter(Generation == 28) %>% Breeders_funct()
h50 <- pheno %>% filter(Generation == 50) %>% Breeders_funct()
h58 <- pheno %>% filter(Generation == 58) %>% Breeders_funct()
herit <- rbind(h0, h18, h28, h50, h58)
herit$Generation <- c(rep(0, 6), rep(18, 6), rep(28, 6), rep(50, 6), rep(58,6))
write_delim(herit, "trait_heritability_per_generation.tsv", "\t")
# Heritability over time

ggplot(herit, aes(Generation, H2, color = trait)) +
  geom_point() +
  geom_smooth(alpha = .4) +
  labs(y = "Broad Sense Heritability",
      title = "Generational Change in Broad Sense Heritability")
ggsave("/bigdata/koeniglab/jmarz001/Ag-Competition/results/H2_over_gens.png", width = 14, height = 10)

# Genetic Variance over time

a <- ggplot(herit, aes(Generation, genetic_var, color = trait)) +
  geom_point() +
  geom_smooth() +
  labs(y = "Genetic Variance",
      title = "Generational Change in Genetic Variance") +
  ylim(0, 2.5)

# Phenotypic Variance over time

a1 <- ggplot(herit, aes(Generation, total_var, color = trait)) +
  geom_point() +
  geom_smooth() +
  ylim(0, 2.5) +
  labs(y = "Phenotypic Variance",
      title = "Generational Change in Phenotypic Variance")

y <- ggarrange(a, a1, ncol =2, nrow = 1)
ggsave("/bigdata/koeniglab/jmarz001/Ag-Competition/results/genetic_phenotypic_var_over_gens.png", y, width = 22, height = 12)

#herit2 <- herit %>% group_by(trait) %>% summarise(across(where(is.numeric), .fns = mean, na.rm = T)) %>% select(-c("Generation")) %>% ungroup()







# Creating BLUP dataframe
blup_output <- tibble("Genotype" = unique(pheno$Genotype), )

# model traits and extract random effect of genotype
for (i in 6:ncol(pheno)){
      model <- lmer(as.numeric(unlist(pheno[,i])) ~ Exp_year + (1|Genotype), data = pheno)
      summary(model)
      model_ranefs <- ranef(model)
      geno_blups <- model_ranefs[[1]][[1]]
      geno_blups <- tibble(geno_blups)
      geno_blups$Genotype <- rownames(model_ranefs[[1]])
      colnames(geno_blups)[1] <- paste(colnames(pheno[,i]),"blup",sep="_")
      blup_output <- full_join(blup_output, geno_blups)
}

write_delim(blup_output, "trait_BLUPs.tsv", "\t")









    # Response to Selection

    subtract <- function(x){
      x[2] - x[1]
    }

    # P to F58

    preserved <- phenos %>% filter(Generation == 0 | Generation == 58)
    P_F58 <- preserved %>% group_by(Generation) %>% summarise(across(.cols = c(SEED_WEIGHT_100, FECUNDITY, ABS_FITNESS, TOTAL_MASS, FT), mean, na.rm = T)) %>% ungroup()
    P_F58 <- P_F58 %>% summarise(across(.cols = c(SEED_WEIGHT_100, FECUNDITY, ABS_FITNESS, TOTAL_MASS, FT), subtract)) %>% ungroup()
    P_F58 <- P_F58 %>% mutate(SEED_WEIGHT_100 = SEED_WEIGHT_100/58,
                              FECUNDITY = FECUNDITY/58,
                              ABS_FITNESS = ABS_FITNESS/58,
                              FT = FT/58,
                              TOTAL_MASS = TOTAL_MASS/58,
                              response_years = "P_F58")
    # P to F18

    preserved <- phenos %>% filter(Generation == 0 | Generation == 18)
    P_F18 <- preserved %>% group_by(Generation) %>% summarise(across(.cols = c(SEED_WEIGHT_100, FECUNDITY, ABS_FITNESS, TOTAL_MASS, FT), mean, na.rm = T)) %>% ungroup()
    P_F18 <- P_F18 %>% summarise(across(.cols = c(SEED_WEIGHT_100, FECUNDITY, ABS_FITNESS, TOTAL_MASS, FT), subtract)) %>% ungroup()
    P_F18 <- P_F18 %>% mutate(SEED_WEIGHT_100 = SEED_WEIGHT_100/18,
                              FECUNDITY = FECUNDITY/18,
                              ABS_FITNESS = ABS_FITNESS/18,
                              FT = FT/18,
                              TOTAL_MASS = TOTAL_MASS/18,
                              response_years = "P_F18")

    # F18 to F58

    preserved <- phenos %>% filter(Generation == 18 | Generation == 58)
    F18_F58 <- preserved %>% group_by(Generation) %>% summarise(across(.cols = c(SEED_WEIGHT_100, FECUNDITY, ABS_FITNESS, TOTAL_MASS, FT), mean, na.rm = T)) %>% ungroup()
    F18_F58 <- F18_F58 %>% summarise(across(.cols = c(SEED_WEIGHT_100, FECUNDITY, ABS_FITNESS, TOTAL_MASS, FT), subtract)) %>% ungroup()
    F18_F58 <- F18_F58 %>% mutate(SEED_WEIGHT_100 = SEED_WEIGHT_100/40,
                              FECUNDITY = FECUNDITY/40,
                              ABS_FITNESS = ABS_FITNESS/40,
                              FT = FT/40,
                              TOTAL_MASS = TOTAL_MASS/40,
                              response_years = "F18_F58")

    rts_df <- rbind(P_F18, P_F58, F18_F58)

    a1 <- ggplot(rts_df, aes(response_years, TOTAL_MASS)) +
      geom_point() +
      ylim(-.03, .030)

    a2 <- ggplot(rts_df, aes(response_years, SEED_WEIGHT_100)) +
      geom_point() +
      ylim(-.03, .03)

    a3 <- ggplot(rts_df, aes(response_years, FT)) +
      geom_point() +
      ylim(-.03, .03)

    a4 <- ggplot(rts_df, aes(response_years, FECUNDITY)) +
      geom_point() +
      ylim(-.03, .03)

    a5 <- ggplot(rts_df, aes(response_years, ABS_FITNESS)) +
      geom_point() +
      ylim(-.03, .03)

    g <- ggarrange(a1, a2, a3, a4, a5, top = "Response to Selection", nrow = 1, ncol = 5)
    ggsave("/bigdata/koeniglab/jmarz001/Ag-Competition/results/Response_to_Selection_Single.png", g, width = 22, height = 10)

    # Selection

    rts_df <- rts_df %>% mutate(H2_FT = herit$H2[4],
                                H2_100 = herit$H2[1],
                                H2_tw = herit$H2[5],
                                H2_fec = herit$H2[3],
                                H2_fit = herit$H2[2])

    rts_df <- rts_df %>% mutate(si_100 = SEED_WEIGHT_100/ H2_100,
                                si_fec = FECUNDITY / H2_fec,
                                si_ft = FT/H2_FT,
                                si_tw = TOTAL_MASS/H2_tw,
                                si_fit = ABS_FITNESS/H2_fit)

    selection_intensity_single <- rts_df %>% select(starts_with("si"), response_years)

    # 100 SW SI

    a <- ggplot(selection_intensity_single, aes(response_years, si_100)) +
      geom_point() +
      labs(y = "100 Seed Weight") +
      ylim(-.04,.3)

    # TW

    a1 <- ggplot(selection_intensity_single, aes(response_years, si_tw)) +
      geom_point() +
      labs(y = "Total Weight") +
      ylim(-.04,.3)

    # FT

    a2 <- ggplot(selection_intensity_single, aes(response_years, si_ft)) +
      geom_point() +
      labs(y = "Flowering Time") +
      ylim(-.04,.3)

    # Fit

    a3 <- ggplot(selection_intensity_single, aes(response_years, si_fit)) +
      geom_point() +
      labs(y = "Absolute Fitness") +
      ylim(-.04,.3)

    # fec

    a4 <- ggplot(selection_intensity_single, aes(response_years, si_fec)) +
      geom_point() +
      labs(y = "Fecundity") +
      ylim(-.04,.3)

    y <- ggarrange(a, a1, a2, a3, a4, top = "Selection Intensity", nrow = 1, ncol =5)
    ggsave("/bigdata/koeniglab/jmarz001/Ag-Competition/results/selection_intensity_single.png", y, width = 24, height = 12)

}

