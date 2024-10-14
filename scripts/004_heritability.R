#!/usr/bin/env Rscript
#SBATCH --mem=50G
#SBATCH --time=02:00:00
#SBATCH --output=/rhome/jmarz001/bigdata/Ag-Competition/scripts/003_change_over_gens.stdout
#SBATCH -p short

library(lme4)
library(methods)
library(gridExtra)

setwd("/rhome/jmarz001/bigdata/Ag-Competition/data")
source("../scripts/CUSTOM_FNS.R")

# Loading Data
PHENO_FULL <- read_delim("FT_FITNESS.tsv")

# Creating Single and Mixed dataframes for BLUP analysis
tmp <- PHENO_FULL %>% select(-c("Haplotype", "TOTAL_SEED_COUNT", "REL_FITNESS", "Plants", "SURVIVAL"))
tmp <- relocate(tmp, Genotype, Generation, Condition, Replicate, Exp_year, FT, TOTAL_MASS, SEED_WEIGHT_100, FECUNDITY, ABS_FITNESS)

# Standardize data

standardize <- function(x){
  (x - mean(x, na.rm = T))/sd(x, na.rm = T)
}

tmp$FT <- standardize(tmp$FT)
tmp$TOTAL_MASS <- standardize(tmp$TOTAL_MASS)
tmp$FEC <- tmp$FEC/sd(tmp$FEC, na.rm = T)
tmp$ABS_FITNESS <- tmp$ABS_FITNESS/sd(tmp$ABS_FITNESS, na.rm = T)
tmp$SEED_WEIGHT_100 <- standardize(tmp$SEED_WEIGHT_100)

# Breeders Function for calculating Heritability

# P

Breeders_funct_0 <- function(x) {
  var_tab <- tibble(trait= character(), genetic_var=double(), total_var=double())

  for(i in c(6:ncol(x))){
    df <- x %>%
      filter(Generation == 0) %>%
      select(c(all_of(i), Genotype, Replicate, Generation, Exp_year))
    # calculate mixed model
    genotype_mod <- lmer(df[,1][[1]] ~ Exp_year + (1| Genotype), data = df)
    ## summary of the model
    geno_mod_summ <- summary(genotype_mod)
    data.frame(geno_mod_summ$varcor)
    ## genotypic variance
    genotypic_variance <- (data.frame(geno_mod_summ$varcor)$sdcor[1]^2)
    ## total variance
    total_variance <- genotypic_variance + (data.frame(geno_mod_summ$varcor)$sdcor[2]^2)
    new_row <- tibble(trait= colnames(PHENO_SINGLE[,i]), genetic_var=genotypic_variance, total_var=total_variance)
    var_tab <- bind_rows(var_tab, new_row)
  }
  ## ratio of Vg and Vp
  var_tab$H2 <- var_tab$genetic_var / var_tab$total_var
  return(var_tab)
}



# F18

Breeders_funct_18 <- function(x) {
  var_tab <- tibble(trait= character(), genetic_var=double(), total_var=double())

  for(i in c(6:ncol(x))){
    df <- x %>%
      filter(Generation == 18) %>%
      select(c(all_of(i), Genotype, Replicate, Generation, Exp_year))
    # calculate mixed model
    genotype_mod <- lmer(df[,1][[1]] ~ Exp_year + (1| Genotype), data = df)
    ## summary of the model
    geno_mod_summ <- summary(genotype_mod)
    data.frame(geno_mod_summ$varcor)
    ## genotypic variance
    genotypic_variance <- (data.frame(geno_mod_summ$varcor)$sdcor[1]^2)
    ## total variance
    total_variance <- genotypic_variance + (data.frame(geno_mod_summ$varcor)$sdcor[2]^2)
    new_row <- tibble(trait= colnames(PHENO_SINGLE[,i]), genetic_var=genotypic_variance, total_var=total_variance)
    var_tab <- bind_rows(var_tab, new_row)
  }
  ## ratio of Vg and Vp
  var_tab$H2 <- var_tab$genetic_var / var_tab$total_var
  return(var_tab)
}

#F28

Breeders_funct_28 <- function(x) {
  var_tab <- tibble(trait= character(), genetic_var=double(), total_var=double())

  for(i in c(6:ncol(x))){
    df <- x %>%
      filter(Generation == 28) %>%
      select(c(all_of(i), Genotype, Replicate, Generation, Exp_year))
    # calculate mixed model
    genotype_mod <- lmer(df[,1][[1]] ~ Exp_year + (1| Genotype), data = df)
    ## summary of the model
    geno_mod_summ <- summary(genotype_mod)
    data.frame(geno_mod_summ$varcor)
    ## genotypic variance
    genotypic_variance <- (data.frame(geno_mod_summ$varcor)$sdcor[1]^2)
    ## total variance
    total_variance <- genotypic_variance + (data.frame(geno_mod_summ$varcor)$sdcor[2]^2)
    new_row <- tibble(trait= colnames(PHENO_SINGLE[,i]), genetic_var=genotypic_variance, total_var=total_variance)
    var_tab <- bind_rows(var_tab, new_row)
  }
  ## ratio of Vg and Vp
  var_tab$H2 <- var_tab$genetic_var / var_tab$total_var
  return(var_tab)
}

#F50

Breeders_funct_50 <- function(x) {
  var_tab <- tibble(trait= character(), genetic_var=double(), total_var=double())

  for(i in c(6:ncol(x))){
    df <- x %>%
      filter(Generation == 50) %>%
      select(c(all_of(i), Genotype, Replicate, Generation, Exp_year))
    # calculate mixed model
    genotype_mod <- lmer(df[,1][[1]] ~ Exp_year + (1| Genotype), data = df)
    ## summary of the model
    geno_mod_summ <- summary(genotype_mod)
    data.frame(geno_mod_summ$varcor)
    ## genotypic variance
    genotypic_variance <- (data.frame(geno_mod_summ$varcor)$sdcor[1]^2)
    ## total variance
    total_variance <- genotypic_variance + (data.frame(geno_mod_summ$varcor)$sdcor[2]^2)
    new_row <- tibble(trait= colnames(PHENO_SINGLE[,i]), genetic_var=genotypic_variance, total_var=total_variance)
    var_tab <- bind_rows(var_tab, new_row)
  }
  ## ratio of Vg and Vp
  var_tab$H2 <- var_tab$genetic_var / var_tab$total_var
  return(var_tab)
}

#F58

Breeders_funct_58 <- function(x) {
  var_tab <- tibble(trait= character(), genetic_var=double(), total_var=double())

  for(i in c(6:ncol(x))){
    df <- x %>%
      filter(Generation == 58) %>%
      select(c(all_of(i), Genotype, Replicate, Generation, Exp_year))
    # calculate mixed model
    genotype_mod <- lmer(df[,1][[1]] ~ Exp_year + (1| Genotype), data = df)
    ## summary of the model
    geno_mod_summ <- summary(genotype_mod)
    data.frame(geno_mod_summ$varcor)
    ## genotypic variance
    genotypic_variance <- (data.frame(geno_mod_summ$varcor)$sdcor[1]^2)
    ## total variance
    total_variance <- genotypic_variance + (data.frame(geno_mod_summ$varcor)$sdcor[2]^2)
    new_row <- tibble(trait= colnames(PHENO_SINGLE[,i]), genetic_var=genotypic_variance, total_var=total_variance)
    var_tab <- bind_rows(var_tab, new_row)
  }
  ## ratio of Vg and Vp
  var_tab$H2 <- var_tab$genetic_var / var_tab$total_var
  return(var_tab)
}


### SINGLE

PHENO_SINGLE <- tmp %>% filter(Condition == "single")

# Creating BLUP dataframe

blup_output <- data.frame(matrix(vector(), 255, 1, dimnames = list(c(), c("Genotype"))))
blup_output <- unique(PHENO_SINGLE[, 1])
blup_output <- add_generation(blup_output)


for (i in 6:ncol(PHENO_SINGLE)){
  model <- lmer(as.numeric(unlist(PHENO_SINGLE[,i])) ~ Exp_year + (1|Genotype), data = PHENO_SINGLE)
  summary(model)
  model_ranefs <- ranef(model)
  geno_blups <- model_ranefs[[1]][[1]]
  geno_blups <- tibble(geno_blups)
  geno_blups$Replicate <- rownames(model_ranefs[[1]])
  colnames(geno_blups) <- c(paste(colnames(PHENO_SINGLE[,i]),"blup",sep="_"), "Genotype")
  assign(x = paste("blup_df",i,sep="_"), value = geno_blups)
}

blup_data_single <-  full_join(blup_output, blup_df_6)
blup_data_single <- full_join(blup_data_single, blup_df_7)
blup_data_single <- full_join(blup_data_single, blup_df_8)
blup_data_single <- full_join(blup_data_single, blup_df_9)
blup_data_single <- full_join(blup_data_single, blup_df_10)

# Extracting Broad Sense Heritability (SINGLE)

# P

h0 <- Breeders_funct_0(PHENO_SINGLE)
h0 <- h0 %>% mutate(Generation = 0)
h0 <- h0 %>% mutate(H2 = round(h0$H2, 2),
                    genetic_var = round(h0$genetic_var, 2),
                    total_var = round(h0$total_var, 2))

# F18

h18 <- Breeders_funct_18(PHENO_SINGLE)
h18 <- h18 %>% mutate(Generation = 18)
h18 <- h18 %>% mutate(H2 = round(h18$H2, 2),
                    genetic_var = round(h18$genetic_var, 2),
                    total_var = round(h18$total_var, 2))

# F28

h28 <- Breeders_funct_28(PHENO_SINGLE)
h28 <- h28 %>% mutate(Generation = 28)
h28 <- h28 %>% mutate(H2 = round(h28$H2, 2),
                    genetic_var = round(h28$genetic_var, 2),
                    total_var = round(h28$total_var, 2))

# F50

h50 <- Breeders_funct_50(PHENO_SINGLE)
h50 <- h50 %>% mutate(Generation = 50)
h50 <- h50 %>% mutate(H2 = round(h50$H2, 2),
                    genetic_var = round(h50$genetic_var, 2),
                    total_var = round(h50$total_var, 2))

# F58

h58 <- Breeders_funct_58(PHENO_SINGLE)
h58 <- h58 %>% mutate(Generation = 58)
h58 <- h58 %>% mutate(H2 = round(h58$H2, 2),
                    genetic_var = round(h58$genetic_var, 2),
                    total_var = round(h58$total_var, 2))

heritability_data_single <- rbind(h0, h18, h28, h50, h58)

# Heritability over time (Single)

ggplot(heritability_data_single, aes(Generation, H2, color = trait)) +
  geom_point() +
  geom_smooth(alpha = .4) +
  labs(y = "Broad Sense Heritability",
       title = "Generational Change in Broad Sense Heritability (Single)")
ggsave("/bigdata/koeniglab/jmarz001/Ag-Competition/results/broad_sense_over_time.png", width = 14, height = 10)

# Genetic Variance over time (Single)

a <- ggplot(heritability_data_single, aes(Generation, genetic_var, color = trait)) +
  geom_point() +
  geom_smooth() +
  labs(y = "Genetic Variance",
       title = "Generational Change in Genetic Variance (Single)") +
  ylim(0, 2.5)

# Phenotypic Variance over time (Single)

a1 <- ggplot(heritability_data_single, aes(Generation, total_var, color = trait)) +
  geom_point() +
  geom_smooth() +
  ylim(0, 2.5) +
  labs(y = "Phenotypic Variance",
       title = "Generational Change in Phenotypic Variance")

y <- arrangeGrob(a, a1, ncol =2, nrow = 1)
ggsave("/bigdata/koeniglab/jmarz001/Ag-Competition/results/scatterplot_genetic_phenotypic_var_single.png", y, width = 22, height = 12)

heritability_data_single <- heritability_data_single %>% group_by(trait) %>% summarise(across(where(is.numeric), .fns = mean, na.rm = T)) %>% select(-c("Generation")) %>% ungroup()

# Response to Selection (Single)

subtract <- function(x, y){
  x[2] - x[1]
}

# P to F58

preserved <- PHENO_SINGLE %>% filter(Generation == 0 | Generation == 58)
P_F58 <- preserved %>% group_by(Generation) %>% summarise(across(.cols = c(SEED_WEIGHT_100, FECUNDITY, ABS_FITNESS, TOTAL_MASS, FT), mean, na.rm = T)) %>% ungroup()
P_F58 <- P_F58 %>% summarise(across(.cols = c(SEED_WEIGHT_100, FECUNDITY, ABS_FITNESS, TOTAL_MASS, FT), subtract)) %>% ungroup()
P_F58 <- P_F58 %>% mutate(SEED_WEIGHT_100 = SEED_WEIGHT_100/58,
                          FECUNDITY = FECUNDITY/58,
                          ABS_FITNESS = ABS_FITNESS/58,
                          FT = FT/58,
                          TOTAL_MASS = TOTAL_MASS/58,
                          response_years = "P_F58")
# P to F18

preserved <- PHENO_SINGLE %>% filter(Generation == 0 | Generation == 18)
P_F18 <- preserved %>% group_by(Generation) %>% summarise(across(.cols = c(SEED_WEIGHT_100, FECUNDITY, ABS_FITNESS, TOTAL_MASS, FT), mean, na.rm = T)) %>% ungroup()
P_F18 <- P_F18 %>% summarise(across(.cols = c(SEED_WEIGHT_100, FECUNDITY, ABS_FITNESS, TOTAL_MASS, FT), subtract)) %>% ungroup()
P_F18 <- P_F18 %>% mutate(SEED_WEIGHT_100 = SEED_WEIGHT_100/18,
                          FECUNDITY = FECUNDITY/18,
                          ABS_FITNESS = ABS_FITNESS/18,
                          FT = FT/18,
                          TOTAL_MASS = TOTAL_MASS/18,
                          response_years = "P_F18")

# F18 to F58

preserved <- PHENO_SINGLE %>% filter(Generation == 18 | Generation == 58)
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

g <- arrangeGrob(a1, a2, a3, a4, a5, top = "Response to Selection (Single)", nrow = 1, ncol = 5)
ggsave("/bigdata/koeniglab/jmarz001/Ag-Competition/results/Response_to_Selection_Single.png", g, width = 22, height = 10)

# Calculating Selection Intensity (single)

rts_df <- rts_df %>% mutate(H2_FT = heritability_data_single$H2[4],
                            H2_100 = heritability_data_single$H2[1],
                            H2_tw = heritability_data_single$H2[5],
                            H2_fec = heritability_data_single$H2[3],
                            H2_fit = heritability_data_single$H2[2])

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

y <- arrangeGrob(a, a1, a2, a3, a4, top = "Selection Intensity (Single)", nrow = 1, ncol =5)
ggsave("/bigdata/koeniglab/jmarz001/Ag-Competition/results/selection_intensity_single.png", y, width = 24, height = 12)

### MIXED

PHENO_MIXED <- tmp %>% filter(Condition == "mixed")

for (i in 6:ncol(PHENO_MIXED)){
  model <- lmer(as.numeric(unlist(PHENO_MIXED[,i])) ~ Exp_year + (1|Genotype), data = PHENO_MIXED)
  summary(model)
  model_ranefs <- ranef(model)
  geno_blups <- model_ranefs[[1]][[1]]
  geno_blups <- tibble(geno_blups)
  geno_blups$Replicate <- rownames(model_ranefs[[1]])
  colnames(geno_blups) <- c(paste(colnames(PHENO_MIXED[,i]),"blup_df",sep="_"), "Genotype")
  assign(x = paste("blup_mixed_df",i,sep="_"), value = geno_blups)
}

# Creating blup dataframe for mixed

blup_output <- data.frame(matrix(vector(), 255, 1, dimnames = list(c(), c("Genotype"))))
blup_output <- unique(PHENO_MIXED[, 1])

blup_output <- add_generation(blup_output)


# Joining blups for all traits into a mixed blup dataframe

blup_data_mixed <- full_join(blup_output, blup_mixed_df_6)
blup_data_mixed <- full_join(blup_data_mixed, blup_df_7)
blup_data_mixed <- full_join(blup_data_mixed, blup_df_8)
blup_data_mixed <- full_join(blup_data_mixed, blup_df_9)
blup_data_mixed <- full_join(blup_data_mixed, blup_df_10)

# Extracting Broad Sense Heritability (MIXED)

# P

h0 <- Breeders_funct_0(PHENO_MIXED)
h0 <- h0 %>% mutate(Generation = 0)
h0 <- h0 %>% mutate(H2 = round(h0$H2, 2),
                      genetic_var = round(h0$genetic_var, 2),
                      total_var = round(h0$total_var, 2))

# F18

h18 <- Breeders_funct_18(PHENO_MIXED)
h18 <- h18 %>% mutate(Generation = 18)
h18 <- h18 %>% mutate(H2 = round(h18$H2, 2),
                      genetic_var = round(h18$genetic_var, 2),
                      total_var = round(h18$total_var, 2))

# F28

h28 <- Breeders_funct_28(PHENO_MIXED)
h28 <- h28 %>% mutate(Generation = 28)
h28 <- h28 %>% mutate(H2 = round(h28$H2, 2),
                      genetic_var = round(h28$genetic_var, 2),
                      total_var = round(h28$total_var, 2))
# F50

h50 <- Breeders_funct_50(PHENO_MIXED)
h50 <- h50 %>% mutate(Generation = 50)
h50 <- h50 %>% mutate(H2 = round(h50$H2, 2),
                      genetic_var = round(h50$genetic_var, 2),
                      total_var = round(h50$total_var, 2))

# F58

h58 <- Breeders_funct_58(PHENO_MIXED)
h58 <- h58 %>% mutate(Generation = 58)
h58 <- h58 %>% mutate(H2 = round(h58$H2, 2),
                      genetic_var = round(h58$genetic_var, 2),
                      total_var = round(h58$total_var, 2))

heritability_data_mixed <- rbind(h0, h18, h28, h50, h58)

# Broad Sense Heritability Over Time (Mixed)

ggplot(heritability_data_mixed, aes(Generation, H2, color = trait)) +
  geom_point() +
  geom_smooth() +
  labs(y = "Broad Sense Heritability",
       title = "Generational Change in Broad Sense Heritability (Mixed)")
ggsave("/bigdata/koeniglab/jmarz001/Ag-Competition/data/broad_sense_mixed.png", width = 14, height = 10)

# Genetic Variance Over Time (Mixed)

a <- ggplot(heritability_data_mixed, aes(Generation, genetic_var, color = trait)) +
  geom_point() +
  geom_smooth() +
  labs(y = "Genetic Variance",
       title = "Generational Change in Genetic Variance (Mixed)") +
  ylim(0, 1.2)

# Phenotypic Variance Over Time (Mixed)

a1 <- ggplot(heritability_data_mixed, aes(Generation, total_var, color = trait)) +
  geom_point() +
  geom_smooth() +
  labs(y = "Phenotypic Variance",
       title = "Generational Change in Phenotypic Variance (Mixed)") +
  ylim(0, 1.2)

y <- arrangeGrob(a, a1, ncol = 2, nrow =1)
ggsave("/bigdata/koeniglab/jmarz001/Ag-Competition/results/scatterplots_genetic_phenotypic_var_mixed.png", y, width = 20, height = 12)

heritability_data_mixed <- heritability_data_mixed %>% group_by(trait) %>% summarise(across(where(is.numeric), .fns = mean, na.rm = T)) %>% select(-c("Generation")) %>% ungroup()

# Response to Selection (MIXED)

# P to F58

preserved <- PHENO_MIXED %>% filter(Generation == 0 | Generation == 58)
P_F58 <- preserved %>% group_by(Generation) %>% summarise(across(.cols = c(SEED_WEIGHT_100, FECUNDITY, ABS_FITNESS, TOTAL_MASS, FT), mean, na.rm = T)) %>% ungroup()
P_F58 <- P_F58 %>% summarise(across(.cols = c(SEED_WEIGHT_100, FECUNDITY, ABS_FITNESS, TOTAL_MASS, FT), subtract)) %>% ungroup()
P_F58 <- P_F58 %>% mutate(SEED_WEIGHT_100 = SEED_WEIGHT_100/58,
                          FECUNDITY = FECUNDITY/58,
                          ABS_FITNESS = ABS_FITNESS/58,
                          FT = FT/58,
                          TOTAL_MASS = TOTAL_MASS/58,
                          response_years = "P_F58")
# P to F18

preserved <- PHENO_MIXED %>% filter(Generation == 0 | Generation == 18)
P_F18 <- preserved %>% group_by(Generation) %>% summarise(across(.cols = c(SEED_WEIGHT_100, FECUNDITY, ABS_FITNESS, TOTAL_MASS, FT), mean, na.rm = T)) %>% ungroup()
P_F18 <- P_F18 %>% summarise(across(.cols = c(SEED_WEIGHT_100, FECUNDITY, ABS_FITNESS, TOTAL_MASS, FT), subtract)) %>% ungroup()
P_F18 <- P_F18 %>% mutate(SEED_WEIGHT_100 = SEED_WEIGHT_100/18,
                          FECUNDITY = FECUNDITY/18,
                          ABS_FITNESS = ABS_FITNESS/18,
                          FT = FT/18,
                          TOTAL_MASS = TOTAL_MASS/18,
                          response_years = "P_F18")

# F18 to F58

preserved <- PHENO_MIXED %>% filter(Generation == 18 | Generation == 58)
F18_F58 <- preserved %>% group_by(Generation) %>% summarise(across(.cols = c(SEED_WEIGHT_100, FECUNDITY, ABS_FITNESS, TOTAL_MASS, FT), mean, na.rm = T)) %>% ungroup()
F18_F58 <- F18_F58 %>% summarise(across(.cols = c(SEED_WEIGHT_100, FECUNDITY, ABS_FITNESS, TOTAL_MASS, FT), subtract)) %>% ungroup()
F18_F58 <- F18_F58 %>% mutate(SEED_WEIGHT_100 = SEED_WEIGHT_100/40,
                              FECUNDITY = FECUNDITY/40,
                              ABS_FITNESS = ABS_FITNESS/40,
                              FT = FT/40,
                              TOTAL_MASS = TOTAL_MASS/40,
                              response_years = "F18_F58")

# Response to Selection (Mixed)

rts_df_mix <- rbind(P_F18, P_F58, F18_F58)

a1 <- ggplot(rts_df_mix, aes(response_years, TOTAL_MASS)) +
  geom_point() +
  ylim(0, .0075)

a2 <- ggplot(rts_df_mix, aes(response_years, SEED_WEIGHT_100)) +
  geom_point() +
  ylim(0, .0075)

a3 <- ggplot(rts_df_mix, aes(response_years, FT)) +
  geom_point() +
  ylim(0, .0075)

a4 <- ggplot(rts_df_mix, aes(response_years, FECUNDITY)) +
  geom_point() +
  ylim(0, .0075)

a5 <- ggplot(rts_df_mix, aes(response_years, ABS_FITNESS)) +
  geom_point() +
  ylim(0, .0075)

y <- arrangeGrob(a1, a2, a3, a4,a5, top = "Response to Selection (Mixed)", nrow = 1, ncol = 5)
ggsave("/bigdata/koeniglab/jmarz001/Ag-Competition/results/response_to_selection_mix.png", y, width = 24, height = 14)


# Calculating Selection Intensity (Mixed)

rts_df_mix <- rts_df_mix %>% mutate(H2_FT = heritability_data_mixed$H2[4],
                            H2_100 = heritability_data_mixed$H2[1],
                            H2_tw = heritability_data_mixed$H2[5],
                            H2_fec = heritability_data_mixed$H2[3],
                            H2_fit = heritability_data_mixed$H2[2])

rts_df_mix <- rts_df_mix %>% mutate(si_100 = SEED_WEIGHT_100/ H2_100,
                            si_fec = FECUNDITY / H2_fec,
                            si_ft = FT/H2_FT,
                            si_tw = TOTAL_MASS/H2_tw,
                            si_fit = ABS_FITNESS/H2_fit)

selection_intensity_mix <- rts_df_mix %>% select(starts_with("si"), response_years)

# 100 SW SI

a <- ggplot(selection_intensity_mix, aes(response_years, si_100)) +
  geom_point() +
  labs(y = "100 Seed Weight") +
  ylim(0,2.8)

# TW

a1 <- ggplot(selection_intensity_mix, aes(response_years, si_tw)) +
  geom_point() +
  labs(y = "Total Weight") +
  ylim(0,2.8)

# FT

a2 <- ggplot(selection_intensity_mix, aes(response_years, si_ft)) +
  geom_point() +
  labs(y = "Flowering Time") +
  ylim(0,2.8)

# Fit

a3 <- ggplot(selection_intensity_mix, aes(response_years, si_fit)) +
  geom_point() +
  labs(y = "Absolute Fitness") +
  ylim(0,2.8)

# fec

a4 <- ggplot(selection_intensity_mix, aes(response_years, si_fec)) +
  geom_point() +
  labs(y = "Fecundity") +
  ylim(0,2.8)

y <- arrangeGrob(a, a1, a2, a3, a4, top = "Selection Intensity (Mixed)", nrow = 1, ncol =5)
ggsave("/bigdata/koeniglab/jmarz001/Ag-Competition/results/selection_intensity_mix.png", y, width = 24, height = 12)
