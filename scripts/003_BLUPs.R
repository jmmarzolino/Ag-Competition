#!/usr/bin/env Rscript
#SBATCH --mem=50G
#SBATCH --time=02:00:00
#SBATCH --output=/rhome/jmarz001/bigdata/Ag-Competition/scripts/003_BLUPs.stdout
#SBATCH -p koeniglab

library(tidyverse)
library(lme4)
#library(methods)
library(ggpubr)

setwd("/rhome/jmarz001/bigdata/Ag-Competition/data")
source("../scripts/CUSTOM_FNS.R")

# Loading Data
pheno <- read_delim("DERIVED_PHENOTYPES.tsv")
# select columns and scale phenotype data
pheno <- pheno %>% filter(Condition == "single") %>% select(c(Genotype, Generation, Condition, Replicate, Exp_year, FT, TOTAL_MASS, GERMINATION, SEED_WEIGHT_100, FECUNDITY, FITNESS))
pheno_scaled <- pheno %>% mutate(across(-c(Genotype, Generation, Condition, Replicate, Exp_year), ~(scale(.) %>% as.vector))) 

# check model fits for different traits 
#for(i in 6:ncol(pheno)) {
#  print(i)
#  print(colnames(pheno[,i]))
#  print(AIC(lmer(as.numeric(unlist(pheno[,i])) ~ Exp_year + Replicate + Generation + (1|Genotype), data=pheno)))
#  print(AIC(lmer(as.numeric(unlist(pheno[,i])) ~ Exp_year + Replicate + (1|Genotype), data=pheno)))
#  print(AIC(lmer(as.numeric(unlist(pheno[,i])) ~ Exp_year + (1|Exp_year:Replicate) + (1|Genotype), data=pheno)))
#  print(AIC(lmer(as.numeric(unlist(pheno[,i])) ~ Exp_year + Generation + (1|Genotype), data=pheno)))
#  print(AIC(lmer(as.numeric(unlist(pheno[,i])) ~ Exp_year + (1|Generation:Genotype), data=pheno)))
#  print(AIC(lmer(as.numeric(unlist(pheno[,i])) ~ Exp_year + (1+Generation|Genotype), data=pheno)))
#  print(AIC(lmer(as.numeric(unlist(pheno[,i])) ~ Exp_year + (1|Genotype), data=pheno)))
#}
# Exp_year + (1|Genotype) is a winner across the board


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

Breeders_funct <- function(pheno){
    # mixed model loop to calculate heritability

    # set up empty df for storing results
    var_tab <- tibble(trait= character(), genetic_var=double(), total_var=double())

    # calculate for each trait
    for (i in 6:ncol(pheno)){
    df <- pheno %>%
        select(c(all_of(i), Genotype, Exp_year))
      # calculate mixed model
      genotype_mod <- lmer(df[,1][[1]] ~ Exp_year + (1| Genotype), data = df)
      ## summary of the model
      geno_mod_summ <- summary(genotype_mod)

      ## genotypic variance
      genotypic_variance <- (data.frame(geno_mod_summ$varcor)$sdcor[1]^2)
      ## total variance
      total_variance <- genotypic_variance + (data.frame(geno_mod_summ$varcor)$sdcor[2]^2)
      new_row <- tibble(trait= colnames(pheno[,i]), genetic_var=genotypic_variance, total_var=total_variance)
      var_tab <- bind_rows(var_tab, new_row)
    }

    ## heritability = ratio of Vg and Vp
    var_tab$H2 <- var_tab$genetic_var / var_tab$total_var
    return(var_tab)
}

heritability <- Breeders_funct(pheno)
write_delim(heritability, "trait_heritability.tsv", "\t")


# what's the difference between variance of averaged phenotypes vs of extracted blups?
pheno %>% group_by(Genotype) %>% summarise(across(where(is.numeric), mean)) %>% ungroup() %>% group_by(Generation) %>% summarise(across(where(is.numeric), var))
add_generation(blup_output) %>% group_by(Generation) %>% summarise(across(where(is.numeric), var))

# number of samples per generation 
pheno_avg <- pheno %>% group_by(Genotype) %>% summarise(across(where(is.numeric), mean)) %>% ungroup() 
pheno_avg %>% group_by(Generation) %>% summarise(n())



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
g <- ggplot(herit, aes(Generation, H2, color = trait)) +
  geom_point() +
  geom_smooth(alpha = .4) +
  theme_bw() +
  labs(y = "Broad Sense Heritability",
      title = "Change in Broad Sense Heritability over Generations")
ggsave("/bigdata/koeniglab/jmarz001/Ag-Competition/results/H2_over_gens.png", g)

# Genetic Variance over time
a <- ggplot(herit, aes(Generation, genetic_var, color = trait)) +
  geom_point() +
  geom_smooth() +
  theme_bw() +
  labs(y = "Genetic Variance",
      title = "Change in Genetic Variance over Generations")
ggsave("/bigdata/koeniglab/jmarz001/Ag-Competition/results/genetic_var_over_gens.png", a)

# Phenotypic Variance over time
a1 <- ggplot(herit, aes(Generation, total_var, color = trait)) +
  geom_point() +
  geom_smooth() +
  theme_bw() +
  labs(y = "Phenotypic Variance",
      title = "Change in Phenotypic Variance over Generations")
ggsave("/bigdata/koeniglab/jmarz001/Ag-Competition/results/phenotypic_var_over_gens.png", a1)
#herit2 <- herit %>% group_by(trait) %>% summarise(across(where(is.numeric), .fns = mean, na.rm = T)) %>% select(-c("Generation")) %>% ungroup()



## scaled data's heritability
heritability_scaled <- Breeders_funct(pheno_scaled)

# a few plots with combined factors of variance explained
# using scaled phenotypes
h0 <- pheno_scaled %>% filter(Generation == 0) %>% Breeders_funct()
h18 <- pheno_scaled %>% filter(Generation == 18) %>% Breeders_funct()
h28 <- pheno_scaled %>% filter(Generation == 28) %>% Breeders_funct()
h50 <- pheno_scaled %>% filter(Generation == 50) %>% Breeders_funct()
h58 <- pheno_scaled %>% filter(Generation == 58) %>% Breeders_funct()
herit <- rbind(h0, h18, h28, h50, h58)
herit$Generation <- c(rep(0, 6), rep(18, 6), rep(28, 6), rep(50, 6), rep(58,6))


herit_piv <- herit %>% pivot_longer(cols=c(genetic_var, total_var), names_to="variance_source")
g1 <- ggplot(herit_piv, aes(Generation, value, color=trait, shape=variance_source)) +
  geom_smooth() +
  geom_point(size=3) +
  theme_bw() +
  labs(y="variance", title="Total and Genetic Variance over Generations")
ggsave("/bigdata/koeniglab/jmarz001/Ag-Competition/results/genetic_phenotypic_var_over_gens_scaled.png", g1)

g2 <- ggplot(herit_piv, aes(Generation, value, color=trait, shape=variance_source)) +
  geom_smooth() +
  geom_point(size=3) +
  theme_bw() +
  labs(y="variance", title="Total and Genetic Variance over Generations") +
  facet_wrap(~trait)
ggsave("/bigdata/koeniglab/jmarz001/Ag-Competition/results/genetic_phenotypic_var_over_gens_trait_facet_scaled.png", g2)

herit_piv2 <- herit %>% pivot_longer(cols=c(genetic_var, total_var, H2), names_to="variance_source")
g3 <- ggplot(herit_piv2, aes(Generation, value, color=trait, shape=variance_source)) +
  geom_smooth(aes(linetype=variance_source)) +
  scale_linetype_manual(values=c("solid", "dashed", "solid")) +
  geom_point(size=3) +
  theme_bw() +
  labs(y="variance", title="Total, Genetic, and Heritable Variance over Generations") +
  facet_wrap(~trait)
ggsave("/bigdata/koeniglab/jmarz001/Ag-Competition/results/genetic_phenotypic_var_H2_over_gens_trait_facet_scaled.png", g3)

gsumm <- ggarrange(g1, g2, g3)
ggsave("/bigdata/koeniglab/jmarz001/Ag-Competition/results/Vg_Vp_H2_over_gens_all_traits_scaled.png", gsumm)




# calculate response and selection

## response
responses <- pheno %>% 
      group_by(Genotype) %>%
      summarise(across(where(is.numeric), \(x) mean(x, na.rm = T))) %>% 
      ungroup() %>%
      group_by(Generation) %>% 
      summarise(across(.cols = c(FT, TOTAL_MASS, GERMINATION, SEED_WEIGHT_100, FECUNDITY, FITNESS), \(x) mean(x, na.rm = T))) %>% 
      ungroup()

P_18 <-  (responses[2,] - responses[1,])/(18-0)
F18_58 <-  (responses[5,] - responses[2,])/(58-18)

resp_join <- bind_rows(P_18, F18_58) %>%
    mutate(Generation = factor(c("Parents_to_F18", "F18_to_F58"), levels = c("Parents_to_F18", "F18_to_F58"))) 
rts_df <- resp_join %>% pivot_longer(cols=-c(Generation), names_to='trait', values_to='response')

a1 <- ggplot(rts_df, aes(Generation, response)) +
  geom_point() +
  facet_wrap(~trait, scales="free") +
  labs(title="Response between generations") +
  theme_bw()
ggsave("/bigdata/koeniglab/jmarz001/Ag-Competition/results/response.png", a1)


## response w scaled phenotypes
responses_scaled <- pheno_scaled %>% group_by(Genotype) %>%
      summarise(across(where(is.numeric), \(x) mean(x, na.rm = T))) %>% 
      ungroup() %>%
      group_by(Generation) %>% 
      summarise(across(.cols = c(FT, TOTAL_MASS, GERMINATION, SEED_WEIGHT_100, FECUNDITY, FITNESS), \(x) mean(x, na.rm = T))) %>% 
      ungroup()

P_18 <-  (responses_scaled[2,] - responses_scaled[1,])/(18-0)
F18_58 <-  (responses_scaled[5,] - responses_scaled[2,])/(58-18)

resp_join2 <- bind_rows(P_18, F18_58) %>%
    mutate(Generation = factor(c("Parents_to_F18", "F18_to_F58"), levels = c("Parents_to_F18", "F18_to_F58"))) 
rts_df2 <- resp_join2 %>% pivot_longer(cols=-c(Generation), names_to='trait', values_to='response')

a2 <- ggplot(rts_df2, aes(Generation, response)) +
  geom_point() +
  facet_wrap(~trait) +
  labs(title="Scaled Response between generations", subtitle="change per generation in standard deviations") +
  theme_bw()
ggsave("/bigdata/koeniglab/jmarz001/Ag-Competition/results/response_scaled.png", a2)


a3 <- ggplot(rts_df2, aes(Generation, response, group=trait, color=trait)) +
  geom_point() +
  #facet_wrap(~trait) +
  labs(title="Scaled Response between generations", subtitle="change per generation in standard deviations") +
  theme_bw()
ggsave("/bigdata/koeniglab/jmarz001/Ag-Competition/results/response_scaled_combined.png", a3)


## join response tables w normal & standard deviation units
rts_join <- full_join(rts_df, rts_df2, by=c("Generation", "trait"), suffix=c("", "_scaled"))
rts_join <- rts_join[order(rts_join$trait),]
write_delim(rts_join, "/bigdata/koeniglab/jmarz001/Ag-Competition/results/response_table_traitunits_and_scaled.tsv")




## selection
herit_mini <- heritability %>% select(c(trait, H2))
herit_response <- full_join(herit_mini, rts_join)
herit_response$selection_est <- herit_response$response / herit_response$H2
herit_response$selection_est_scaled <- herit_response$response_scaled / herit_response$H2

write_delim(herit_response, "trait_selection_ests.tsv", "\t")

a <- ggplot(herit_response, aes(Generation, selection_est)) +
  geom_point() +
  facet_wrap(~trait) +
  labs(y="selection estimate", x="time span") +
  theme_bw()

ggsave("/bigdata/koeniglab/jmarz001/Ag-Competition/results/selection.png", a)

a2 <- ggplot(herit_response, aes(Generation, selection_est)) +
  geom_point() +
  facet_wrap(~trait, scales="free_y") +
  labs(y="selection estimate", x="time span") +
  theme_bw()
ggsave("/bigdata/koeniglab/jmarz001/Ag-Competition/results/selection_freescale.png", a2)


a3 <- ggplot(herit_response, aes(Generation, selection_est_scaled, color=trait)) +
  geom_point(position = position_dodge(0.2)) +
  labs(y="selection estimate", x="time span", title="Selection between Generations for All Traits") +
  theme_bw()

ggsave("/bigdata/koeniglab/jmarz001/Ag-Competition/results/selection.png", a3)

