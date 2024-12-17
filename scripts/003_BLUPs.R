#!/usr/bin/env Rscript
#SBATCH --mem=50G
#SBATCH --time=02:00:00
#SBATCH --output=/rhome/jmarz001/bigdata/Ag-Competition/scripts/003_BLUPs.stdout
#SBATCH -p koeniglab

library(tidyverse)
library(lme4)
#library(methods)
library(ggpubr)
library(inti)

setwd("/rhome/jmarz001/bigdata/Ag-Competition/data")
source("../scripts/CUSTOM_FNS.R")

# Loading Data
pheno <- read_delim("DERIVED_PHENOTYPES.tsv")

pheno <- pheno %>% 
    filter(Condition == "single") %>% 
    select(c(Genotype, Generation, Replicate, Exp_year, FT, TOTAL_MASS, GERMINATION, SEED_WEIGHT_100, FECUNDITY, FITNESS)) %>% 
    mutate(across(-c(Genotype, Generation, Replicate, Exp_year, GERMINATION), ~(scale(.) %>% as.vector))) 



#lmer(FT ~ 1 + (1|Genotype) + (1|Genotype:Exp_year), pheno_year_scale) %>% plot

#mod <-  lmer(FT ~ 1 + Exp_year + (1|Genotype), pheno_year_scale)
#pheno_year_scale$resids <- resid(mod)
#pheno_year_scale$fitted <- fitted(mod)
#ggplot(pheno_year_scale, aes(fitted, resids, color=GERMINATION)) + geom_point() + geom_smooth() + geom_hline(yintercept=0)


#lmer(TOTAL_MASS ~ 1 + Exp_year + GERMINATION + Replicate + (1|Genotype) + (1|Genotype:Exp_year), pheno_year_scale) 
#lmer(TOTAL_MASS ~ 1 + Exp_year + FITNESS + Generation + Replicate + (1|Genotype) + (1|Genotype:Exp_year), pheno_year_scale)


## not fitable
#trait = "GERMINATION"


#mod <- lmer(FECUNDITY ~ 1 + GERMINATION + FT + (1|Genotype) + (1|Genotype:Exp_year), pheno_year_scale)
#pheno_year_scale$resids <- resid(mod)
#pheno_year_scale$fits <- fitted(mod)E
#ggplot(pheno_year_scale, aes(fitted, resids)) + geom_point() + geom_smooth() + geom_hline(yintercept=0)


#mod <- lmer(FITNESS ~ 1 + FT + (1|Genotype) + (1|Genotype:Exp_year), pheno_year_scale)



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


#### testing
# what's the difference between variance of averaged phenotypes vs of extracted blups?
pheno %>% group_by(Genotype) %>% summarise(across(where(is.numeric), mean)) %>% ungroup() %>% group_by(Generation) %>% summarise(across(where(is.numeric), var))
add_generation(blup_output) %>% group_by(Generation) %>% summarise(across(where(is.numeric), var))

# number of samples per generation 
pheno_avg <- pheno %>% group_by(Genotype) %>% summarise(across(where(is.numeric), mean)) %>% ungroup() 
pheno_avg %>% group_by(Generation) %>% summarise(n())








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

## scaled data's heritability
heritability_scaled <- Breeders_funct(pheno_scaled)
write_delim(heritability_scaled, "trait_heritability_scaledphenos.tsv", "\t")




# Extracting Broad Sense Heritability per generation
# using scaled phenotypes
h0 <- pheno_scaled %>% filter(Generation == 0) %>% Breeders_funct()
h18 <- pheno_scaled %>% filter(Generation == 18) %>% Breeders_funct()
h28 <- pheno_scaled %>% filter(Generation == 28) %>% Breeders_funct()
h50 <- pheno_scaled %>% filter(Generation == 50) %>% Breeders_funct()
h58 <- pheno_scaled %>% filter(Generation == 58) %>% Breeders_funct()
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


# plots with combined factors of variance explained
herit_piv2 <- herit %>% pivot_longer(cols=c(genetic_var, total_var, H2), names_to="variance_source")
g3 <- ggplot(herit_piv2, aes(Generation, value, color=trait, shape=variance_source)) +
  geom_smooth(aes(linetype=variance_source)) +
  scale_linetype_manual(values=c("solid", "dashed", "solid")) +
  geom_point(size=3) +
  theme_bw() +
  labs(y="variance", title="Total, Genetic, and Heritable Variance over Generations") +
  facet_wrap(~trait)
ggsave("/bigdata/koeniglab/jmarz001/Ag-Competition/results/genetic_phenotypic_var_H2_over_gens_trait_facet.png", g3)








# calculate response and selection

## response w scaled phenotypes
responses_scaled <- pheno_scaled %>% group_by(Genotype) %>%
      summarise(across(where(is.numeric), \(x) mean(x, na.rm = T))) %>% 
      ungroup() %>%
      group_by(Generation) %>% 
      summarise(across(.cols = c(FT, TOTAL_MASS, GERMINATION, SEED_WEIGHT_100, FECUNDITY, FITNESS), \(x) mean(x, na.rm = T))) %>% 
      ungroup()

P_18 <-  (responses_scaled[2,] - responses_scaled[1,])/(18-0)
F18_58 <-  (responses_scaled[5,] - responses_scaled[2,])/(58-18)

resp_join <- bind_rows(P_18, F18_58) %>%
    mutate(Generation = factor(c("Parents_to_F18", "F18_to_F58"), levels = c("Parents_to_F18", "F18_to_F58"))) 
write_delim(resp_join, "/bigdata/koeniglab/jmarz001/Ag-Competition/results/response_table.tsv")


rts_df2 <- resp_join %>% pivot_longer(cols=-c(Generation), names_to='trait', values_to='response')

a2 <- ggplot(rts_df2, aes(Generation, response)) +
  geom_point() +
  facet_wrap(~trait) +
  labs(title="Scaled Response between generations", subtitle="change per generation in standard deviations") +
  theme_bw()
ggsave("/bigdata/koeniglab/jmarz001/Ag-Competition/results/response.png", a2)

a3 <- ggplot(rts_df2, aes(Generation, response, group=trait, color=trait)) +
  geom_point() +
  geom_line() +
  labs(title="Scaled Response between generations", subtitle="change per generation in standard deviations") +
  theme_bw()
ggsave("/bigdata/koeniglab/jmarz001/Ag-Competition/results/response_combined.png", a3)







## selection
herit_mini <- heritability %>% select(c(trait, H2))
herit_response <- full_join(herit_mini, rts_df2)
herit_response$selection_est <- herit_response$response / herit_response$H2

write_delim(herit_response, "trait_selection_ests.tsv", "\t")

a <- ggplot(herit_response, aes(Generation, selection_est)) +
  geom_point() +
  facet_wrap(~trait) +
  labs(y="selection estimate", x="time span") +
  theme_bw()

ggsave("/bigdata/koeniglab/jmarz001/Ag-Competition/results/selection_facet.png", a)

a2 <- ggplot(herit_response, aes(Generation, selection_est)) +  
  geom_point() +
  facet_wrap(~trait, scales="free_y") +
  labs(y="selection estimate", x="time span") +
  theme_bw()
ggsave("/bigdata/koeniglab/jmarz001/Ag-Competition/results/selection_facet_freescale.png", a2)


a3 <- ggplot(herit_response, aes(Generation, selection_est, color=trait, group=trait)) +
  geom_point() +
  geom_line() +
  labs(y="selection estimate", x="time span", title="Selection between Generations for All Traits") +
  theme_bw()

ggsave("/bigdata/koeniglab/jmarz001/Ag-Competition/results/selection.png", a3)

