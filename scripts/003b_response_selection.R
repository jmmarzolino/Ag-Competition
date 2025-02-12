#!/usr/bin/env Rscript
#SBATCH --mem=50G
#SBATCH --time=02:00:00
#SBATCH --output=/rhome/jmarz001/bigdata/Ag-Competition/scripts/003b_response_selection.stdout
#SBATCH -p koeniglab

library(tidyverse)
library(inti)
library(data.table)

setwd("/rhome/jmarz001/bigdata/Ag-Competition/data")
source("../scripts/CUSTOM_FNS.R")

# Loading Data
df <- fread("DERIVED_PHENOTYPES.tsv")
df <- df %>% 
        filter(Condition == "single") %>% select(-c(Condition, Replicate, SEED_COUNT)) %>%
        group_by(Genotype, Exp_year) %>% summarise(across(where(is.numeric), mean)) %>%
        ungroup() %>% select(-Exp_year) %>% 
        group_by(Genotype) %>% summarise(across(where(is.numeric), mean)) %>%
        ungroup() 

dfb <- fread("trait_BLUPs.tsv")
df <- full_join(df, dfb, by=c('Genotype'), suffix=c("", "_blup"))
df <- add_generation(df)
df <- df %>% select(c('Generation', 'FT', 'TOTAL_MASS', 'SEED_WEIGHT_100', 'Plants', 'FECUNDITY', 'MASS_PER_PLANT', 'FT_blup', 'TOTAL_MASS_blup', 'Plants_blup', 'SEED_WEIGHT_100_blup', 'FECUNDITY_blup', 'MASS_PER_PLANT_blup'))



# calculate response and selection

## response w scaled phenotypes
responses_scaled <- pheno_scaled %>% group_by(Genotype) %>%
      summarise(across(where(is.numeric), \(x) mean(x, na.rm = T))) %>% 
      ungroup() %>%
      group_by(Generation) %>% 
      summarise(across(.cols = c(FT, TOTAL_MASS, Plants, SEED_WEIGHT_100, FECUNDITY, MASS_PER_PLANT), \(x) mean(x, na.rm = T))) %>% 
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

