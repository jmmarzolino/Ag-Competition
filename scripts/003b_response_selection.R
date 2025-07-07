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

# calculate response and selection


# Loading Data
#df <- fread("DERIVED_PHENOTYPES.tsv")
#df <- df %>% 
#        select(-c(Replicate, SEED_COUNT)) %>%
#        group_by(Genotype, Exp_year) %>% 
#        summarise(across(where(is.numeric), mean)) %>%
#        ungroup() %>% select(-Exp_year) %>% 
#        group_by(Genotype) %>% summarise(across(where(is.numeric), mean)) %>%
#        ungroup() 

df <- fread("trait_BLUPs.tsv")
#df <- full_join(df, dfb, by=c('Genotype'), suffix=c("", "_blup"))
df <- add_generation(df)
#df <- df %>% select(c('Genotype', 'Generation', 'FT', 'TOTAL_MASS', 'SEED_WEIGHT_100', 'Germination', 'FECUNDITY', 'MASS_PER_PLANT', 'FT_blup', 'TOTAL_MASS_blup', 'Germination_blup', 'SEED_WEIGHT_100_blup', 'FECUNDITY_blup', 'MASS_PER_PLANT_blup'))



# RESPONSE
## CALCULATE RESPONSE IN SD AND NORMAL UNITS
response <- df %>%
      group_by(Generation) %>% 
      summarise(across(where(is.numeric), \(x) mean(x, na.rm = T))) %>% 
      ungroup()

response_scaled <- df %>%
      mutate(across(-c(Genotype, Generation), ~(scale(.) %>% as.vector))) %>%
      group_by(Generation) %>% 
      summarise(across(where(is.numeric), \(x) mean(x, na.rm = T))) %>% 
      ungroup()

P_18 <-  (response[2,] - response[1,])/(18-0)
F18_58 <-  (response[5,] - response[2,])/(58-18)

P_18_sd <-  (response_scaled[2,] - response_scaled[1,])/(18-0)
F18_58_sd <-  (response_scaled[5,] - response_scaled[2,])/(58-18)

resp_join <- bind_rows(P_18, F18_58) %>%
    mutate(Generation = factor(c("Parents_to_F18", "F18_to_F58"), levels = c("Parents_to_F18", "F18_to_F58"))) 
resp_join_sd <- bind_rows(P_18_sd, F18_58_sd) %>%
    mutate(Generation = factor(c("Parents_to_F18_sd", "F18_to_F58_sd"), levels = c("Parents_to_F18_sd", "F18_to_F58_sd"))) 

responses_joined <- full_join(resp_join, resp_join_sd)

write_delim(responses_joined, "/bigdata/koeniglab/jmarz001/Ag-Competition/results/response_table.tsv")


rts_df2 <- resp_join_sd %>% pivot_longer(cols=-c(Generation), names_to='trait', values_to='response')

# tidy trait name text before plotting
rts_df2$trait <- tidy_text_substitution(rts_df2$trait)
# and tidy generation text
rts_df2$Generation <- gsub("Parents_to_F18_sd", "Parents to F18", rts_df2$Generation)
rts_df2$Generation <- gsub("F18_to_F58_sd", "F18 to F58", rts_df2$Generation)

# re-order generation indicators 
rts_df2$Generation <- factor(rts_df2$Generation, levels = c("Parents to F18", "F18 to F58"))


a2 <- ggplot(rts_df2, aes(Generation, response)) +
  geom_point() +
  geom_hline(aes(yintercept=0), color="grey") +
  facet_wrap(~trait) +
  labs(title="Scaled Response between Generations", subtitle="average change per generation in standard deviations", x="time span") +
  theme_bw(base_size = 14) 

ggsave("/bigdata/koeniglab/jmarz001/Ag-Competition/results/response.png", a2, height=7*2, width=7*3, units="in")

a3 <- ggplot(rts_df2, aes(Generation, response, group=trait, color=trait)) +
  geom_point() +
  geom_line() +
  labs(title="Scaled Response between Generations", subtitle="average change per generation in standard deviations", x="time span") +
  theme_bw(base_size = 14)
ggsave("/bigdata/koeniglab/jmarz001/Ag-Competition/results/response_combined.png", a3, height=7, width=7, units="in")




## selection
herit <- fread("trait_heritability.tsv")

# write out one table with response & selection calculated for trait units & stand. dev.s
rts <- responses_joined %>% pivot_longer(cols=-c(Generation), names_to='trait', values_to='response')
herit_response <- full_join(herit, rts)

herit_response$selection_est <- herit_response$response / herit_response$H2
write_delim(herit_response, "trait_selection_ests.tsv", "\t")

# plots comparing response & selection estimates with sd only
rts <- resp_join_sd %>% pivot_longer(cols=-c(Generation), names_to='trait', values_to='response')
herit_response <- full_join(herit, rts)
herit_response$selection_est <- herit_response$response / herit_response$H2

# tidy trait name text before plotting
herit_response$trait <- tidy_text_substitution(herit_response$trait)
# and tidy generation text
herit_response$Generation <- gsub("Parents_to_F18_sd", "Parents to F18", herit_response$Generation)
herit_response$Generation <- gsub("F18_to_F58_sd", "F18 to F58", herit_response$Generation)
# re-order generation indicators 
herit_response$Generation <- factor(herit_response$Generation, levels = c("Parents to F18", "F18 to F58"))


a <- ggplot(herit_response, aes(Generation, selection_est)) +
  geom_point() +
  geom_hline(aes(yintercept=0), color="grey") +
  facet_wrap(~trait) +
  labs(y="selection estimate", x="time span") +
  theme_bw(base_size = 14) 

ggsave("/bigdata/koeniglab/jmarz001/Ag-Competition/results/selection_facet.png", a, height=7*2, width=7*3, units="in")

a3 <- ggplot(herit_response, aes(Generation, selection_est, color=trait, group=trait)) +
  geom_point() +
  geom_line() +
  labs(y="selection estimate", x="time span", title="Selection between Generations") +
  theme_bw(base_size = 14) 

ggsave("/bigdata/koeniglab/jmarz001/Ag-Competition/results/selection.png", a3, height=7, width=7, units="in")

