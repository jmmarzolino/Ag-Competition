#!/usr/bin/env Rscript
#SBATCH --mem=20G
#SBATCH --time=01:00:00
#SBATCH --output=/rhome/jmarz001/bigdata/Ag-Competition/scripts/haplotype-resampled-population/004b_variance_test_hapreppop.stdout
#SBATCH -p koeniglab

setwd("/rhome/jmarz001/bigdata/Ag-Competition/data")
source("../scripts/CUSTOM_FNS.R")
library(tidyverse)
library(data.table)
library(car)

# load data
df <- fread("trait_BLUPs_HapRepPop.tsv") %>% tibble
parents <- fread("BLUPs.tsv")
parents <- add_generation(parents)
parents <- parents %>% filter(Generation == 0) %>% select(-Genotype)
df <- bind_rows(df, parents)

# calculate variance per generation
trait_var <- df %>% group_by(Generation) %>% 
              summarise(across(where(is.numeric), \(x) var(x, na.rm=T))) %>% 
              pivot_longer(cols = -Generation, names_to = 'traits') %>% 
              pivot_wider(names_from = Generation, values_from = value, names_prefix="Gen") 
fwrite(trait_var, "haprep_gen_vars.tsv")

trait_var <- trait_var %>% 
              select(c('traits', 'Gen0', 'Gen18', 'Gen58'))

# difference between generational variance
trait_var$P_18_diff <- trait_var$Gen18 - trait_var$Gen0
trait_var$F18_58_diff <- trait_var$Gen58 - trait_var$Gen18
trait_var$P_18_diff_direction <- as.character(sign(trait_var$P_18_diff))
trait_var$F18_58_diff_direction <- as.character(sign(trait_var$F18_58_diff))

direction_cols <- grep("direction", colnames(trait_var))
for(d in direction_cols) {
  trait_var[which(trait_var[,d] > 0), d] <- "increase"
  trait_var[which(trait_var[,d] < 0), d] <- "decrease"
  trait_var[which(trait_var[,d] == 0), d] <- "no change"
}

print("direction of change from P to F18")
table(trait_var$P_18_diff_direction)
print("direction of change from F18 to F58")
table(trait_var$F18_58_diff_direction)

# what traits inc in var?
print("traits w inc in var")
trait_var %>% filter(P_18_diff_direction == "increase" | F18_58_diff_direction == "increase")

df$Generation <- as.factor(df$Generation)
p_values <- tibble(trait=NULL, P_18_pval=NULL, F18_F58_pval=NULL)

for(j in 1:(ncol(df)-1)) {
  trait <- (colnames(df[,j]))

  #comparing variance of parent gen to F18, and F18 to F58
  P_18_df <- df %>% filter(Generation == 0 | Generation == 18)
  F18_F58_df <- df %>% filter(Generation == 18 | Generation == 58)

  #### calculate and save p-values
  x <- (leveneTest(data=P_18_df, unlist(P_18_df[, j]) ~ Generation))[1,3]
  y <- (leveneTest(data=F18_F58_df, unlist(F18_F58_df[, j]) ~ Generation))[1,3]

  ## store the results of Levene tests
  new_row <- tibble(trait=trait, P_18_pval=x, F18_F58_pval=y)
  p_values <- bind_rows(data=p_values, new_row)
}

# add a significance col for each set of generations
p_values$P_18_sig <- (p_values$P_18_pval < 0.05)
p_values$F18_F58_sig <- (p_values$F18_F58_pval < 0.05)

# combine variance change w p-values
variance_table <- full_join(trait_var, p_values, by=c("traits"="trait"))

write_delim(variance_table, "variance_test_pvals.tsv", "\t")