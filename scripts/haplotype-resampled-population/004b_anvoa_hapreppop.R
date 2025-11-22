#!/usr/bin/env Rscript
#SBATCH --mem=32G
#SBATCH --time=01:00:00
#SBATCH --output=/rhome/jmarz001/bigdata/Ag-Competition/scripts/haplotype-resampled-population/004b_anova_hapreppop.stdout
#SBATCH -p koeniglab

setwd("/rhome/jmarz001/bigdata/Ag-Competition/data")
library(tidyverse)
library(data.table)

### statistically test for trait changes between generations
# with phenotypes sampled to more accurately represent their proportion in the orignal population
hap <- fread("trait_BLUPs_HapRepPop.tsv")

hap2 <- hap %>% mutate(across(-c(Generation), ~(scale(.) %>% as.vector))) %>% select(-Generation)
colnames(hap2) <- gsub("(.*)", "\\1_scaled", colnames(hap2))
hap <- bind_cols(hap, hap2)


####  GENERATION AVERAGES    
# Make table of generational means & size of generation differences 
gen_avgs <- hap %>% select(-contains("FECUNDITY")) %>% group_by(Generation) %>% summarise(across(where(is.numeric), \(x) mean(x, na.rm=T))) 
gen_avgs <- gen_avgs %>% select(c("Generation", "FT", "FT_scaled", "SEED_WEIGHT_100", "SEED_WEIGHT_100_scaled"))
# calculate relative fecundity to add to avg table
# (seed#/mean seed # across all BLUPs including parents)
#fec <- hap %>% select(c("Generation", contains("FECUNDITY"))) 
#fec %>% mutate()

#%>% group_by(Generation) %>% summarise(across(where(is.numeric), \(x) mean(x, na.rm=T))) 


fwrite(gen_avgs, "haprep_gen_avgs.tsv")


###   make list of traits w sig diff var ~ gens
# remove traits w no variance & select relevant cols
variance <- read_delim("variance_test_pvals.tsv") %>% filter(!Gen0 == 0) %>% select(c(traits, ends_with("sig")))

print('traits w sig var diff, by intergen period')
variance %>% filter(P_18_sig)
variance %>% filter(F18_F58_sig)

# divide traits into two groups
# those w equal variance across gens and those without
# check to see if there are any generations w equal var
if(sum(rowSums(variance[,2:3]) < length(variance[,2]))){
  equal_var <- variance[which(rowSums(variance[,2:3]) == 0), ]
  unequal_var <- variance[-(which(variance$traits %in% equal_var$traits)), ]
} else {
  unequal_var <- variance
}

###    DIFFERENCES BETWEEN GENERATIONS
# calculate each generation's trait average
gen_means <- hap %>%
  group_by(Generation) %>%
  summarise(across(where(is.numeric), \(x) mean(x, na.rm=T), .names="{.col}"))

# reshape the data so each generation mean is a column
gen_means2 <- gen_means %>%
           pivot_longer(cols = -Generation, names_to = 'traits') %>%
           pivot_wider(names_from = Generation, values_from = value)
# renaming columns for clarity
colnames(gen_means2) <- c("trait", "Gen0_mean", "Gen18_mean", "Gen28_mean", "Gen50_mean", "Gen58_mean")


# find difference btwn gen avgs
gen_means2$P_18_effect <- gen_means2$Gen18_mean - gen_means2$Gen0_mean
gen_means2$F18_F58_effect <- gen_means2$Gen58_mean - gen_means2$Gen18_mean

# add column indicating direction of inter-generational change
gen_means2$P_18_direction <- as.character(sign(gen_means2$P_18_effect))
gen_means2$F18_F58_direction <- as.character(sign(gen_means2$F18_F58_effect))

direction_cols <- grep("direction", colnames(gen_means2))

for(d in direction_cols) {
  gen_means2[which(gen_means2[,d] > 0), d] <- "increase"
  gen_means2[which(gen_means2[,d] < 0), d] <- "decrease"
  gen_means2[which(gen_means2[,d] == 0), d] <- "no change"
}


####    test for sig change over gens 
##  calculate p-values for each hap & bind into table
anova_pval <- function(x){
    z <- x %>%
        select(-Generation) %>%
        map(~ summary(aov(.x ~ as.factor(x$Generation))))

  p_vals <- NA
  length(p_vals) <- length(z)

  for(i in 1:length(z)){
      p_vals[i] <- unlist(z[i])[9]
      traits <- names(z)
      p_values <- data.frame(p_vals, traits) 
    }
  return(p_values)
}

##  parametric test
# filter to trait and generation col
if(length(equal_var$traits) >0){
  equal_var_traits <- hap %>% select(c(Generation, all_of(equal_var$traits)))
  p_values <- anova_pval(equal_var_traits)
  p_values$sig <- (p_values$p_vals < 0.05)
  p_values <- p_values[order(p_values$p_vals, decreasing = F), ]
  print('ANOVA results')
  print(summary(p_values))
  print(p_values)

  p_values$equal_var <- T
}


##  non-parametric test
unequal_trait_df <- hap %>% select(c(Generation, all_of(unequal_var$traits)))

kruskal_fits <- tibble('trait' = colnames(unequal_trait_df)[2:ncol(unequal_trait_df)], 'kruskal_pval'=NA)
for(i in 2:ncol(unequal_trait_df)){
        print(colnames(unequal_trait_df[,..i]))
        print(result <- kruskal.test(unlist(unequal_trait_df[,..i]) ~ Generation,
                    data = unequal_trait_df))

        kruskal_fits[i-1, 2] <- result$p.value
}

kruskal_fits$kruskal_sig <- kruskal_fits$kruskal_pval < 0.05


##  combine the data sets
kruskal_fits$equal_var <- F

if(exists(quote(p_values))){
  x <- full_join(p_values, kruskal_fits, by=c('p_vals'='kruskal_pval', 'traits'='trait', 'sig'='kruskal_sig', 'equal_var'))
  colnames(x) <- c("anova_pval", "trait", "unequal_means", "equal_vars")
} else {
  x <- kruskal_fits %>% select(c("kruskal_pval", "trait", "kruskal_sig", "equal_var"))
  colnames(x) <- c("anova_pval", "trait", "unequal_means", "equal_vars")
}

write_delim(x, "../results/anova.tsv", "\t")
