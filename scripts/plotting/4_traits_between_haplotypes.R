#!/usr/bin/env Rscript

#SBATCH --ntasks=1
#SBATCH --mem=30G
#SBATCH --time=02:00:00
#SBATCH --output=/rhome/jmarz001/bigdata/Ag-Competition/scripts/plotting/4_traits_between_haplotypes.stdout
#SBATCH -p koeniglab


library(tidyverse)
library(ggpubr)
library(ggplot2)

setwd("/bigdata/koeniglab/jmarz001/Ag-Competition/data")
source("../scripts/CUSTOM_FNS.R")
df <- read_delim("weighted_generation_trait_avgs.tsv")


## plot weighted trait averages over generations

# format data longways
df_lng <- df %>% pivot_longer(cols=colnames(df)[2:7], names_to="trait", values_to="value")
df_lng$trait <- tidy_text_substitution(df_lng$trait)


g <- ggplot(df_lng, aes(x= generation, y = value)) +
      geom_point() + geom_line() +
      facet_wrap(~trait, scales="free_y") + 
      labs(title = "Trait Average over Generations", subtitle="Haplotype-Frequency Weighted", x="Generation", y="") +
      theme_bw(base_size = 16) 
ggsave("../results/trait_avg_over_gens_hap_weighted.png", g, width = 16, height = 12)




#### plot weighted averages using all available haplotypes
df <- read_delim("weighted_generation_trait_avgs_AvgdOverAllHaps.tsv")

# format data longways
df_lng <- df %>% pivot_longer(cols=colnames(df)[2:7], names_to="trait", values_to="value")
df_lng$trait <- tidy_text_substitution(df_lng$trait)


g <- ggplot(df_lng, aes(x= generation, y = value)) +
      geom_point() + geom_line() +
      facet_wrap(~trait, scales="free_y") + 
      labs(title = "Trait Average over Generations", subtitle="Haplotype-Frequency Weighted, Haplotype Values Averaged Over Generations", x="Generation", y="") +
      theme_bw(base_size = 16) 
ggsave("../results/trait_avg_over_gens_hap_weighted_AvgdOverAllHaps.png", g, width = 16, height = 12)













### Atlas

Average_Haplo_rep <- Average_Haplo_rep %>%  mutate(Atlas_Avg_Fec = 2363.51,
                  Atlas_Avg_Fitness = 21347.22,
                  Atlas_Avg_Total_Weight = 126.8267)
Average_Haplo_rep$Numbers <- ifelse(Average_Haplo_rep$Condition == "mixed", 1, 0)


### Adding Averaged Atlas values into the table and adding columns for centered data

Atlas_tbl <- sw_avg %>% filter(Genotype == "48_5") %>%
  mutate(Atlas_Avg_Fec = mean(FECUNDITY),
         Atlas_Avg_Fit = mean(FITNESS),
         Atlas_Avg_TW = mean(TOTAL_MASS))

sw_avg <- sw_avg %>% mutate(Atlas_Avg_Fec = 1413.799,
            Atlas_Avg_Fit = 12880.63,
            Atlas_Avg_TW = 64.55571)

sw_avg <- sw_avg %>% mutate(Centered_Fit = FITNESS - mean(FITNESS, na.rm = TRUE),
Centered_FT = FT - mean(FT, na.rm = TRUE),
Centered_Fec = FECUNDITY - mean(FECUNDITY, na.rm = TRUE),
Centered_TW = TOTAL_MASS - mean(TOTAL_MASS))

### Calculating Contribution of each seed to phenotype

### FECUNDITY
sw_avg$Exp_Fec_Per_Plant <- ifelse(sw_avg$Condition == "mixed",
Exp_Fec_Mixed(sw_avg$FEC),
Exp_Single(sw_avg$FEC))

### FITNESS
sw_avg$Exp_Fit_Per_Plant <- ifelse(sw_avg$Condition == "mixed",
Exp_Fit_Mixed(sw_avg$FITNESS),
Exp_Single(sw_avg$FITNESS))

### Total Weight
sw_avg$Exp_TW_Per_Plant <- ifelse(sw_avg$Condition == "mixed",
Exp_TW_mix(sw_avg$TOTAL_MASS),
Exp_Single(sw_avg$TOTAL_MASS))
