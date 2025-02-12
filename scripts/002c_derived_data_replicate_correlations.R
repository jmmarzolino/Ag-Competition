#!/usr/bin/env Rscript
#SBATCH --mem=50G
#SBATCH --time=02:00:00
#SBATCH --output=/rhome/jmarz001/bigdata/Ag-Competition/scripts/002c_derived_data_replicate_correlations.stdout
#SBATCH -p short

library(tidyverse)
#library(ggpubr)
library(data.table)
#library(ggrepel)


setwd("/rhome/jmarz001/bigdata/Ag-Competition/data")
source("../scripts/CUSTOM_FNS.R")
# load data
df <- fread("DERIVED_PHENOTYPES.tsv")


replicate_df <- df %>%
  pivot_longer(cols=c(FT, TOTAL_MASS, SEED_WEIGHT_100, Plants, SEED_COUNT, FECUNDITY, MASS_PER_PLANT), names_to='PHENOTYPE', values_to="VALUE") %>%
  pivot_wider(names_from=Replicate, names_prefix="REP_", values_from='VALUE')

# another way of summarizing trait replicate correlations by year, condition, and phenotype
replicate_df %>% group_by(Exp_year, PHENOTYPE) %>% summarise('correlation'=cor(REP_1, REP_2, use="pairwise.complete.obs", method="pearson"))

# across years
replicate_df %>% group_by(PHENOTYPE) %>% summarise('correlation'=cor(REP_1, REP_2, use="pairwise.complete.obs"))
# correlation across years results in some weird results - particularly in total weight, and seed weight
# and with expected but near 0 results for 'Plants'
# 'Plants' shouldn't necessarily have a correlation so that's fine


# check out what's causing the miniscule correlations for total-mass (which is the source of low correlations in fec and fit as well)

replicate_df %>% filter(PHENOTYPE=="TOTAL_MASS") %>%
ggplot(aes(x=REP_1, y=REP_2))+
  geom_point(aes(color=as.factor(Generation), shape=as.factor(Exp_year), group=as.factor(Generation)), size=2) +
  geom_smooth(method='lm') +
  labs(title="total mass") +
  theme_bw()

# 
replicate_df %>% filter(PHENOTYPE=="TOTAL_MASS") %>%
ggplot(aes(x=REP_1, y=REP_2))+
  geom_point() +
  geom_smooth(method="lm") +
  facet_wrap(Exp_year~Generation, nrow=2)


replicate_df %>%
ggplot(aes(x=REP_1, y=REP_2))+
  geom_point() +
  geom_smooth(method="lm") +
  facet_wrap(Exp_year~PHENOTYPE, scales="free")
ggsave("derived_trait_replicate_correlation_plots.png")





## how to filter...
## based on a verylow correlation of total mass between replicates within-year
# calculate difference between R1 and R2
replicate_df <- replicate_df %>% filter(PHENOTYPE=="TOTAL_MASS"|PHENOTYPE=="SEED_WEIGHT_100")
replicate_df$rep_diff <- replicate_df$REP_1 - replicate_df$REP_2

#replicate_df %>% group_by(Exp_year, PHENOTYPE) %>% summarise(mean(rep_diff, na.rm=T))
#replicate_df %>% group_by(Exp_year, PHENOTYPE) %>% reframe(range(rep_diff, na.rm=T))

# plot the distribution of those differences
## total mass
replicate_df %>% filter(PHENOTYPE=="TOTAL_MASS") %>%
ggplot() +
  geom_histogram(aes(rep_diff)) +
  facet_wrap(~Exp_year, nrow=2)
# vals beyond +/-100 seems like a good cutoff point for total mass

# filter the largest differences
tm_filt <- replicate_df %>% filter(PHENOTYPE=="TOTAL_MASS") %>% filter(rep_diff<100 & rep_diff>-100)
problem_reps_tm <- replicate_df %>% filter(PHENOTYPE=="TOTAL_MASS") %>% filter(rep_diff>100 | rep_diff< -100)


## 100-seed weight
replicate_df %>% filter(PHENOTYPE=="SEED_WEIGHT_100") %>%
ggplot() +
  geom_histogram(aes(rep_diff)) +
  facet_wrap(~Exp_year, nrow=2)
# vals beyond +/-2 seems like a good cutoff point for total mass

# filter the largest differences
sw_filt <- replicate_df %>% filter(PHENOTYPE=="SEED_WEIGHT_100") %>% filter(rep_diff<2 & rep_diff>-2)
problem_reps_sw <- replicate_df %>% filter(PHENOTYPE=="SEED_WEIGHT_100") %>% filter(rep_diff>2 | rep_diff< -2)




# check correlation of filtered values
tm_filt %>% group_by(Exp_year) %>% summarise('correlation'=cor(REP_1, REP_2, use="pairwise.complete.obs"))

sw_filt %>% group_by(Exp_year) %>% summarise('correlation'=cor(REP_1, REP_2, use="pairwise.complete.obs"))
# that makes a big difference...
# let's try proceeding with that filter to start

## effect across years
tm_filt %>% summarise('correlation'=cor(REP_1, REP_2, use="pairwise.complete.obs"))

sw_filt %>% summarise('correlation'=cor(REP_1, REP_2, use="pairwise.complete.obs"))


