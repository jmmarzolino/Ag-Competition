#!/usr/bin/env Rscript
#SBATCH --mem=30G
#SBATCH --time=02:00:00
#SBATCH --output=/rhome/jmarz001/bigdata/Ag-Competition/scripts/002_derived_traits_and_QC.stdout
#SBATCH -p short

library(tidyverse)
library(ggpubr)
library(car)
library(dunn.test)

setwd("/bigdata/koeniglab/jmarz001/Ag-Competition")
source("scripts/CUSTOM_FNS.R")


df <- read_delim("data/JOINED_PHENOTYPES.tsv")
# remove any rows without genotype
df <- df %>% 
    filter(!is.na(Genotype)) %>%
    select(-c(BED, ROW))

# fn adds generation col based on 'Genotype' col
df <- add_generation(df)


# filter for outlier values
## 100 seed weight
upper <- median(df$SEED_WEIGHT_100, na.rm=T) + (2 * IQR(df$SEED_WEIGHT_100, na.rm=T))
lower <- median(df$SEED_WEIGHT_100, na.rm=T) - (2 * IQR(df$SEED_WEIGHT_100, na.rm=T))

g <- ggplot(df, aes(SEED_WEIGHT_100)) + geom_histogram() + 
  geom_vline(aes(xintercept=upper), color = "#F31919") + 
  geom_vline(aes(xintercept=lower), color = "#F31919") + 
  geom_vline(aes(xintercept=median(SEED_WEIGHT_100, na.rm=T)), color = "#F31919", linetype="dashed") + 
  theme_bw()

ggsave("results/seed_weight_outlier_distribution.png", g)
# record outlier 100 seed weight
print('SEED_WEIGHT_100')
summary(df$SEED_WEIGHT_100)
df[which(df$SEED_WEIGHT_100 > upper | df$SEED_WEIGHT_100 < lower),]
# remove extreme 100 seed weight values
df <- df[which(df$SEED_WEIGHT_100 <= round(upper)),]
df <- df[which(df$SEED_WEIGHT_100 >= round(lower)),]
summary(df$SEED_WEIGHT_100)


## total seed weight
upper <- median(df$TOTAL_MASS, na.rm=T) + (2 * IQR(df$TOTAL_MASS, na.rm=T))
lower <- median(df$TOTAL_MASS, na.rm=T) - (2 * IQR(df$TOTAL_MASS, na.rm=T))

g <- ggplot(df, aes(TOTAL_MASS)) + geom_histogram() + 
  geom_vline(aes(xintercept=upper), color = "#F31919") + 
  geom_vline(aes(xintercept=lower), color = "#F31919") + 
  geom_vline(aes(xintercept=median(TOTAL_MASS, na.rm=T)), color = "#F31919", linetype="dashed") + 
  theme_bw()

ggsave("results/total_weight_outlier_distribution.png", g)
# record outlier total weight
print('TOTAL_MASS')
summary(df$TOTAL_MASS)
df[which(df$TOTAL_MASS > upper | df$TOTAL_MASS < lower),]

df <- df[which(df$TOTAL_MASS < upper),]
df <- df[which(df$TOTAL_MASS > lower),]
summary(df$TOTAL_MASS)


## flowering time
upper <- median(df$FT, na.rm=T) + (2 * IQR(df$FT, na.rm=T))
lower <- median(df$FT, na.rm=T) - (3 * IQR(df$FT, na.rm=T))
# manually increasing the lower FT limit to keep early flowering lines

g <- ggplot(df, aes(FT)) + geom_histogram() + 
  geom_vline(aes(xintercept=upper), color = "#F31919") + 
  geom_vline(aes(xintercept=lower), color = "#F31919") + 
  geom_vline(aes(xintercept=median(FT, na.rm=T)), color = "#F31919", linetype="dashed") + 
  theme_bw() +
  facet_wrap(~Exp_year)


ggsave("results/flowering_time_outlier_distribution.png", g)
# record outlier total weight
print('FT')
summary(df$FT)
df[which(df$FT > upper | df$FT < lower),]

df <- df[which(df$FT < upper),]
df <- df[which(df$FT >= lower),]
summary(df$FT)




# calculate derived phenotypes
## germination = (plants/10) 
over_sown <- df %>% filter(Plants > 10)
# for plots have 11 or 12 seeds planted, 
# assume max number planted is the same
over_sown$GERMINATION <- 1
# plants / 11 or 12 will be 100% germination

df2 <- df %>% filter(Plants <= 10)
df2 <- df2 %>% mutate(GERMINATION = Plants / 10)

df3 <- full_join(df2, over_sown) 
#, by=c('Genotype', 'Condition', 'Exp_year', 'Plants', 'FT', 'TOTAL_MASS', 'SEED_WEIGHT_100', 'Generation', 'GERMINATION'))



# fecundity
## fecundity = seed produced per plot
## total seed weight / seed-weight-100/100
## fitness = germination * fecundity

df3 <- df3 %>% 
    mutate(SEED_COUNT = TOTAL_MASS / (SEED_WEIGHT_100/100)) %>%
    mutate(FECUNDITY = SEED_COUNT / Plants) %>% 
    select(-c(Plants)) %>% 
    mutate(FITNESS = GERMINATION * FECUNDITY) 



#problem w fecundity!
#values over 500 up to 2500
df3[which(df3$FECUNDITY > 1000),] 
# 20 rows w extreme fecundity, all w germination of 0.1-0.2

## for the MOST part, rows w outlier fecundity have germination between 1-5 out of 10 seeds (0.1 - 0.5)
# one occurance of high germination 0.9
df3[which(df3$FECUNDITY > 500),] -> x
summary(x$GERMINATION)
summary(df3$GERMINATION)
# so fecundity cutoff of 500 is probably too low on its own

# either ... filter outlier fecundity values that have germination below 50%...
# or filter all low germination values
# or filter all outlier fecundity values...
# the 0.9 germination row does have a huge total weight, like double any others in the set

# how many rows removed if we filter out low germination?
df3 %>% filter(GERMINATION < 0.3 )
(xy <- df3 %>% filter(GERMINATION < 0.3 ) %>% nrow)
xy / nrow(df3)
# removes about 2% of plots

(yy <- df3 %>% filter(GERMINATION <= 0.3 ) %>% nrow)
yy / nrow(df3)
# only 3% fitlered if lower suvival threshhold to 0.3 

#upper <- median(df3$FECUNDITY) + 2*IQR(df3$FECUNDITY)
#df3 <- df3 %>% filter(FECUNDITY <= upper)

df3 <- df3 %>% filter(GERMINATION > 0.2 & FECUNDITY < 750) 
write_delim(df3, "data/DERIVED_PHENOTYPES.tsv", "\t")
