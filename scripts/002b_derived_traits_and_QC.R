#!/usr/bin/env Rscript
#SBATCH --mem=30G
#SBATCH --time=02:00:00
#SBATCH --output=/rhome/jmarz001/bigdata/Ag-Competition/scripts/002b_derived_traits_and_QC.stdout
#SBATCH -p short

## filter outliers, calcualte derived-trait values (ie. fitness), & filter extreme derived-trait values & low plant-number plots
library(tidyverse)
library(ggpubr)
library(car)
library(dunn.test)

setwd("/bigdata/koeniglab/jmarz001/Ag-Competition")
source("scripts/CUSTOM_FNS.R")


df <- read_delim("data/JOINED_PHENOTYPES.tsv")
# remove any rows without genotype
df <- df %>% 
    filter(Condition=="single") %>%
    select(-c(BED, ROW, Condition)) %>%
    select(c(Genotype, Exp_year, Replicate, Plants, FT, TOTAL_MASS, SEED_WEIGHT_100))

# calculate upper and lower bounds for each trait, in each year separately
up_bnds <- df %>% 
  group_by(Exp_year) %>% 
  summarise(across(-c(Genotype, Replicate), \(x) median(x, na.rm=T) + (2* IQR(x, na.rm=T))))

  df %>% 
  group_by(Exp_year, Replicate) %>% 
  summarise(across(-c(Genotype), \(x) median(x, na.rm=T) + (2* IQR(x, na.rm=T))))


lw_bnds <- df %>% 
  group_by(Exp_year) %>% 
  summarise(across(-c(Genotype, Replicate), \(x) median(x, na.rm=T) - (2* IQR(x, na.rm=T))))
  
df %>% 
  group_by(Exp_year, Replicate) %>% 
  summarise(across(-c(Genotype), \(x) median(x, na.rm=T) - (2* IQR(x, na.rm=T))))

# filter for outlier values

# open png to print distribution plots
pdf("results/trait_outlier_distributions.pdf")

for(year in c(2022,2023)) {
  x <- df %>% filter(Exp_year == year)

  for(trait in 4:7){

    match <- colnames(x)[trait]
    # find the upper and lower bounds for trait & year
    upp <- up_bnds %>% filter(Exp_year==year) %>% select(all_of(match))
    dwnn <- lw_bnds %>% filter(Exp_year==year) %>% select(all_of(match))

    # plot trait distribution
    g1 <- ggplot(x) + geom_histogram(aes(unlist(x[,trait]))) + 
      geom_vline(aes(xintercept=upp[[1]]), color = "#F31919") + 
      geom_vline(aes(xintercept=dwnn[[1]]), color = "#F31919") + 
      theme_bw() +
      labs(title=match, x=match, subtitle=year) 
      
      print(g1)


    # filter outlier values
    # don't remove the WHOLE ROW
    # just replace the outlier value with NA
    x[which(x[,trait] > upp[[1]]), trait] <- NA
    x[which(x[,trait] < dwnn[[1]]), trait] <- NA

    g2 <- ggplot(x) + geom_histogram(aes(unlist(x[,trait]))) + 
      geom_vline(aes(xintercept=upp[[1]]), color = "#F31919") + 
      geom_vline(aes(xintercept=dwnn[[1]]), color = "#F31919") + 
      theme_bw() +
      labs(title=match, x=match, subtitle=year) 

    print(g2)
  }

  # save the year's traits
  assign(paste0("ayear_", year), x)

}

dev.off()

# re-join data from the two years of experiment
df2 <- full_join(ayear_2022, ayear_2023)

# any rows that are entirely NA?
df2[which(rowSums(is.na(df2))>4),]







#############################################
# fecundity
## fecundity = seed produced per plot
# fec is also our proxy for fitness then
## total seed weight / seed-weight-100/100

df3 <- df2 %>% 
    mutate(MASS_PER_PLANT = TOTAL_MASS/Plants) %>% 
    mutate(SEED_COUNT = TOTAL_MASS / (SEED_WEIGHT_100/100)) %>%
    mutate(FECUNDITY = SEED_COUNT / Plants) 

write_delim(df3, "data/DERIVED_PHENOTYPES.tsv", "\t")











ggplot(df3) + geom_histogram(aes(FECUNDITY))

 df3 %>% select(c(Exp_year, Replicate, FECUNDITY)) %>% group_by(Exp_year) %>% summarise(across(-c(Replicate), \(x) median(x, na.rm=T) - (2* IQR(x, na.rm=T))))
df3 %>% select(c(Exp_year, Replicate, FECUNDITY)) %>% group_by(Exp_year) %>% summarise(across(-c(Replicate), \(x) median(x, na.rm=T) + (2* IQR(x, na.rm=T))))

#problem w fecundity....
#values over 600...
df3[which(df3$FECUNDITY > 400),] 









































####################
# check for large differences between replicates
# specifically for traits with low replicate-correlation values: 100-seed weight & total mass

#### total seed mass / yield
# calculate total mass relative to plant # and THEN checking range of value differences
df_tm <- df2 %>% select(c(Genotype, Exp_year, Replicate, Plants, TOTAL_MASS)) %>% 
  mutate("mass_per_plant"=TOTAL_MASS/Plants) %>% select(-c(Plants, TOTAL_MASS)) 


up_bnds <- df_tm %>% 
  group_by(Exp_year) %>% 
  summarise(across(-c(Genotype, Replicate), \(x) median(x, na.rm=T) + (2* IQR(x, na.rm=T))))
lw_bnds <- df_tm %>% 
  group_by(Exp_year) %>% 
  summarise(across(-c(Genotype, Replicate), \(x) median(x, na.rm=T) - (2* IQR(x, na.rm=T))))


for(year in c(2022,2023)) {
  x <- df_tm %>% filter(Exp_year == year)
  trait <- 4
  match <- colnames(x)[trait]

  # find the upper and lower bounds for trait & year
  upp <- up_bnds %>% filter(Exp_year==year) %>% select(all_of(match))
  dwnn <- lw_bnds %>% filter(Exp_year==year) %>% select(all_of(match))

  # plot trait distribution
  g1 <- ggplot(x) + geom_histogram(aes(unlist(x[,trait]))) + 
    geom_vline(aes(xintercept=upp[[1]]), color = "#F31919") + 
    geom_vline(aes(xintercept=dwnn[[1]]), color = "#F31919") + 
    theme_bw() +
    labs(title=match, x=match, subtitle=year) 
    
    print(g1)


  # filter outlier values
  # just replace the outlier value with NA
  x[which(x[,trait] > upp[[1]]), trait] <- NA
  x[which(x[,trait] < dwnn[[1]]), trait] <- NA

  g2 <- ggplot(x) + geom_histogram(aes(unlist(x[,trait]))) + 
    geom_vline(aes(xintercept=upp[[1]]), color = "#F31919") + 
    geom_vline(aes(xintercept=dwnn[[1]]), color = "#F31919") + 
    theme_bw() +
    labs(title=match, x=match, subtitle=year) 

  print(g2)

  # save the year's traits
  assign(paste0("ayear_", year), x)
  }


# re-join data from the two years of experiment
join_tm <- full_join(ayear_2022, ayear_2023)
  
join_tm %>% 
  pivot_wider(names_from=Replicate, names_prefix = "REP_", values_from='mass_per_plant') %>% group_by(Exp_year) %>% summarise('correlation'=cor(REP_1, REP_2, use="pairwise.complete.obs"))




## filter large differences
df_tm2$diff <- df_tm2$REP_1 - df_tm2$REP_2
ggplot(df_tm2) +
  geom_histogram(aes(abs(diff))) +
  facet_wrap(~Exp_year, nrow=2)
# is more than 15s of grams of difference in mass per plant reasonable?
# I'm raising it from 10g, and lowering it from 20-23g
median(df_tm2$diff, na.rm=T) + (2* IQR(df_tm2$diff, na.rm=T))
# calculated outlier limit is ~14, so I'll set it to 15g difference which also feels reflective of 
up_bnds <- df_tm2 %>%
  group_by(Exp_year) %>%
  summarise(across(diff, \(x) median(x, na.rm=T) + (2* IQR(x, na.rm=T))))

  upp <- up_bnds %>% filter(Exp_year==year) %>% select(all_of(match))
  dwnn <- lw_bnds %>% filter(Exp_year==year) %>% select(all_of(match))



df_tm[which(df_tm$Exp_year==2022 & abs(df_tm$diff) > 15),]
df_tm[which(df_tm$Exp_year==2023 & abs(df_tm$diff) > 15),] 
# and that'll only remove 1 & 4 values (respective to year)




ttt <- round((nrow(df_tm[which(abs(df_tm$diff) > 10),])/nrow(df_tm))*100)

print("removing samples with more than 10g seed-yield per-plant difference between replicates")
print(paste0("removes ", (ttt), " percent of remaining samples"))


df_tm_filt <- df_tm %>% filter(abs(diff) < 10)
ggplot(df_tm_filt) +
  geom_histogram(aes(diff)) +
  facet_wrap(~Exp_year, nrow=2)
df_tm_filt %>% group_by(Exp_year) %>% summarise('correlation'=cor(REP_1, REP_2, use="pairwise.complete.obs"))

## how do you combine this filter...back into the rest of the data frame
# if you filtered the huge differences, it was a combination of genotype & year, not replicate number
kept <- df_tm_filt %>% select(c(Genotype, Exp_year))
# now...filter joined...
smth <- right_join(df2, kept, by=c("Genotype", "Exp_year"))
# that's a list of...the genotype-year combos to keep...
# but for the other genotype-year combos in df2... you just need to remove the total mass measurements for it

df2[which(Genotype %in% geno_list & Exp_year == 2022), which(colnames(df2=="TOTAL_MASS"))] <- NA





#### 100-seed mass
df_sw <- df2 %>% select(c(Genotype, Exp_year, Replicate, SEED_WEIGHT_100))

df_sw <- df_sw %>% 
  pivot_wider(names_from=Replicate, names_prefix = "REP_", values_from='SEED_WEIGHT_100')
df_sw %>% group_by(Exp_year) %>% summarise('correlation'=cor(REP_1, REP_2, use="pairwise.complete.obs"))

# calculate difference between R1 and R2
df_sw$rep_diff <- df_sw$REP_1 - df_sw$REP_2

# plot the distribution of those differences
ggplot(df_sw) +
  geom_histogram(aes(abs(rep_diff))) +
  facet_wrap(~Exp_year, nrow=2)

# more than 2 g difference between replicates ...is that reasonable?
# 2 g per 100 seeds...
# I don't think I'll filter that any further


