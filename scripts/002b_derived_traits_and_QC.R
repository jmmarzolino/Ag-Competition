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
    select(c(Genotype, Exp_year, Replicate, Germination, FT, TOTAL_MASS, SEED_WEIGHT_100))

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


# how many values were removed from each trait with the outlier filtering?
print("percent missing data for traits, after outlier filtering")
df22 <- df2 %>% filter(Exp_year==2022)
df23 <- df2 %>% filter(Exp_year==2023)

print('number of rows per year')
nrow(df22)
nrow(df23)

print("amount of missing data for the 3 measured traits")
df2 %>% group_by(Exp_year) %>% reframe(across(c(Germination, FT, TOTAL_MASS, SEED_WEIGHT_100), \(x) sum(is.na(x))))

print("percent of missing data per year")
df22 %>% summarise(across(c(Germination, FT, TOTAL_MASS, SEED_WEIGHT_100), \(x) (sum(is.na(x))/508)*100))

df23 %>% summarise(across(c(Germination, FT, TOTAL_MASS, SEED_WEIGHT_100), \(x) (sum(is.na(x))/510)*100))



#############################################
# fecundity
## fecundity = seed produced per plot
# fec is also our proxy for fitness then
## total seed weight / seed-weight-100/100

df3 <- df2 %>% 
    mutate(MASS_PER_PLANT = TOTAL_MASS/Germination) %>% 
    mutate(SEED_COUNT = TOTAL_MASS / (SEED_WEIGHT_100/100)) %>%
    mutate(FECUNDITY = SEED_COUNT / Germination) 

write_delim(df3, "data/DERIVED_PHENOTYPES.tsv", "\t")


# how many values exist from each derived trait?
print("percent missing data for derived traits")
df22 <- df3 %>% filter(Exp_year==2022)
df23 <- df3 %>% filter(Exp_year==2023)

print('number of rows per year')
nrow(df22)
nrow(df23)

print("amount of missing data for traits")
df3 %>% group_by(Exp_year) %>% reframe(across(where(is.numeric), \(x) sum(is.na(x))))

print("percent of missing data per year")
df22 %>% summarise(across(where(is.numeric), \(x) (sum(is.na(x))/508)*100))

df23 %>% summarise(across(where(is.numeric), \(x) (sum(is.na(x))/508)*100))
