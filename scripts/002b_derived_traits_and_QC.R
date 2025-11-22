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

png("results/RAW_Germination_vs_TOTAL_MASS.png")
plot(df$Germination, df$TOTAL_MASS)
dev.off()
# maybe what we can say is not that I counted germination right or wrong, but that the discrepency between expected and observed relationship btwn germination and total mass when germination equals 1...
# is *probably* that those seeds take longer to germinate, and were (for those w total mass > 100?) probably germinated more individuals after they were scored and thus actually belong in a different category

# for rows w low germination, remove the seed mass/number related traits but leave flowering time record
df[which(df$Germination == 1), 6] <- NA
# try removing just total mass... because the 100-seed weights is still consistent for the 6 lines w germination of 1...


#####################  calculate derived phenotypes  ########################
# fecundity
## fecundity = seed produced per plot
# fec is also our proxy for fitness then
## total seed weight / seed-weight-100/100
df2 <- df %>% 
    #mutate(MASS_PER_PLANT = TOTAL_MASS/Germination) %>% 
    mutate(SEED_COUNT = TOTAL_MASS / (SEED_WEIGHT_100/100)) %>%
    mutate(FECUNDITY = SEED_COUNT / Germination) 
write_delim(df2, "data/DERIVED_PHENOTYPES_FULL.tsv", "\t")


df2 %>% filter(Germination < 4)
df2 %>% filter(Germination < 4) %>% filter(FECUNDITY > 1000)


df3 <- df2 %>%
    select(-c(Germination, TOTAL_MASS, SEED_COUNT))



# calculate upper and lower bounds for each trait, in each year separately
up_bnds <- df3 %>% 
  group_by(Exp_year) %>% 
  summarise(across(-c(Genotype, Replicate), \(x) median(x, na.rm=T) + (3* IQR(x, na.rm=T))))

#df %>% 
#  group_by(Exp_year, Replicate) %>% 
#  summarise(across(-c(Genotype), \(x) median(x, na.rm=T) + (2* IQR(x, na.rm=T))))


lw_bnds <- df3 %>% 
  group_by(Exp_year) %>% 
  summarise(across(-c(Genotype, Replicate), \(x) median(x, na.rm=T) - (3* IQR(x, na.rm=T))))
  
#df %>% 
#  group_by(Exp_year, Replicate) %>% 
#  summarise(across(-c(Genotype), \(x) median(x, na.rm=T) - (2* IQR(x, na.rm=T))))

# filter for outlier values

# open png to print distribution plots
pdf("results/trait_outlier_distributions.pdf")

for(year in c(2022,2023)) {
  x <- df3 %>% filter(Exp_year == year)

  for(trait in 4:6){

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
    #x[which(x[,trait] > upp[[1]]), trait] <- NA
    #x[which(x[,trait] < dwnn[[1]]), trait] <- NA

    #g2 <- ggplot(x) + geom_histogram(aes(unlist(x[,trait]))) + 
    #  geom_vline(aes(xintercept=upp[[1]]), color = "#F31919") + 
    #  geom_vline(aes(xintercept=dwnn[[1]]), color = "#F31919") + 
    #  theme_bw() +
    #  labs(title=match, x=match, subtitle=year) 

    #print(g2)
  }

  # save the year's traits
  #assign(paste0("ayear_", year), x)

}

dev.off()

# re-join data from the two years of experiment
#df2 <- full_join(ayear_2022, ayear_2023)

df3 %>% filter(FECUNDITY > 1000)
df3[which(df3$FECUNDITY > 1000), 5:6] <- NA
write_delim(df3, "data/raw_phenotypes.tsv", "\t")

# any rows that are entirely NA?
df2[which(rowSums(is.na(df2))>4),]


# how many values were removed from each trait with the outlier filtering?
print("percent missing data for traits, after outlier filtering")
df22 <- df3 %>% filter(Exp_year==2022)
df23 <- df3 %>% filter(Exp_year==2023)

print('number of rows per year')
nrow(df22)
nrow(df23)

print("amount of missing data for the 3 traits")
df3 %>% group_by(Exp_year) %>% reframe(across(c(FT, SEED_WEIGHT_100, FECUNDITY), \(x) sum(is.na(x))))

print("percent of missing data per year")
df22 %>% summarise(across(c(FT, SEED_WEIGHT_100, FECUNDITY), \(x) (sum(is.na(x))/508)*100))

df23 %>% summarise(across(c(FT, SEED_WEIGHT_100, FECUNDITY), \(x) (sum(is.na(x))/510)*100))
