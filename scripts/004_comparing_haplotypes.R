#!/usr/bin/env Rscript
#SBATCH --mem=30G
#SBATCH --time=02:00:00
#SBATCH --output=/rhome/jmarz001/bigdata/Ag-Competition/scripts/004_comparing_haplotypes.stdout
#SBATCH -p koeniglab

library(tidyverse)
library(ggpubr)
library(data.table)

setwd("/bigdata/koeniglab/jmarz001/Ag-Competition/data")
source("../scripts/CUSTOM_FNS.R")
hap <- fread("hap_assign.txt")
df <- fread("trait_BLUPs.tsv")

hap$Haplotype <- as.factor(hap$Haplotype)



## calculate haplotype frequency per generation
f18_hap <- hap[grep("^1_\\d+", hap$Genotype),]
f18_hap_table <- data.frame(table(f18_hap$Haplotype))
colnames(f18_hap_table) <- c("Haplotype", "Frequency")
f18_hap_table$fraction <- f18_hap_table$Frequency / sum(f18_hap_table$Frequency)
f18_hap_table <- f18_hap_table %>% filter(Frequency > 0)


f28_hap <- hap[grep("^2_\\d+", hap$Genotype),]
f28_hap_table <- data.frame(table(f28_hap$Haplotype))
colnames(f28_hap_table) <- c("Haplotype", "Frequency")
f28_hap_table$fraction <- f28_hap_table$Frequency / sum(f28_hap_table$Frequency)
f28_hap_table <- f28_hap_table %>% filter(Frequency > 0)


f50_hap <- hap[grep("^3_\\d+", hap$Genotype),]
f50_hap_table <- data.frame(table(f50_hap$Haplotype))
colnames(f50_hap_table) <- c("Haplotype", "Frequency")
f50_hap_table$fraction <- f50_hap_table$Frequency / sum(f50_hap_table$Frequency)
f50_hap_table <- f50_hap_table %>% filter(Frequency > 0)


f58_hap <- hap[grep("^7_\\d+", hap$Genotype),]
f58_hap_table <- data.frame(table(f58_hap$Haplotype))
colnames(f58_hap_table) <- c("Haplotype", "Frequency")
f58_hap_table$fraction <- f58_hap_table$Frequency / sum(f58_hap_table$Frequency)
f58_hap_table <- f58_hap_table %>% filter(Frequency > 0)




## join phenotype data & haplotypes
# cut 3-digit genotype codes down to 2 digits
# to match family lines
df$Genotype <- gsub("(\\d+_\\d+)_\\d", "\\1", df$Genotype)
grep("\\d+_\\d+_\\d+", df$Genotype)

setdiff(df$Genotype, hap$Genotype)
# parents aren't in haplotype file
hap_join <- inner_join(hap, df, by = "Genotype")
# check for missing genotype-haplotype match in table
which(is.na(hap_join$Haplotype))

# how many repeated haplotypes?
table(hap_join$Haplotype)
# not many repeat except 1 haplotype
sort(table(hap_join$Haplotype), decreasing=T)
sort(table(hap_join$Haplotype), decreasing=T) %>% table



## calculate generations' trait averages weighted by haplotype

# subset data by generation
f18 <- hap_join[grep("^1_\\d+", hap_join$Genotype)]
f28 <- hap_join[grep("^2_\\d+", hap_join$Genotype)]
f50 <- hap_join[grep("^3_\\d+", hap_join$Genotype)]
f58 <- hap_join[grep("^7_\\d+", hap_join$Genotype)]


# average phenotype for each haplotype
f18_hap_trait_avg <- f18 %>% group_by(Haplotype) %>% summarise(across(where(is.numeric), mean))
f28_hap_trait_avg <- f28 %>% group_by(Haplotype) %>% summarise(across(where(is.numeric), mean))
f50_hap_trait_avg <- f50 %>% group_by(Haplotype) %>% summarise(across(where(is.numeric), mean))
f58_hap_trait_avg <- f58 %>% group_by(Haplotype) %>% summarise(across(where(is.numeric), mean))


# join haplotype, generation frequency, and haplotype phenotypes
f18_hap_join <- full_join(f18_hap_trait_avg, f18_hap_table)
f28_hap_join <- full_join(f28_hap_trait_avg, f28_hap_table)
f50_hap_join <- full_join(f50_hap_trait_avg, f50_hap_table)
f58_hap_join <- full_join(f58_hap_trait_avg, f58_hap_table)





# how many haplotypes that were present in the generations' sampling
# do not have phenotype measurements? (ie. weren't included in Ag experiment)
print("how many haplotypes don't have trait measurements?")
print("F18")
nrow(f18_hap_join[which(is.na(f18_hap_join$FT)), ])
print("F28")
nrow(f28_hap_join[which(is.na(f28_hap_join$FT)), ])
print("F50")
nrow(f50_hap_join[which(is.na(f50_hap_join$FT)), ])
print("F58")
nrow(f58_hap_join[which(is.na(f58_hap_join$FT)), ])


print("what fraction of the generation's population doesn't have measurements?")
print("F18")
sum(f18_hap_join[which(is.na(f18_hap_join$FT)), ]$fraction)
print("F28")
sum(f28_hap_join[which(is.na(f28_hap_join$FT)), ]$fraction)
print("F50")
sum(f50_hap_join[which(is.na(f50_hap_join$FT)), ]$fraction)
print("F58")
sum(f58_hap_join[which(is.na(f58_hap_join$FT)), ]$fraction)

# weighted trait averages table
gens <- c(18, 28, 50, 58)
trait <- colnames(f18_hap_join)[2:7]
wta <- tibble("generation" = sort(rep(gens, 6)), "trait"=rep(trait, 4), "value"=as.numeric(NA))


for(j in c('f18_hap_join', 'f28_hap_join', 'f50_hap_join', 'f58_hap_join')){
  print(j)
  x <- get(j)
  gn <- as.numeric(gsub("F(\\d+)_.*", "\\1", toupper(j)))

  for(i in c(2:7)){
    print(colnames(x)[i])
    trt <- colnames(x)[i]
    print(weighted.mean(unlist(x[,i]), x$fraction, na.rm=T))

    wta[which(wta$generation == gn & wta$trait==trt), 3] <- weighted.mean(unlist(x[,i]), x$fraction, na.rm=T)
  }
}


wta2 <- wta %>% pivot_wider(names_from="trait") 
write_delim(wta2, "weighted_generation_trait_avgs.tsv", "\t")
