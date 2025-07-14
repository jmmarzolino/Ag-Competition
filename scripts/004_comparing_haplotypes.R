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

# set haplotypes to factors instead of numbers
hap$Haplotype <- as.factor(hap$Haplotype)
# cut 3-digit genotype codes down to 2 digits
# to match family lines
df$Genotype <- gsub("(\\d+_\\d+)_\\d", "\\1", df$Genotype)
grep("\\d+_\\d+_\\d+", df$Genotype)


# how common was each haplotype in the 4 sampled generations?
## calculate haplotype frequency per generation
# select line IDs from each generation
f18_hap <- hap[grep("^1_\\d+", hap$Genotype),]
# quantify how common each haplotype was in that generation
f18_hap_table <- data.frame(table(f18_hap$Haplotype))
colnames(f18_hap_table) <- c("Haplotype", "Frequency")
# divide haplotype occurance by the total number of haplotypes in the generation
f18_hap_table$fraction <- f18_hap_table$Frequency / sum(f18_hap_table$Frequency)
# filter out hap numbers not in the generation (were included in table b/c of hap levels)
f18_hap_table <- f18_hap_table %>% filter(Frequency > 0)

# repeat for other generations
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

# combine genotype, phenotype, and haplotype
hap_join <- inner_join(hap, df, by = "Genotype")

# how many repeated haplotypes?
print("table of repeated haplotypes")
sort(table(hap_join$Haplotype), decreasing=T)
# not many repeat except 1 haplotype
print("table of haplotype re-occurance")
sort(table(hap_join$Haplotype), decreasing=T) %>% table


# average phenotype for each haplotype
hap_trait_avg <- hap_join %>% group_by(Haplotype) %>% summarise(across(where(is.numeric), mean)) 

# join haplotype, generation frequency, and haplotype phenotypes
f18_hap_join <- right_join(hap_trait_avg, f18_hap_table, by="Haplotype")
f28_hap_join <- right_join(hap_trait_avg, f28_hap_table, by="Haplotype")
f50_hap_join <- right_join(hap_trait_avg, f50_hap_table, by="Haplotype")
f58_hap_join <- right_join(hap_trait_avg, f58_hap_table, by="Haplotype")


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




# copy rows to reflect their generation frequency number
for(j in c("f18_hap_join", "f28_hap_join", "f50_hap_join", "f58_hap_join")) {
  x <- get(j)
  out_df <- tibble(.rows=sum(x$Frequency))
  for(i in c(2:4)) {
    p <- rep((x[,i][[1]]), x$Frequency)
    out_df <- cbind(out_df, p)
  }
  colnames(out_df) <- colnames(x)[c(2:4)]
  assign(paste0(j, "_poprepd"), out_df)
}


# add col indicating phenotype's generation
f18_hap_join_poprepd$Generation <- 18
f28_hap_join_poprepd$Generation <- 28
f50_hap_join_poprepd$Generation <- 50
f58_hap_join_poprepd$Generation <- 58
## represent parental/F1 generation phenotypes at frequency
## with one entry for each haplotype you have
#table(hap_trait_avg$Haplotype)
f0_hap_join_poprepd <- hap_trait_avg %>% select(-Haplotype)
f0_hap_join_poprepd$Generation <- 1

## join generations together
x1 <- rbind(f18_hap_join_poprepd, f28_hap_join_poprepd)
x2 <- rbind(f50_hap_join_poprepd, f58_hap_join_poprepd)
joined_happops <- rbind(x1,x2)
joined_happops <- rbind(joined_happops,f0_hap_join_poprepd)

# remove NA columns?
joined_happops <- joined_happops[-which(is.na(joined_happops$FT)),]
# write out haplotype-frequency-ajusted data
write_delim(joined_happops, "trait_BLUPs_HapRepPop.tsv")
