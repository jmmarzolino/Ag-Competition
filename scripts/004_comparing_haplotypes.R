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
# cut 3-digit genotype codes down to 2 digits
# to match family lines
df$Genotype <- gsub("(\\d+_\\d+)_\\d", "\\1", df$Genotype)
grep("\\d+_\\d+_\\d+", df$Genotype)


# how common was each haplotype in the 4 sampled generations?
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

# combine genotype, phenotype, and haplotype
hap_join <- inner_join(hap, df, by = "Genotype")

# check for missing genotype-haplotype match in table
which(is.na(hap_join$Haplotype))

# how many repeated haplotypes?
table(hap_join$Haplotype)
# not many repeat except 1 haplotype
sort(table(hap_join$Haplotype), decreasing=T)
sort(table(hap_join$Haplotype), decreasing=T) %>% table

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






## calculate generations' trait averages weighted by haplotype
# average phenotype for each haplotype
hap_trait_avg <- hap_join %>% group_by(Haplotype) %>% summarise(across(where(is.numeric), mean)) 

# join haplotype, generation frequency, and haplotype phenotypes
f18_hap_join <- right_join(hap_trait_avg, f18_hap_table, by="Haplotype")
f28_hap_join <- right_join(hap_trait_avg, f28_hap_table, by="Haplotype")
f50_hap_join <- right_join(hap_trait_avg, f50_hap_table, by="Haplotype")
f58_hap_join <- right_join(hap_trait_avg, f58_hap_table, by="Haplotype")


# copy rows to reflect their generation frequency number
for(j in c("f18_hap_join", "f28_hap_join", "f50_hap_join", "f58_hap_join")) {
  x <- get(j)
  out_df <- tibble(.rows=sum(x$Frequency))
  for(i in c(2:7)) {
    p <- rep((x[,i][[1]]), x$Frequency)
    out_df <- cbind(out_df, p)
  }
  colnames(out_df) <- colnames(x)[c(2:7)]
  assign(paste0(j, "_poprepd"), out_df)
}


# add col indicating phenotype's generation
f18_hap_join_poprepd$Generation <- 18
f28_hap_join_poprepd$Generation <- 28
f50_hap_join_poprepd$Generation <- 50
f58_hap_join_poprepd$Generation <- 58
## represent parental generation phenotypes at frequency
## with one entry for each haplotype
#table(hap_trait_avg$Haplotype)
f0_hap_join_poprepd <- hap_trait_avg %>% select(-Haplotype)
f0_hap_join_poprepd$Generation <- 0

## join generations together
x1 <- rbind(f18_hap_join_poprepd, f28_hap_join_poprepd)
x2 <- rbind(f50_hap_join_poprepd, f58_hap_join_poprepd)
joined_happops <- rbind(x1,x2)
joined_happops <- rbind(joined_happops,f0_hap_join_poprepd)

# write out haplotype-frequency-ajusted data
write_delim(joined_happops, "trait_BLUPs_HapRepPop.tsv")
