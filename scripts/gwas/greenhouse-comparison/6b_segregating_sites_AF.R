#!/usr/bin/env Rscript
#SBATCH --job-name=GWAS
#SBATCH -o /rhome/jmarz001/bigdata/Ag-Competition/scripts/gwas/greenhouse-comparison/6b_segregating_sites_AF.stdout
#SBATCH --ntasks=1
#SBATCH --mem=80gb
#SBATCH -t 10:00:00
#SBATCH -p koeniglab

# produce list of sites segregating in IPK sample
library(pacman)
p_load(tidyverse, data.table, gridExtra, ggsci, Cairo)

setwd("/rhome/jmarz001/bigdata/Ag-Competition/results/gwas")

# genotypes based on
##FORMAT=<ID=GT,Number=1,Type=String,Description="Best Guessed Genotype with posterior probability threshold of 0.9">
df <- fread("COMBINED_AFS.txt")
colnames(df) <- c("CHR", "BP", "A1", "A2", #"F0A1", "F0A2", 
"F18A1", "F18A2", "F28A1", "F28A2", "F50A1", "F50A2", "F58A1", "F58A2")

## filter to population subset's sites
#gwas_sites <- fread("prog_sites.txt")
#colnames(gwas_sites) <- c("CHR", "BP")
#print(paste0("number of sites used in gwas: ", dim(gwas_sites)[1]))
#print(paste0("number of sites in vcf: ", dim(df)[1]))


############ UGH WHAT THE ACTUAL FUCK
## damn vcfs
#tmp <- full_join(gwas_sites, df, by=c('CHR', 'BP'))
#tmp[which(rowSums(is.na(tmp[,6:14])) == 9),]

# filter full df to gwas sites
#df <- inner_join(gwas_sites, df, by=c('CHR', 'BP'))
#print(paste0("number of sites retained after join: ", dim(df)[1]))

# calculate allele frequencies
#df$F0_AF <- df$F0A1 / (df$F0A1 + df$F0A2)
df$F18_AF <- df$F18A1 / (df$F18A1 + df$F18A2)
df$F28_AF <- df$F28A1 / (df$F28A1 + df$F28A2)
df$F50_AF <- df$F50A1 / (df$F50A1 + df$F50A2)
df$F58_AF <- df$F58A1 / (df$F58A1 + df$F58A2)
# not making it polarized to parent minor allele since this is only for filtering segregating sites. so I just need to know if any frequencies are 0/100
AF <- df %>% select(c(CHR, BP, A1, A2, #F0_AF, 
F18_AF, F28_AF, F50_AF, F58_AF))

# filter sites that don't segregate in progeny
print("summary of un-filtered progeny allele counts")
print(summary(AF[, 6:ncol(AF)]))
# there's definietly sites at 0 and 1...

# how many sites are fixed in F18?
AF[which(AF$F18_AF == 0),]
AF[which(AF$F18_AF == 0),] %>% nrow
# but not all sites at 0 in 18 are at 0 in other generations... site not fully fixed?

AF$summ <- AF$F18_AF + AF$F28_AF + AF$F50_AF + AF$F58_AF
AF$avg <- AF$summ/4
# sites with progeny AF sums or 0 or 4 are fixed (0*4 or 1*4)
summary(AF$summ)
summary(AF$avg)
# values range from 0 to *nearly* but not quite 4

# definitely filter all 0s
print(paste0("number of sites fixed at 0 in progeny: ", length(which(AF$summ==0))))
AF <- AF[-which(AF$summ==0),]

# set a MAF of 0.05? or lower? (0.01*4 = 0.04), (0.99*4=3.96)
# worth investigating further
# sites w average progeny AF < rounded 1% of samples
#AF[which((AF$summ)/4 < 0.004),]
# print(paste0("number of sites nearly 0 in progeny... ", length(which((AF$summ)/4 < 0.004))))
print(paste0("number of sites nearly fixed at 0 in progeny (x < 0.0015): ", length(which(AF$avg < 0.0015))))
AF <- AF[-which(AF$avg < 0.0015),]
#AF[order(AF$summ, decreasing=F),]

#
print(paste0("number of sites fixed at nearly 1 (x > 0.995) in progeny: ", length(which(AF$avg > 0.995))))
#AF[order(AF$summ, decreasing=T),]
#AF <- AF[-which(AF$avg > 0.995),]
# I'm leaving sites that maybe don't segregate because the values will be different in the larger population

AF <- AF %>% select(-c(summ, avg))
print(paste0("number of sites remaining: ", nrow(AF)))
write_delim(AF, "A1_allele_freqs.tsv")
