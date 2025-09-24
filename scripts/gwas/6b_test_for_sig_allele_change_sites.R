#!/usr/bin/env Rscript

#SBATCH --job-name='test for allele change'
#SBATCH -o /rhome/jmarz001/bigdata/Ag-Competition/scripts/gwas/6b_test_for_sig_allele_change.stdout
#SBATCH --ntasks=1
#SBATCH --mem=80gb
#SBATCH -t 10:00:00
#SBATCH -p koeniglab

# which sites have significant change in allele freq (via allele counts)
library(pacman)
p_load(tidyverse, data.table, gridExtra, ggsci, ggpubr, Cairo)

setwd("/rhome/jmarz001/bigdata/Ag-Competition/results/gwas")

# start by using the allele counts to determine sites w sig AC change
df <- fread("COMBINED_AFS.txt")
colnames(df) <- c("CHR", "BP", "A1", "A2", "F0A1", "F0A2", "F18A1", "F18A2", "F28A1", "F28A2", "F50A1", "F50A2", "F58A1", "F58A2")
# filter to only parents, F18, and F58 counts
df <- df %>% select(c(CHR, BP, A1, A2, F0A1, F0A2, F18A1, F18A2, F58A1, F58A2))

# set default chi-square test as 1
# then run loop to format matrix of counts and test for sig diffs
df$chi_p_18_58 <- 1
for(i in 1:nrow(df)){
  F18A1 <- df[i, F18A1]
  F18A2 <- df[i, F18A2]
  F58A1 <- df[i, F58A1]
  F58A2 <- df[i, F58A2]

  x <- matrix(c(F18A1, F18A2, F58A1, F58A2), nrow=2)
  dimnames(x) <- list(allele=c("A1", "A2"),
                      generation=c("F18", "F58"))

  y <- chisq.test(x)$p.value
  #y <- fisher.test(x)$p.value
  df[i, 'chi_p_18_58'] <- y
}

# test for signifcant difference in allele counts between parents and F18 generation
# add another p-value column
df$chi_p_0_18 <- 1
for(i in 1:nrow(df)){
  F0A1 <- df[i, F0A1]
  F0A2 <- df[i, F0A2]
  F18A1 <- df[i, F18A1]
  F18A2 <- df[i, F18A2]

  x <- matrix(c(F0A1, F0A2, F18A1, F18A2), nrow=2)
  dimnames(x) <- list(allele=c("A1", "A2"),
                      generation=c("F0", "F18"))

  y <- chisq.test(x)$p.value
  #y <- fisher.test(x)$p.value
  df[i, 'chi_p_0_18'] <- y
}


print("fraction of sites w sig changed allele frequencies")
print("between F18 and F58")
print((nrow(df[which(df$chi_p_18_58 <= 0.05),])) / nrow(df))
print("between parents and F18")
print((nrow(df[which(df$chi_p_0_18 <= 0.05),])) / nrow(df))

print("fraction of sites w sig changed allele frequencies")
print("between F18 and F58")
print((nrow(df[which(df$chi_p_18_58 <= (0.05/nrow(df))),]))/nrow(df))
print("between parents and F18")
print((nrow(df[which(df$chi_p_0_18 <= (0.05/nrow(df))),]))/nrow(df))



# join chi-sq results w polarized allele frequencies
# and update allele frequency orientations...
allAF <- fread("all_sites_allele_freqs.tsv")

joind <- full_join(df, allAF, by=c("CHR", "BP", "A1"="A0", "A2"="A1"))
# write out results of chi-sq tests w allele frequencies
fwrite(df, "chisq_allele_freq_change_results.tsv")

# filter lines to sites w significant AF change 
thresh <- 0.05/(nrow(df))
# make one file for each significant generation period
df2 <- df %>% filter(chi_p_0_18 <= thresh | chi_p_18_58 <= thresh)
print('number of sites where either time period has significant changes:'); print(nrow(df2))

joindearly <- df %>% filter(chi_p_0_18 <= thresh)
joindlate <- df %>% filter(chi_p_18_58 <= thresh)

# and write out
fwrite(joindearly, "chisq_allele_freq_change_results_sig_earlygens.tsv")
fwrite(joindlate, "chisq_allele_freq_change_results_sig_lategens.tsv")


joinddearly <- df %>% filter(chi_p_0_18 <= thresh & chi_p_18_58 <= thresh)

fwrite(joinddearly, "chisq_allele_freq_change_results_sig_bothgens.tsv")



## read in list of suggestive sites and filter to see if any/all of them have chi-sq test significant allele freq change over gens
# list of top sites
top <- fread("gwas_top_sites.tsv")

tiptop <- left_join(top, df, by=c("chr"="CHR", "ps"="BP"))
tiptop$P_18_58_sig <- tiptop$chi_p_18_58 <= thresh 
tiptop$P_0_18_sig <- tiptop$chi_p_0_18 <= thresh 

tiptop %>% filter(sig_or_suggestive == "sig_site")