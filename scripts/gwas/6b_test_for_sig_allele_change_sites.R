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



df$chi_p_18_58 <- 1
for(i in 1:nrow(df)){
  F18A1 <- df[[i, 5]]
  F18A2 <- df[[i, 6]]
  F58A1 <- df[[i, 7]]
  F58A2 <- df[[i, 8]]

  x <- matrix(c(F18A1, F18A2, F58A1, F58A2), nrow=2)
  dimnames(x) <- list(allele=c("A1", "A2"),
                      generation=c("F18", "F58"))

  y <- chisq.test(x)$p.value
  #y <- fisher.test(x)$p.value
  df[[i, ncol(df)]] <- y
}

# test for signifcant difference in allele counts between parents and F18 generation
# add another p-value column
df$chi_p_0_18 <- 1
for(i in 1:nrow(df)){
  F0A1 <- df[[i, 5]]
  F0A2 <- df[[i, 6]]
  F18A1 <- df[[i, 7]]
  F18A2 <- df[[i, 8]]

  x <- matrix(c(F0A1, F0A2, F18A1, F18A2), nrow=2)
  dimnames(x) <- list(allele=c("A1", "A2"),
                      generation=c("F0", "F18"))

  y <- chisq.test(x)$p.value
  #y <- fisher.test(x)$p.value
  df[[i, ncol(df)]] <- y
}




print("number of sites w significant (p<0.05) allele change")
print(nrow(df[which(df$chi_p <= 0.05),]))
print("fraction of sites w sig changed allele frequencies")
print((nrow(df[which(df$chi_p <= 0.05),])) / nrow(df))

print("number of sites w significant (p< bonferroni) allele change")
print(nrow(df[which(df$chi_p <= (0.05/nrow(df))),]))
print("fraction of sites w sig changed allele frequencies")
print((nrow(df[which(df$chi_p <= (0.05/nrow(df))),]))/nrow(df))




### polarize progeny frequencies to parental minor allele frequency
# split df by parental minor allele 1 or 2
MAFA1 <- df2[which(df2$F0_AF < 0.5), ]
print("num of minor allele 1")
print(dim(MAFA1)[1])
MAFA2 <- df2[which(df2$F0_AF > 0.5), ]
print("num of minor allele 2")
print(dim(MAFA2)[1])

# invert AF to alt allele frequency if A1 is major
MAFA2[,11:13] <- 1-MAFA2[,11:13]
# label which allele is parental minor
MAFA1$minor_allele <- 1
MAFA2$minor_allele <- 2
# put df back together, in order
minor <- bind_rows(MAFA1, MAFA2)
minor <- minor[order(minor$CHR, minor$BP),]


# add change in AF to df
minor$DELTA_P18 <- minor$F18_AF - minor$F0_AF
minor$DELTA_1858 <- minor$F58_AF - minor$F18_AF


write_delim(minor, "MAF.tsv")



# calculate allele frequencies for each generation 
changed <- minor %>% filter(chi_p < (0.005/nrow(minor)))

print("sig changed AF sites: P to F18 allele freq change distribution")
print(summary(changed$DELTA_P18))
print(summary(abs(changed$DELTA_P18)))

print("sig changed AF sites: F18 to F58 allele freq change distribution")
print(summary(changed$DELTA_1858))
print(summary(abs(changed$DELTA_1858)))





df$F0_AF <- df$F0A1/(df$F0A1 + df$F0A2)
df$F18_AF <- df$F18A1/(df$F18A1+df$F18A2)
df$F58_AF <- df$F58A1/(df$F58A1+df$F58A2)
# calculate allele frequency change between gens
df$DELTA_AF_P <- df$F18_AF - df$F0_AF
df$DELTA_AF <- df$F58_AF - df$F18_AF

#### allele freq change stats
changed <- df %>% filter(chi_p <= (0.05/nrow(df)))


#### plotting
# plot distribution of allele frequency changes over generations
print('all sites: delta AF distribution')
print(quantile(df$DELTA_AF))
#hist(df$DELTA_AF)
g <- ggplot(df, aes(DELTA_AF)) +
        geom_histogram(binwidth=0.05) +
        labs(x="", title="Change in Allele Frequency F18 to F58") +
        theme_bw(base_size=20)


g2 <- ggplot(df, aes(DELTA_AF)) +
        geom_histogram(binwidth=0.1) +
        labs(x="", title="Change in Allele Frequency F18 to F58") +
        theme_bw(base_size=20)

g3 <- ggplot(df, aes(abs(DELTA_AF))) +
        geom_histogram(binwidth=0.05) +
        labs(x="", title="Change in Allele Frequency F18 to F58", subtitle="Absolute Values") +
        theme_bw(base_size=20)

print('sig changed AF sites: delta AF distribution')
print(quantile(changed$DELTA_AF))
#hist(df3$DELTA_AF)
h <- ggplot(changed, aes(DELTA_AF)) +
        geom_histogram(binwidth=0.05) +
        labs(x="", title="Change in Allele Frequency F18 to F58", subtitle="Sites with Significantly Changed Frequency") +
        theme_bw(base_size=20)
h2 <- ggplot(changed, aes(abs(DELTA_AF))) +
        geom_histogram(binwidth=0.05) +
        labs(x="", title="Change in Allele Frequency F18 to F58", subtitle="Sites with Significantly Changed Frequency") +
        theme_bw(base_size=20)


gh <- ggarrange(g, g3, h, h2)
ggsave(gh, filename="deltaAF_histograms.png", height=(7*2)+2, width=(7*2)+2, units="in", dpi="print")







################################# define which allele is beneficial / associated w increasing frequency over generations
#### beneficial alleles
### for sig changed sites, determine beneficial allele
changed$inc_or_dec <- sign(changed$DELTA_AF)
sum(changed$inc_or_dec)
table(changed$inc_or_dec)

# if freq inc, A1 is the beneficial allele
changed$beneficial_allele <- "A1"
# so if inc_or_dec col is = -1, change "A1" to "A2"
changed[which(changed$inc_or_dec==-1), ncol(changed)] <- "A2"

# join sig changed sites with neutral sites
m2 <- full_join(df, changed)
write_delim(m2, "sig_AF_change_test.tsv")


# filter to relevant sites
AF <- fread("A1_allele_freqs.tsv") %>% select(-c(F28_AF, F50_AF))
df <- right_join(df, AF, by=c("CHR", "BP", "A1", "A2"))
