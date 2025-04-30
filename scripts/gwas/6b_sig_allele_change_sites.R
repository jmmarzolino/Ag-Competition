#!/usr/bin/env Rscript

#SBATCH --job-name=GWAS
#SBATCH -o /rhome/jmarz001/bigdata/Ag-Competition/scripts/gwas/6b_sig_allele_change_sites.stdout
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
colnames(df) <- c("CHR", "BP", "A1", "A2", #"F0A1", "F0A2",
 "F18A1", "F18A2", "F28A1", "F28A2", "F50A1", "F50A2", "F58A1", "F58A2")
# filter to only parents, F18, and F58 counts
df <- df %>% select(c(CHR, BP, A1, A2, #F0A1, F0A2, 
F18A1, F18A2, F58A1, F58A2))
# filter to relevant sites
#AF <- fread("A1_allele_freqs.tsv") %>% select(-c(F28_AF, F50_AF))
#df2 <- right_join(df, AF, by=c("CHR", "BP", "A1", "A2"))
df2 <- df

#(Xsq <- chisq.test(M))  # Prints test summary
#Xsq$observed   # observed counts (same as M)
#Xsq$expected   # expected counts under the null

df2$chi_p <- 1
for(i in 1:nrow(df2)){
  F18A1 <- df2[[i, 5]]
  F18A2 <- df2[[i, 6]]
  F58A1 <- df2[[i, 7]]
  F58A2 <- df2[[i, 8]]

  x <- matrix(c(F18A1, F18A2, F58A1, F58A2), nrow=2)
  dimnames(x) <- list(allele=c("A1", "A2"),
                      generation=c("F18", "F58"))

  y <- chisq.test(x)$p.value
  #y <- fisher.test(x)$p.value
  df2[[i, ncol(df2)]] <- y
}


print("number of sites w significant (p<0.05) allele change")
print(nrow(df2[which(df2$chi_p < 0.05),]))
print("fraction of sites w sig changed allele frequencies")
print((nrow(df2[which(df2$chi_p < 0.05),])) / nrow(df2))

print("number of sites w significant (p< bonferroni) allele change")
print(nrow(df2[which(df2$chi_p < (0.05/nrow(df2))),]))
print("fraction of sites w sig changed allele frequencies")
print((nrow(df2[which(df2$chi_p < (0.05/nrow(df2))),]))/nrow(df2))




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


#### allele freq change stats
changed <- minor %>% filter(chi_p < (0.005/nrow(minor)))

print("sig changed AF sites: P to F18 allele freq change distribution")
print(summary(changed$DELTA_P18))
print(summary(abs(changed$DELTA_P18)))

print("sig changed AF sites: F18 to F58 allele freq change distribution")
print(summary(changed$DELTA_1858))
print(summary(abs(changed$DELTA_1858)))




#### plotting
# plot distribution of allele frequency changes over generations
print('all sites: delta AF distribution')
print(quantile(minor$DELTA_1858))
#hist(df2$DELTA_1858)
g <- ggplot(minor, aes(DELTA_1858)) +
        geom_histogram(binwidth=0.05) +
        labs(x="", title="Change in Allele Frequency F18 to F58") +
        theme_bw(base_size=20)

g2 <- ggplot(minor, aes(abs(DELTA_1858))) +
        geom_histogram(binwidth=0.05) +
        labs(x="", title="Change in Allele Frequency F18 to F58") +
        theme_bw(base_size=20)


print('sig changed AF sites: delta AF distribution')
print(quantile(changed$DELTA_1858))
#hist(df3$DELTA_1858)
h <- ggplot(changed, aes(DELTA_1858)) +
        geom_histogram(binwidth=0.05) +
        labs(x="", title="Change in Allele Frequency F18 to F58", subtitle="Sites with Significantly Changed Frequency") +
        theme_bw(base_size=20)
h2 <- ggplot(changed, aes(abs(DELTA_1858))) +
        geom_histogram(binwidth=0.05) +
        labs(x="", title="Change in Allele Frequency F18 to F58", subtitle="Sites with Significantly Changed Frequency") +
        theme_bw(base_size=20)


gh <- ggarrange(g, h, nrow=2)
ggsave(gh, filename="sigsites_deltaAF_histograms.png", height=9.5, width=8, units="in", dpi="print")

gh2 <- ggarrange(g2, h2, nrow=2)
ggsave(gh2, filename="sigsites_deltaAF_histograms_absval.png", height=9.5, width=8, units="in", dpi="print")






#### beneficial alleles
### for sig changed sites, determine beneficial allele
changed$inc_or_dec <- sign(changed$DELTA_1858)
sum(changed$inc_or_dec)
table(changed$inc_or_dec, changed$minor_allele)

# frequencies are polarized to parent minor, so if freq inc, the minor is the beneficial allele
changed$beneficial_allele <- as.numeric()
# and if the freq of the minor dec, its NOT the beneficial allele
# if freq inc, copy number of minor allele
changed[which(changed$inc_or_dec==1 & changed$minor_allele==1), 19] <- 1
changed[which(changed$inc_or_dec==1 & changed$minor_allele==2), 19] <- 2

# if freq dec, NOT minor allele number
# if minor is A1, change beneficial num to 2
changed[which(changed$inc_or_dec==-1 & changed$minor_allele==1), 19] <- 2
# if minor is A2, change beneficial num to 1
changed[which(changed$inc_or_dec==-1 & changed$minor_allele==2), 19] <- 1
table(changed$beneficial_allele, changed$minor_allele)
table(changed$inc_or_dec, changed$beneficial_allele)


# re-join sig changed sites with neutral sites
changed <- changed %>% select(-c(inc_or_dec))
m2 <- full_join(minor, changed, by = join_by(CHR, BP, A1, A2, F0A1, F0A2, F18A1, F18A2, F58A1, F58A2, F0_AF, F18_AF, F58_AF, chi_p, minor_allele, DELTA_P18, DELTA_1858))



### join list of traits associated with sites (and thus, beneficial alleles)
tr <- fread("../raw_gwas_assoc_sig_sites.tsv")
tr$trait <- as.numeric(gsub("trait_(\\d+)", "\\1", tr$trait))

num <- fread("../num_sig_sites_per_trait_assoc_file.tsv")
# filter to relevant trait nums & names
num <- num %>% filter(no_sig_sites > 0) %>% select(c(trait_name, trait_num))

jn <- left_join(tr, num, by=c('trait'='trait_num')) %>% select(-'trait')
write_delim(jn, "gwas_sites_to_traits.tsv")


# combine list of traits that significantly changed w sig changed sites, change in AF over generations, and gwas site associated traits
## save minor allele freq w beneficial allele col + associated traits
m2$rs <- paste(m2$CHR, m2$BP, sep="_")
jn <- jn %>% select(-c(beta, p_lrt, allele0, allele1))
jnm <- right_join(jn, m2, by=c("chr"="CHR", "ps"="BP", "rs"))
write_delim(jnm, "MAF_BeneficialAllele_TraitAssoc.tsv")





# list of significantly changed traits
sig <- fread("../../002A_change_over_time_derived_traits_significant.tsv")
st <- unlist(sig$traits)


# How many sites w sig allele-freq change are associated with a) significantly changed traits, b) non-significant-changed traits, or c) no trait?
# look at significantly changed sites only
ham <- jnm %>% filter(chi_p < (0.005/nrow(minor)))
print("number of sites w significantly changed AF:")
print(length(unique(ham$rs)))


print("num of (unique) sig-changed sites associated w ANY traits")
ham %>% filter(!is.na(trait_name)) %>% select(rs) %>% unique %>% nrow %>% print


# 1 = sites associated w significantly changed trait
s1 <- ham %>% filter(trait_name %in% st) %>% select(c('chr', 'rs', 'ps', 'DELTA_1858', 'trait_name'))
print("num of sig-changed sites associated w sig-changed traits")
unique(s1$rs) %>% length %>% print


# 2 = sites associated w any other traits
s2 <- ham %>% filter(!is.na(trait_name)) %>% filter(!(trait_name %in% st)) %>% select(c('chr', 'rs', 'ps', 'DELTA_1858', 'trait_name'))
print("num of sig changed sites associated w NON-sig changed traits")
unique(s2$rs) %>% length %>% print


# 3 = all other sites
s3 <- ham %>% filter(is.na(trait_name)) %>% select(c('chr', 'rs', 'ps', 'DELTA_1858', 'trait_name'))
print("num of sig changed sites associated w NO trait")
unique(s3$rs) %>% length %>% print
