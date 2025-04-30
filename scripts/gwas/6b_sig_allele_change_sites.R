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





# allele frequency and change between gens
df2$F18_AF <- df2$F18A1/(df2$F18A1+df2$F18A2)
df2$F58_AF <- df2$F58A1/(df2$F58A1+df2$F58A2)
df2$DELTA_AF <- df2$F58_AF - df2$F18_AF

#### allele freq change stats
changed <- df2 %>% filter(chi_p < (0.005/nrow(df2)))

## fixed sites
#df2[which(is.na(df2$chi_p)),]
# more ways... delta 18-58 is 0 & af18 & f58 are both 0 or both 1


#### plotting
# plot distribution of allele frequency changes over generations
print('all sites: delta AF distribution')
print(quantile(df2$DELTA_AF))
#hist(df2$DELTA_AF)
g <- ggplot(df2, aes(DELTA_AF)) +
        geom_histogram(binwidth=0.05) +
        labs(x="", title="Change in Allele Frequency F18 to F58") +
        theme_bw(base_size=20)


g2 <- ggplot(df2, aes(DELTA_AF)) +
        geom_histogram(binwidth=0.1) +
        labs(x="", title="Change in Allele Frequency F18 to F58") +
        theme_bw(base_size=20)

g3 <- ggplot(df2, aes(abs(DELTA_AF))) +
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
m2 <- full_join(df2, changed)
write_delim(m2, "progeny_AF_change.tsv")















##################################################################
############################################
############################################
######################
######################
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
s1 <- ham %>% filter(trait_name %in% st) %>% select(c('chr', 'rs', 'ps', 'DELTA_AF', 'trait_name'))
print("num of sig-changed sites associated w sig-changed traits")
unique(s1$rs) %>% length %>% print


# 2 = sites associated w any other traits
s2 <- ham %>% filter(!is.na(trait_name)) %>% filter(!(trait_name %in% st)) %>% select(c('chr', 'rs', 'ps', 'DELTA_AF', 'trait_name'))
print("num of sig changed sites associated w NON-sig changed traits")
unique(s2$rs) %>% length %>% print


# 3 = all other sites
s3 <- ham %>% filter(is.na(trait_name)) %>% select(c('chr', 'rs', 'ps', 'DELTA_AF', 'trait_name'))
print("num of sig changed sites associated w NO trait")
unique(s3$rs) %>% length %>% print
