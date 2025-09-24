#!/usr/bin/env Rscript

#SBATCH --job-name='beta polarization'
#SBATCH -o /rhome/jmarz001/bigdata/Ag-Competition/scripts/gwas/6b_polarize_sites_to_gwas_beta.stdout
#SBATCH --ntasks=1
#SBATCH --mem=80gb
#SBATCH -t 10:00:00
#SBATCH -p koeniglab

library(pacman)
p_load(tidyverse, data.table, gridExtra, ggsci, ggpubr, Cairo)

setwd("/rhome/jmarz001/bigdata/Ag-Competition/results/gwas")


#### script purpose: 
#       calculate allele frequencies across all sites; 
#       for sites significantly associated w traits via gwas, polarize the allele frequencies to a positive trait change (beta value, calculated for allele change from a0 to a1);
#       plot allele frequencies over generations for gwas sites

#############
# read in data: 
# site allele counts for each generation 
afs <- fread("COMBINED_AFS.txt")
colnames(afs) <- c("CHR", "BP", "A0", "A1", "F0A0", "F0A1", "F18A0", "F18A1", "F28A0", "F28A1", "F50A0", "F50A1", "F58A0", "F58A1")

# calculate allele frequencies for A0/A1
# V5 = A0 count, V6 = A1 count...
afs$F0_AF <- afs$F0A0/(afs$F0A0 + afs$F0A1)
afs$F18_AF <- afs$F18A0/(afs$F18A0 + afs$F18A1)
afs$F28_AF <- afs$F28A0/(afs$F28A0 + afs$F28A1)
afs$F50_AF <- afs$F50A0/(afs$F50A0 + afs$F50A1)
afs$F58_AF <- afs$F58A0/(afs$F58A0 + afs$F58A1)

afs <- afs %>% select(-c( F0A0, F0A1, F18A0, F18A1, F28A0, F28A1, F50A0, F50A1, F58A0, F58A1))
fwrite(afs, "all_sites_allele_freqs.tsv")


### load list of gwas sites & combine w AFs
sig_sites <- fread("gwas_top_sites.tsv")
# filter to relevant cols
sig_sites <- sig_sites %>% select(c(chr, ps, allele1, allele0, beta, associated_trait))
# combine & filter list of all site AFs with list of gwas sites
joined_sites <- left_join(sig_sites, afs, by=c("chr"="CHR", "ps"="BP"))




##########        MAJOR ALLELE POLARIZATION
# check that A0/A1 letters match
# they don't! so we need to .... divide, re-order, rejoin...
x1 <- joined_sites[which(joined_sites$A1 != joined_sites$allele1), ]
x2 <- joined_sites[which(joined_sites$A1 == joined_sites$allele1), ] 
#A0 corresponds to allele0 and A1 corresponds to allele1

# sites w/out alleles matching
## polarize allele freqs to major allele
# switch A0/A1 & allele freq to match the allele0/allele1 orientation 
#(make ref/alt allele pulled from vcf match the allele0/allele1 notation from gemma)
# gemma uses major allele, ref/alt are pulled from vcf...
# aka just polarize A0/A1 and allele freq to major allele instead of ref/alt
x1[, grep("F\\d+_AF", colnames(x1))] <- x1 %>% reframe(across(starts_with("F"), \(x) 1-x))
# and switch allele states to match gwas allele0/allele1 format
x1$A0 <- x1$allele0 
x1$A1 <- x1$allele1

#rejoin data
maj_pol <- bind_rows(x1, x2)
# & check for matching alleles...
maj_pol[which(maj_pol$A1 != maj_pol$allele1), ]
# allele orientaiton for A0/A1 match 
# and all allele frequencies as of now are A0=major allele

# gemma defines beta as the phenotype change associated with allele change from allele0 (major allele) to allele1 (minor allele)




#############    BETA POLARIZATION
#######          polarize to gwas-beta orientation
### polarize change in allele frequency to the snp that causes positive trait change
# and record which allele is associated w positive beta

# divide data into sites w and w/out positive beta value
pos_beta <- maj_pol[which(sign(maj_pol$beta) == 1),]
neg_beta <- maj_pol[which(sign(maj_pol$beta) == -1),]


# polarize allele freq to the direction of gwas' beta

# for positive beta values, 
# record that allele1 increases beta
# keep beta value the same; flip allele frequency 
pos_beta$inc_beta_allele <- pos_beta$allele1
pos_beta[, grep("F\\d+_AF", colnames(pos_beta))] <- pos_beta %>% reframe(across(starts_with("F"), \(x) 1-x))


# for negative beta values, 
# record that allele0 increases beta
# flip beta value; keep allele frequency the same
neg_beta$inc_beta_allele <- neg_beta$allele0
neg_beta$beta <- neg_beta$beta * -1


# check that all allele1==A1 and allele0==A0
# now all allele frequencies are based on the count of the 'inc_beta_allele' out of the total allele counts & you can remove the A1/A0 columns
which(pos_beta$A0 != pos_beta$allele0)
which(neg_beta$A0 != neg_beta$allele0)
# rejoin data and drop A0/A1
rejoined <- bind_rows(neg_beta, pos_beta) %>% select(-c(A0, A1))
fwrite(rejoined, "beta_polarized_sig_sites_AFs.tsv")

# and record the site identified for multiple traits
rejoined$site <- paste0(rejoined$chr, ":", rejoined$ps)
print("any sites identified for multiple traits:")
dup_sites <- names(which(table(rejoined$site)>1))
rejoined %>% filter(site %in% dup_sites) %>% print






######################                  PLOTTING
#######        plot allele freq over gens for sites identified in gwas
rejoined <- rejoined %>% select(-c(chr, ps))
plotting <- rejoined %>% select(c(site, associated_trait, ends_with("AF"))) %>% pivot_longer(cols=c(F0_AF, F18_AF, F28_AF, F50_AF, F58_AF), values_to="allele_frequency", names_to="generation")
plotting$generation <- as.numeric(gsub("F(\\d+)_AF", "\\1", plotting$generation))

# there's one site in common between fecundity and flowering time, 
# so make a temp plotting column for 
plotting$tmp <- paste0(plotting$site, "_", plotting$associated_trait)

g <- ggplot(plotting, aes(x=generation, y=allele_frequency, group=tmp, color=associated_trait)) + 
  geom_point(alpha=0.4) + geom_line(alpha=0.4) +
  theme_bw() +
  labs(x="Generation", y="Allele Frequency", title="Allele Frequency over Generation", subtitle="site frequency polarized to trait increase", color="Associated Trait")
ggsave("trait_assoc_sites_allele_freq_over_gens.png", g, width=8, height=7)


# update the levels of associated trait as a factor to make plotting order of traits consistent
plotting$associated_trait <- factor(plotting$associated_trait,  # Change ordering manually
levels = c("Flowering Time", "100-Seed Weight", "Fecundity"))

g2 <- ggplot(plotting, aes(x=generation, y=allele_frequency, group=site)) + 
 geom_point(alpha=0.4) + geom_line(alpha=0.4) + 
 facet_wrap(~associated_trait) +
 theme_bw() +
 labs(x="Generation", y="Allele Frequency", title="Allele Frequency over Generation", subtitle="site frequency polarized to trait increase")
ggsave("trait_assoc_sites_allele_freq_over_gens_trait_faceted.png", g2, height=7, width=5.5*3)



##    now plot just generations F18 to F58...
#plotting2 <- plotting %>% filter(generation != 0)
#g <- ggplot(plotting2, aes(x=generation, y=allele_frequency, group=tmp, color=associated_trait)) + 
#  geom_point(alpha=0.4) + geom_line(alpha=0.4) +
#  theme_bw() +
#  labs(x="Generation", y="Allele Frequency", title="Allele Frequency over Generation", subtitle="site frequency polarized to trait increase", color="Associated Trait")
#ggsave("trait_assoc_sites_allele_freq_over_gens_18to58.png", g, width=8, height=7)

#g2 <- ggplot(plotting2, aes(x=generation, y=allele_frequency, group=site)) + 
# geom_point(alpha=0.4) + geom_line(alpha=0.4) + 
# facet_wrap(~associated_trait) +
# theme_bw() +
# labs(x="Generation", y="Allele Frequency", title="Allele Frequency over Generation", subtitle="site frequency polarized to trait increase")
#ggsave("trait_assoc_sites_allele_freq_over_gens_trait_faceted_18to58.png", g2, height=7, width=5.5*3)

# reload list of significant sites
sig_sites <- fread("gwas_top_sites.tsv")
# this time, keep only significant sites (remove suggestive sites)
sig_sites <- sig_sites %>% filter(sig_or_suggestive == "sig_site")
# format sites to match plotting...
sig_sites$site <- paste0(sig_sites$chr, ":", sig_sites$ps)

print("sites w p-val below 5e-07")
sig_sites %>% filter(p_lrt <= 5e-07)
print("sites w p-val below 5e-06")
sig_sites %>% filter(p_lrt <= 5e-06)

# select the list of sites from your significant site list
sigs <- sig_sites$site

# filter polarized data set to only significant sites
plot_sigs <- plotting %>% filter(site %in% sigs)
# plot allele freq over time with only significant/stronlgy suggestive? sites
g3 <- ggplot(plot_sigs, aes(x=generation, y=allele_frequency, group=site)) + 
 geom_point(alpha=0.4) + geom_line(alpha=0.4) + 
 facet_wrap(~associated_trait) +
 theme_bw() +
 ylim(c(0,1)) +
 labs(x="Generation", y="Allele Frequency", title="Allele Frequency over Generation", subtitle="site frequency polarized to trait increase")
ggsave("trait_assoc_sites_allele_freq_over_gens_trait_faceted_onlygwassignificantsites.png", g3, height=7, width=5.5*1)
