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


# read in data: 
# site allele counts for each generation 
afs <- fread("COMBINED_AFS.txt")
colnames(afs) <- c("CHR", "BP", "A1", "A2", "F0A1", "F0A2", "F18A1", "F18A2", "F28A1", "F28A2", "F50A1", "F50A2", "F58A1", "F58A2")

# calculate allele frequencies for A1/A2
# V5 = A1 count, V6 = A2 count...
afs$F0_AF <- afs$F0A1/(afs$F0A1 + afs$F0A2)
afs$F18_AF <- afs$F18A1/(afs$F18A1 + afs$F18A2)
afs$F28_AF <- afs$F28A1/(afs$F28A1 + afs$F28A2)
afs$F50_AF <- afs$F50A1/(afs$F50A1 + afs$F50A2)
afs$F58_AF <- afs$F58A1/(afs$F58A1 + afs$F58A2)

afs <- afs %>% select(-c( F0A1, F0A2, F18A1, F18A2, F28A1, F28A2, F50A1, F50A2, F58A1, F58A2))
fwrite(afs, "all_sites_allele_freqs.tsv")


### load list of gwas sites & combine w AFs
sig_sites <- fread("gwas_top_sites.tsv")
# filter to relevant cols
sig_sites <- sig_sites %>% select(c(chr, ps, allele1, allele0, beta, associated_trait))
# combine & filter list of all site AFs with list of gwas sites
joined_sites <- left_join(sig_sites, afs, by=c("chr"="CHR", "ps"="BP"))



#############
# check that A1/A2 letters match
# they don't! so we need to .... divide, re-order, rejoin...
x1 <- joined_sites[which(joined_sites$A2 != joined_sites$allele1), ]
x2 <- joined_sites[which(joined_sites$A2 == joined_sites$allele1), ] #A1 corresponds to allele0 and A2 corresponds to allele1

# sites w/out alleles matching
# switch A1/A2 & allele freq to match the a1/a2 orientation of other sites...
#(make ref/alt allele consistent to A1/A2 pulled from plink (which uses major allele) instead of vcfs a1/a2 allele orientation) 
x1$F0_AF <- 1 - x1$F0_AF
x1$F18_AF <- 1 - x1$F18_AF
x1$F28_AF <- 1 - x1$F28_AF
x1$F50_AF <- 1 - x1$F50_AF
x1$F58_AF <- 1 - x1$F58_AF
# switch allele states to match other set of sites...
a2 <- x1$A2
a1 <- x1$A1

x1$allele0 <- a1
x1$allele1 <- a2

#rejoin data
df4 <- bind_rows(x1, x2)
# & check for matching alleles...
df4[which(df4$A2 != df4$allele1), ]
# allele orientaiton for a1/A2 match that of other sites now...
# which means v3 = ref allele (allele0) & A2 = alt allele (allele1)
# and beta reflects the phenotype change from ref to alt



#######          polarize to gwas-beta orientation
### polarize change in allele frequency to the snp that causes positive trait change
# record which allele is associated w positive beta
pos_beta <- df4[which(sign(df4$beta) == 1),]
neg_beta <- df4[which(sign(df4$beta) == -1),]

# polarize allele freq to the direction of gwas' beta
## then polarize allele frequencies to allele which increases phenotype (positive beta, flip sign for negative beta)
# for positive beta, record that alt allele is beneficial
pos_beta$beneficial_allele <- pos_beta$allele1
neg_beta$beneficial_allele <- neg_beta$allele0

# now, for negative beta values, switch them to positive values
neg_beta$beta <- abs(neg_beta$beta)
# and switch their allele freq orientation
neg_beta$F0_AF <- 1 - neg_beta$F0_AF
neg_beta$F18_AF <- 1 - neg_beta$F18_AF
neg_beta$F28_AF <- 1 - neg_beta$F28_AF
neg_beta$F50_AF <- 1 - neg_beta$F50_AF
neg_beta$F58_AF <- 1 - neg_beta$F58_AF

# rejoin data and you can now drop A1/A2
df5 <- bind_rows(neg_beta, pos_beta) %>% select(-c(A1, A2))
fwrite(df5, "beta_polarized_sig_sites_AFs.tsv")



######################                  PLOTTING
#######        plot allele freq over gens for sites identified in gwas
df5$site <- paste0(df5$chr, ":", df5$ps)
df5 <- df5 %>% select(-c(chr, ps))
plotting <- df5 %>% select(c(site, associated_trait, ends_with("AF"))) %>% pivot_longer(cols=c(F0_AF, F18_AF, F28_AF, F50_AF, F58_AF), values_to="allele_frequency", names_to="generation")
plotting$generation <- as.numeric(gsub("F(\\d+)_AF", "\\1", plotting$generation))

g <- ggplot(plotting, aes(x=generation, y=allele_frequency, group=site, color=associated_trait)) + 
  geom_point(alpha=0.4) + geom_line(alpha=0.4) +
  theme_bw() +
  labs(x="Generation", y="Allele Frequency", title="Allele Frequency over Generation", subtitle="site frequency polarized to trait increase", color="Associated Trait")
ggsave("trait_assoc_sites_allele_freq_over_gens.png", g, width=8, height=7)

g2 <- ggplot(plotting, aes(x=generation, y=allele_frequency, group=site)) + 
 geom_point(alpha=0.4) + geom_line(alpha=0.4) + 
 facet_wrap(~associated_trait) +
 theme_bw() +
 labs(x="Generation", y="Allele Frequency", title="Allele Frequency over Generation", subtitle="site frequency polarized to trait increase")
ggsave("trait_assoc_sites_allele_freq_over_gens_trait_faceted.png", g2, height=7, width=5.5*3)



##    now plot just generations F18 to F58...
plotting2 <- plotting %>% filter(generation != 0)
g <- ggplot(plotting2, aes(x=generation, y=allele_frequency, group=site, color=associated_trait)) + 
  geom_point(alpha=0.4) + geom_line(alpha=0.4) +
  theme_bw() +
  labs(x="Generation", y="Allele Frequency", title="Allele Frequency over Generation", subtitle="site frequency polarized to trait increase", color="Associated Trait")
ggsave("trait_assoc_sites_allele_freq_over_gens_18to58.png", g, width=8, height=7)

g2 <- ggplot(plotting2, aes(x=generation, y=allele_frequency, group=site)) + 
 geom_point(alpha=0.4) + geom_line(alpha=0.4) + 
 facet_wrap(~associated_trait) +
 theme_bw() +
 labs(x="Generation", y="Allele Frequency", title="Allele Frequency over Generation", subtitle="site frequency polarized to trait increase")
ggsave("trait_assoc_sites_allele_freq_over_gens_trait_faceted_18to58.png", g2, height=7, width=5.5*3)
