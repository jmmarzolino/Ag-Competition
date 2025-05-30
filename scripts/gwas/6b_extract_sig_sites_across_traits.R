#!/usr/bin/env Rscript
#SBATCH -o /rhome/jmarz001/bigdata/Ag-Competition/scripts/gwas/6b_extract_sig_sites_across_traits.stdout
#SBATCH --ntasks=1
#SBATCH --mem=64gb
#SBATCH -t 01:00:00
#SBATCH -p koeniglab

library(tidyverse)
library(data.table)

setwd("/rhome/jmarz001/bigdata/Ag-Competition/results/gwas")
source("../../scripts/CUSTOM_FNS.R")


file_list <- fread("association_files_traits.txt")

for(i in 1:nrow(file_list)){

    # read in file
    df <- fread(file_list[[i, 3]])

    # extract the top x% of sites from gwas
    # define bonnferroni threshhold for significance with multiple testing
    bonn <- 0.05 / nrow(df)

    # keep only sites at or below our threshhold
    df <- df %>% filter(p_lrt <= bonn)

    # label associated trait
    df$associated_trait <- file_list[[i, 1]]

    # now save your resulting list of sites
    assign(gsub("(ASSOC_\\d+).assoc.txt", "\\1", file_list[[i, 3]]), df)
}



lm_files <- full_join(ASSOC_6, ASSOC_7)
lm_files <- full_join(lm_files, ASSOC_8)
lm_files <- full_join(lm_files, ASSOC_9)
lm_files <- full_join(lm_files, ASSOC_10)
lm_files <- full_join(lm_files, ASSOC_11)

lm_files$gwas_method <- "lm"


# repeat for lmm model files
# with file name ref col changed & assigned name change
for(i in 1:nrow(file_list)){

    # read in file
    df <- fread(file_list[[i, 4]])

    # extract the top x% of sites from gwas
    # define bonnferroni threshhold for significance with multiple testing
    bonn <- 0.05 / nrow(df)

    # keep only sites at or below our threshhold
    df <- df %>% filter(p_lrt <= bonn)

    # label associated trait
    df$associated_trait <- file_list[[i, 1]]

    # now save your resulting list of sites
    assign(gsub("(ASSOC_\\d+_lmm).assoc.txt", "\\1", file_list[[i, 4]]), df)
}

lmm_files <- full_join(ASSOC_6_lmm, ASSOC_7_lmm)
lmm_files <- full_join(lmm_files, ASSOC_8_lmm)
lmm_files <- full_join(lmm_files, ASSOC_9_lmm)
lmm_files <- full_join(lmm_files, ASSOC_10_lmm)
lmm_files <- full_join(lmm_files, ASSOC_11_lmm)

lmm_files$gwas_method <- "lmm"


df2 <- full_join(lm_files, lmm_files)
df2$associated_trait <- tidy_text_substitution(df2$associated_trait)

write_delim(df2, "all_gwas_sig_sites.tsv")
print(paste0("number of significant sites across traits and gwas models: ", length(df2$rs)))
print(paste0("number of significant sites across traits: ", length(unique(df2$rs))))

print("number of sig sites per trait")
print(table(df2$associated_trait))

print("number of significant sites for more than one trait (if any. traits w the same sig site twice indicate significance with 2 gwas models): ",)
table(df2$rs, df2$associated_trait)[which(rowSums(table(df2$rs, df2$associated_trait)) > 1), ]


print("are there any sites overlapping trait/method?")
print("differing numbers indicate site overlap")
tmp <- df2 %>% select(c(chr, rs, ps, associated_trait, gwas_method))
length(tmp$chr)
length(unique(tmp$rs))


print("are there any sites that are sig for multiple traits?")
print("if numbers match, there are no overlaping sites")
tmp2 <- df2 %>% select(c(chr, rs, ps, associated_trait)) %>% unique
length(tmp2$rs)
length(unique(tmp2$rs))




####   CHANGE IN ALLELE FREQUENCY OVER GENERATIONS
# make list of significant sites
sig_sites <- df2 %>% select(c(chr, rs, ps, allele1, allele0, beta, associated_trait))

# load site allele freuqencies for each generation 
afs <- fread("COMBINED_AFS.txt")
# calculate allele frequencies for A1/A2
# V5 = A1 count, V6 = A2 count...
afs$F0_AF <- afs$V5/(afs$V5 + afs$V6)
afs$F18_AF <- afs$V7/(afs$V7 + afs$V8)
afs$F28_AF <- afs$V9/(afs$V9 + afs$V10)
afs$F50_AF <- afs$V11/(afs$V11 + afs$V12)
afs$F58_AF <- afs$V13/(afs$V13 + afs$V14)
afs <- afs %>% select(-c(V5, V6, V7, V8, V9, V10, V11, V12, V13, V14))

# combine list of sig sites w all sites alllele freqs
df3 <- left_join(sig_sites, afs, by=c("chr"="V1", "ps"="V2")) %>% select(-c("rs"))

# check that A1/A2 letters match
# they don't! so we need to ....
x1 <- df3[which(df3$V4 != df3$allele1), ]
x2 <- df3[which(df3$V4 == df3$allele1), ]


# sites w/out alleles matching
# switch A1/A2 & allele freq to match the a1/a2 orientation of other sites...
#(make ref/alt allele consistent to A1/A2 pulled from plink (which uses major allele) instead of vcfs a1/a2 allele orientation) 
x1$F18_AF <- 1 - x1$F18_AF
x1$F28_AF <- 1 - x1$F28_AF
x1$F50_AF <- 1 - x1$F50_AF
x1$F58_AF <- 1 - x1$F58_AF
a1 <- x1$V4
a2 <- x1$V3

x1$V3 <- a1
x1$V4 <- a2

#rejoin data
df4 <- bind_rows(x1, x2)
# & check for matching alleles...
df4[which(df4$V4 != df4$allele1), ]
# allele orientaiton for v3/v4 match that of other sites now...
# which means v3 = ref allele (allele0) & v4 = alt allele (allele1)
# and beta reflects the phenotype change from ref to alt


# record which allele is associated w positive beta
pos_beta <- df4[which(sign(df4$beta) == 1),]
neg_beta <- df4[which(sign(df4$beta) == -1),]


### polarize change in allele frequency to the snp that causes positive trait change
# polarize allele freq to the direction of gwas' beta
## then polarize allele frequencies to allele which increases phenotype (positive beta, flip sign for negative beta)
# for positive beta, record that alt allele is beneficial
pos_beta$beneficial_allele <- pos_beta$allele1
neg_beta$beneficial_allele <- neg_beta$allele0

# now, for negative beta values, switch them to positive values
neg_beta$beta <- abs(neg_beta$beta)
# and switch their allele freq orientation
neg_beta$F18_AF <- 1 - neg_beta$F18_AF
neg_beta$F28_AF <- 1 - neg_beta$F28_AF
neg_beta$F50_AF <- 1 - neg_beta$F50_AF
neg_beta$F58_AF <- 1 - neg_beta$F58_AF

# rejoin data and you can now drop v3/v4
df5 <- bind_rows(pos_beta, neg_beta) %>% select(-c(V3, V4))



#####     plot allele freq over gens for sites identified in gwas
df5$site <- paste0(df5$chr, ":", df5$ps)
df5 <- df5 %>% select(-c(chr, ps))
plotting <- df5 %>% select(c(site, associated_trait, ends_with("AF"))) %>% pivot_longer(cols=c(F18_AF, F28_AF, F50_AF, F58_AF), values_to="allele_frequency", names_to="generation")
plotting$generation <- as.numeric(gsub("F(\\d+)_AF", "\\1", plotting$generation))

g <- ggplot(plotting, aes(x=generation, y=allele_frequency, group=site, color=associated_trait)) + 
  geom_point(alpha=0.6) + geom_line(alpha=0.6) +
  theme_bw() +
  labs(x="Generation", y="Allele Frequency", title="Allele Frequency over Generation", subtitle="site frequency polarized to trait increase", color="Associated Trait")
ggsave("trait_assoc_sites_allele_freq_over_gens.png", g, width=8, height=7)

g2 <- ggplot(plotting, aes(x=generation, y=allele_frequency, group=site)) + 
 geom_point(alpha=0.6) + geom_line(alpha=0.6) + 
 facet_wrap(~associated_trait) +
 theme_bw() +
 labs(x="Generation", y="Allele Frequency", title="Allele Frequency over Generation", subtitle="site frequency polarized to trait increase")
ggsave("trait_assoc_sites_allele_freq_over_gens_trait_faceted.png", g2, height=7, width=5.5*3)
