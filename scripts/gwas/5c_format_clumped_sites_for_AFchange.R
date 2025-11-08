#!/usr/bin/env Rscript

#SBATCH --job-name=manhattan
#SBATCH -o /rhome/jmarz001/bigdata/Ag-Competition/scripts/gwas/5c_format_clumped_sites_for_AFchange.stdout
#SBATCH --ntasks=1
#SBATCH --mem=20gb
#SBATCH -t 00:30:00
#SBATCH -p short

setwd("/bigdata/koeniglab/jmarz001/Ag-Competition/results/gwas")
library(pacman)
p_load(tidyverse, data.table, ggsci, Cairo, qqman)

#c("ASSOC_6_lmm.assoc.clumped", "ASSOC_7_lmm.assoc.clumped", "ASSOC_8_lmm.assoc.clumped")
#c("ASSOC_6_lmm.assoc.txt", "ASSOC_7_lmm.assoc.txt", "ASSOC_8_lmm.assoc.txt")

# read in file w lead snps list

for(i in 6:8){
    # read in full gwas file
    full <- fread(paste0("ASSOC_", i, "_lmm.assoc.txt"))

    # read in sig sites clumped list
    sig <- fread(paste0("ASSOC_", i, "_lmm.assoc.clumped"))[,1:5]

    # filter full gwas to sig sites
    filt <- full %>% filter(rs %in% sig$SNP) %>% select(c(chr, rs, ps, allele1, allele0, af, beta, p_lrt))

    # label associated trait
    filt$associated_trait <- i
    # save resulting list of sites
    assign(paste0("trait_", i), filt)
}

# join all your sites together into one list to pass to sig AF change analysis
joined <- full_join(trait_6, trait_7)
joined <- full_join(joined, trait_8)

# sub numbers for trait terms
joined$associated_trait <- gsub("6", "Flowering Time", joined$associated_trait)
joined$associated_trait <- gsub("7", "100-Seed Weight", joined$associated_trait)
joined$associated_trait <- gsub("8", "Fecundity", joined$associated_trait)


# mark which sites have a significant p-val from original gwas based on bonferroni threshold
bonn <- 0.05/nrow(full)
joined$sig_or_suggestive <- "suggestive_site"
joined[which(joined$p_lrt <= bonn), ncol(joined)] <- "sig_site"

# now save your resulting list of sites
write_delim(joined, "gwas_top_sites.tsv")

print(table(joined$associated_trait))
print(table(joined$associated_trait, joined$sig_or_suggestive))

repeated <- print(names(which(table(joined$rs)>1)))
joined[which(rs %in% repeated), ]