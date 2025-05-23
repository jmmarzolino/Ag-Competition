#!/usr/bin/env Rscript
#SBATCH --job-name=GWAS
#SBATCH -o /rhome/jmarz001/bigdata/Ag-Competition/scripts/gwas/greenhouse-comparison/6c_plot_top_sites_allele_effect.stdout
#SBATCH --ntasks=1
#SBATCH --mem=80gb
#SBATCH -t 10:00:00
#SBATCH -p koeniglab

# for 'top' snps identified in analysis
# plot phenotype associated with each allele/genotype
library(pacman)
p_load(tidyverse, data.table, gridExtra, ggsci, Cairo)

setwd("/rhome/jmarz001/bigdata/Ag-Competition/results/gwas")



## read in data sets
# line IDs + genotype
gt <- fread("out.GT.FORMAT")
# format that data...



# line ID + phenotype
ph <- fread("all_traits.fam")
# separate phenotype cols?

# combine data sets
join <- full_join(gt, ph)



for (i in 6:11) {

    # filtering file: 'top' snps
    tops <- fread(paste0("ASSOC_", i, "_top_sites.txt"))
    # separate list of top sites for each gwas
    tops_i <- tdf filter by tops...

    inner_join(join, tops)

}







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
AF <- df %>% select(c(CHR, BP, A1, A2, #F0_AF, F18_AF, 
F28_AF, F50_AF, F58_AF))

# filter sites that don't segregate in progeny
print("summary of un-filtered progeny allele counts")
print(summary(AF[, 6:9]))
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
#AF <- AF[-which(AF$avg < 0.0015),]
AF[order(AF$summ, decreasing=F),]

#
print(paste0("number of sites fixed at nearly 1 (x > 0.995) in progeny: ", length(which(AF$avg > 0.995))))
AF[order(AF$summ, decreasing=T),]
#AF <- AF[-which(AF$avg > 0.995),]
# I'm leaving sites that maybe don't segregate because the values will be different in the larger population

AF <- AF %>% select(-c(summ, avg))
print(paste0("number of sites remaining: ", nrow(AF)))
write_delim(AF, "A1_allele_freqs.tsv")
