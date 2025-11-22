#!/usr/bin/env Rscript
#SBATCH --job-name=allele_effect
#SBATCH -o /rhome/jmarz001/bigdata/Ag-Competition/scripts/gwas/6c_plot_allele_effect.stdout
#SBATCH --mem=20gb
#SBATCH -t 02:00:00
#SBATCH -p koeniglab

# plot phenotype associated with allele/genotype background
library(pacman)
p_load(tidyverse, data.table, gridExtra, ggsci, Cairo, EnvStats)

setwd("/rhome/jmarz001/bigdata/Ag-Competition/results/gwas")
source("../../scripts/CUSTOM_FNS.R")

## read in data sets
afs <- fread("progeny_AF_change.tsv")
# line IDs + genotype
gt <- fread("out.GT.FORMAT")
# line ID + phenotype
ph <- fread("all_traits.fam")

# filter genotype file by line IDs in phenotype file
lines <- ph$V1
gt <- gt %>% select(c("CHROM", "POS", any_of(lines)))


# for VRN candidate site, also plot allele effect on fecundity

## read in list of gwas sites
sites <- fread("ASSOC_6_lmm.assoc.clumped")[,1:11]
# filter to site of interest 
afs_i <- right_join(afs, sites, by=c("CHR"="CHR", "BP"="BP")) %>% filter(CHR=="chr4H" & BP==584673252)

# filter genotype data to common sites
gt_filt <- right_join(gt, afs_i, by=c("CHROM"="CHR", "POS"="BP")) 

# divide genotype subset file line IDs by genotype
tgs <- t(gt_filt)
a1 <- names(tgs[grep("0\\|0", tgs),])
a2 <- names(tgs[grep("1\\|1", tgs),])
a1a2 <- c(names(tgs[grep("1\\|0", tgs),]), names(tgs[grep("0\\|1", tgs),]))


a1t <- tibble("A1" = a1)
a1f <- right_join(ph, a1t, by=c("V1"="A1")) %>% select(all_of("V8"))
a1f$genotype <- "A1 homozygous"

a2t <- tibble("A2"=a2)
a2f <- right_join(ph, a2t, by=c("V1"="A2")) %>% select(all_of("V8"))
a2f$genotype <- "A2 homozygous"

jn <- bind_rows(a1f, a2f)

# there may be no heterozygous genotypes
# in which case, you can't join that data
if(length(a1a2)!=0) {
    a1a2t <- tibble("A1A2"=a1a2)
    a1a2f <- right_join(ph, a1a2t, by=c("V1"="A1A2")) %>% select(all_of("V8"))
    a1a2f$genotype <- "heterozygous"

    jn <- bind_rows(jn, a1a2f)

    # set the levels
    jn$genotype <- factor(jn$genotype, levels=c("A1 homozygous", "heterozygous", "A2 homozygous")) 
}

# count sample numbers
sample_counts <- jn %>% count(genotype) 
# difference in average A1 / A2 phenotype
grp_avg <- jn %>% group_by(genotype) %>% summarise(avg = mean(V8, na.rm=T))
a1a2_diff <- grp_avg[1, 2] - grp_avg[nrow(grp_avg), 2]

# plot
g <- ggplot(jn, aes(x=genotype, y= get("V8"))) + geom_boxplot() + stat_n_text() + labs(y="Fecundity", title=paste(afs_i$CHR, afs_i$BP), subtitle=paste0("difference in A1, A2 trait average: ", round(a1a2_diff, 3))) + theme_bw()
ggsave(paste0("allele_effect/sig_sites_FECUNDITY_", afs_i$CHR, "_", afs_i$BP, ".png"), g)
