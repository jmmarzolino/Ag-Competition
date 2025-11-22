#!/usr/bin/env Rscript
#SBATCH --job-name=GWAS
#SBATCH -o /rhome/jmarz001/bigdata/Ag-Competition/scripts/gwas/6c_plot_top_sites_allele_effect.stdout
#SBATCH --ntasks=1
#SBATCH --mem=80gb
#SBATCH -t 10:00:00
#SBATCH -p koeniglab

# for 'top' snps identified in analysis
# plot phenotype associated with each allele/genotype
library(pacman)
p_load(tidyverse, data.table, gridExtra, ggsci, Cairo, EnvStats)

setwd("/rhome/jmarz001/bigdata/Ag-Competition/results/gwas")
source("../../scripts/CUSTOM_FNS.R")

##################################################
### read in data sets
# allele frequencies
afs <- fread("progeny_AF_change.tsv")
# list of top sites
top <- fread("gwas_top_sites.tsv")
# filter afs to just top sites
afs <- left_join(top, afs, by=c("chr"="CHR", "ps"="BP"))


# line IDs + genotype
gt <- fread("top_sites.gt")
gtnm <- fread("combined_filt.gt.names", header=F)
# replace first 28 parent-line names w their KL numbers (45 to 72)
gtnm[1:28, 1] <- 45:72
# pair up genotype names with columns of genotype data
colnames(gt) <- c("CHROM", "POS", gtnm$V1)


# line ID + phenotype
# read in full phenotype file
pheno <- fread("../../data/BLUPs.tsv")

# shorten genotype codes to 2 generations (\\d_\\d), for consistency, IDing genotype families instead of lines
pheno$Genotype <- gsub("(\\d+_\\d+)_\\d+", "\\1", pheno$Genotype)
# add generation column...
pheno <- add_generation(pheno)
# if generation is 0, remove one more _ of genotype ID
pheno[which(pheno$Generation == 0), ]$Genotype <- gsub("(\\d+)_\\d", "\\1", pheno[which(pheno$Generation == 0), ]$Genotype)

# make sure trait names match later
colnames(pheno) <- tidy_text_substitution(colnames(pheno))


# filter genotype file by line IDs in phenotype file
lines <- pheno$Genotype
gt <- gt %>% select(c("CHROM", "POS", any_of(lines)))
#colnames(gt)[3:ncol(gt)]
##################################################


for(i in c(unique(afs$associated_trait))){

    site_list_i <- afs[which(afs$associated_trait == i)]
    phcol <- which(c(unique(afs$associated_trait)) ==i)

    for(m in 1:nrow(site_list_i)) {
        # filter afs data to each top site in turn
        ms <- site_list_i[m, ]
        # filter genotype data to that site & only genotype digits
        gs <- gt %>% filter(CHROM==ms$chr) %>% filter(POS ==ms$ps) %>% select(-c(CHROM, POS))

        # divide line IDs by genotype
        tgs <- t(gs)
        a1 <- names(tgs[grep("0", tgs),])
        a2 <- names(tgs[grep("2", tgs),])
        a1a2 <- names(tgs[grep("1", tgs),])

        a1t <- tibble("A1" = a1)
        a1f <- right_join(pheno, a1t, by=c("Genotype"="A1")) %>% select(all_of(i))
        a1f$genotype <- "A1 homozygous"

        a2t <- tibble("A2"=a2)
        a2f <- right_join(pheno, a2t, by=c("Genotype"="A2")) %>% select(all_of(i))
        a2f$genotype <- "A2 homozygous"

        jn <- bind_rows(a1f, a2f)

        # there may be no heterozygous genotypes
        # in which case, you can't join that data
        if(length(a1a2)!=0) {
                a1a2t <- tibble("A1A2"=a1a2)
                a1a2f <- right_join(pheno, a1a2t, by=c("Genotype"="A1A2")) %>% select(all_of(i))
                a1a2f$genotype <- "heterozygous"

                jn <- bind_rows(jn, a1a2f)

                # need to set the levels...
                jn$genotype <- factor(jn$genotype, levels=c("A1 homozygous", "heterozygous", "A2 homozygous")) 
        }


        ## definitley needs sample numbers, cus there's like 2 heterozygotes...
        # maybe even have a filter step for low numbers of heterozygotes
        sample_counts <- jn %>% count(genotype) 

        # difference in average A1 / A2 phenotype
        grp_avg <- jn %>% group_by(genotype) %>% summarise(avg = mean(get(i), na.rm=T))
        a1a2_diff <- grp_avg[1, 2] - grp_avg[nrow(grp_avg), 2]
        
        g <- ggplot(jn, aes(x=genotype, y= get(i))) + geom_boxplot() +  stat_n_text() + labs(y=tidy_text_substitution(i), title=paste(ms$chr, ms$ps), subtitle=paste0("difference in A1, A2 trait average: ", round(a1a2_diff, 3))) + theme_bw()
        ggsave(paste0("allele_effect/", ms$chr, "_", ms$ps, "_", i, ".png"), g)
    }
}