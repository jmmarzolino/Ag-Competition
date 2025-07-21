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

## read in data sets

# separate phenotype cols?
# run one association file results & phenotype file at a time
#ie. phenotype ft = association 6...
file_traits <- fread("association_files_traits.txt")
afs <- fread("progeny_AF_change.tsv")
# line IDs + genotype
gt <- fread("out.GT.FORMAT")
# line ID + phenotype
ph <- fread("all_traits.fam")

# filter genotype file by line IDs in phenotype file
#colnames(gt_filt)[3:ncol(gt_filt)]
lines <- ph$V1
gt <- gt %>% select(c("CHROM", "POS", any_of(lines)))



for (i in 6:8) {

    traittxt <- file_traits[which(i == file_traits$trait_num), 1][[1]]

    ## read in list of gwas significant & significantly changed allele count sites
    ## these will be the only sites you need
        # filtering file: 'top' snps
    sites <- fread(paste0("ASSOC_", i, "_lmm.assoc.clumped"))[,1:11]

    ### start by filtering allele allele freqs to only those in sites file
    afs_i <- right_join(afs, sites, by=c("CHR"="CHR", "BP"="BP"))
    # cool, great, 7700 sites is much better to deal w for testing especially...

    # format that data...
    # filter genotype data to common sites
    gt_filt <- right_join(gt, afs_i, by=c("CHROM"="CHR", "POS"="BP")) 

    # ph filter to V6 - V11 colnames that match number in file numbers
    phcol <- paste0("V", i)

    #pdf(paste0("top_sites_allele_effects_", i, ".pdf"))


    for(m in 1:nrow(afs_i)) {
        ms <- afs_i[m, ]
        gs <- gt_filt[m, ]


        # divide genotype subset file line IDs by genotype
        tgs <- t(gs)
        a1 <- names(tgs[grep("0\\|0", tgs),])
        a2 <- names(tgs[grep("1\\|1", tgs),])
        a1a2 <- c(names(tgs[grep("1\\|0", tgs),]), names(tgs[grep("0\\|1", tgs),]))


        a1t <- tibble("A1" = a1)
        a1f <- right_join(ph, a1t, by=c("V1"="A1")) %>% select(all_of(phcol))
        a1f$genotype <- "A1 homozygous"

        a2t <- tibble("A2"=a2)
        a2f <- right_join(ph, a2t, by=c("V1"="A2")) %>% select(all_of(phcol))
        a2f$genotype <- "A2 homozygous"

        jn <- bind_rows(a1f, a2f)

        # there may be no heterozygous genotypes
        # in which case, you can't join that data
        if(length(a1a2)!=0) {
                a1a2t <- tibble("A1A2"=a1a2)
                a1a2f <- right_join(ph, a1a2t, by=c("V1"="A1A2")) %>% select(all_of(phcol))
                a1a2f$genotype <- "heterozygous"

                jn <- bind_rows(jn, a1a2f)

                # need to set the levels...
                jn$genotype <- factor(jn2$genotype, levels=c("A1 homozygous", "heterozygous", "A2 homozygous")) 
        }


        ## definitley needs sample numbers, cus there's like 2 heterozygotes...
        # maybe even have a filter step for low numbers of heterozygotes
        sample_counts <- jn %>% count(genotype) 

        g <- ggplot(jn, aes(x=genotype, y= get(phcol))) + geom_boxplot() +  stat_n_text() + labs(y=tidy_text_substitution(traittxt), subtitle=paste0("site past ", ms$top, " threshhold"), title=paste(ms$CHR, ms$BP)) + theme_bw()
        ggsave(paste0("allele_effect/", "Rplot_trait", i, "_row", m, ".png"), g)
    }

    #dev.off()
}
