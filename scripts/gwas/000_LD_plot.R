#!/usr/bin/env Rscript
#SBATCH -o /rhome/jmarz001/bigdata/Ag-Competition/scripts/gwas/000_LD_plot.stdout
#SBATCH --ntasks=1
#SBATCH --mem=100gb
#SBATCH -t 01:30:00
#SBATCH -p short

library(tidyverse)
library(data.table)

setwd("/rhome/jmarz001/bigdata/Ag-Competition/results/gwas")
source("../../scripts/CUSTOM_FNS.R")

ld <- fread("LD_10kbwin.ld.gz")
sites <- fread("top_sites.txt")$SNP

# format LD files SNP cols to match site list
ld$SNP_A <- gsub("(chr\\dH):(\\d+),", "\\1_\\2", ld$SNP_A)
ld$SNP_B <- gsub("(chr\\dH):(\\d+),", "\\1_\\2", ld$SNP_B)

# loop over each top site and plot the LD
for(i in 1:length(sites)){

    x <- sites[i]

    # filter ld file to the site
    ld_filt <- ld %>% filter(SNP_A %in% x)

    # print average LD over chromosome
    print(mean(ld_filt$R2))
    print(quantile(ld_filt$R2))

    # make heatmap plot
    tit <- paste0(x, " LD")
    g <- ggplot(ld_filt, aes(x=SNP_A, y=SNP_B, fill=R2)) + geom_tile() + labs(title=tit) + 
    scale_fill_gradient2(low="lightblue", high="darkblue", 
        midpoint=0.5, 
        breaks=seq(0,1,0.25), #breaks in the scale bar
        limits=c(0, 1))

    # write out
    ggsave(paste0("LD_heatmap_", i, ".png"), g)

    # make distribution of LD values?
    h <- ggplot(ld_filt, aes(R2)) + geom_histogram()
    ggsave(paste0("LD_hist_", i, ".png"), h)
}




## filter to top sites all at once
ld_filt <- ld %>% filter(SNP_A %in% sites)

print(mean(ld_filt$R2))
print(quantile(ld_filt$R2))

m <- ggplot(ld_filt, aes(x=SNP_A, y=SNP_B, fill=R2)) + 
    geom_tile() + 
    labs(title=tit) + 
    scale_fill_gradient2(low="lightblue", high="darkblue", 
        midpoint=0.5, 
        breaks=seq(0,1,0.25), #breaks in the scale bar
        limits=c(0, 1))

ggsave("LD_heatmap_all.png", m)
