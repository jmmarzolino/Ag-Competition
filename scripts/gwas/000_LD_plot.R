#!/usr/bin/env Rscript
#SBATCH -o /rhome/jmarz001/bigdata/Ag-Competition/scripts/gwas/000_LD_plot.stdout
#SBATCH --ntasks=1
#SBATCH --mem=100gb
#SBATCH -t 01:30:00
#SBATCH -p short

library(tidyverse)
library(gridExtra)
library(data.table)

setwd("/rhome/jmarz001/bigdata/Ag-Competition/results/gwas")
source("../../scripts/CUSTOM_FNS.R")
# manhattan plots, zoomed in on specific peak, and points colored by LD with peak site


# load LD measures
ld <- fread("LD_10kbwin.ld.gz")

# load top sites list & their gwas info
sites <- fread("top_sites.txt")$SNP
sites2 <- fread("gwas_top_sites.tsv")
#setdiff(sites, sites2$rs) # check for diffs in list of clump sites

# read in file connecting file name/number and trait
files <- fread("association_files_traits.txt")
files$trait_names <- tidy_text_substitution(files$trait_names)


# loop over each top site and plot the LD
addPlot <- function(site){

    # filter LD file to R2 value of each site vs focal site
    tmp <- ld %>% filter(SNP_A == site)

    # now that all data is related to that focal site, remove the "A" site columns
    tmp <- tmp %>% select(-c("CHR_A", "BP_A", "SNP_A"))
    # and replace col names
    colnames(tmp) <- gsub("(\\w+)_B", "\\1", colnames(tmp))

    # pull the trait associated w the top site
    gwas_data <- sites2 %>% filter(site == rs)
    gwas_file <- unlist(files[which(files$trait_names == gwas_data$associated_trait), 4])
    full_gwas <- fread(gwas_file)

    # filter full gwas file to sites +/- 500kb of top site
    site_chr <- gsub("(chr\\dH):\\d+", "\\1", site)
    site_bp <- as.numeric(gsub("chr\\dH:(\\d+)", "\\1", site))
    max_bp <- site_bp + 500000
    min_bp <- site_bp - 500000

    site_region_gwas <- full_gwas %>% filter(chr == site_chr & ps <= max_bp & ps >= min_bp)

    ## combine R2 values w gwas sites
    jnt <- full_join(site_region_gwas, tmp, by=c("chr"="CHR", "rs"="SNP", "ps"="BP"))

    ## manhattan plot
  # following code adapted from:
    # https://www.r-graph-gallery.com/wp-content/uploads/2018/02/Manhattan_plot_in_R.html
    # format for plotting
    result <- jnt %>%
      # Compute chromosome size
      group_by(chr) %>%
      summarise(chr_len = max(ps)) %>%
      # Calculate cumulative position of each chromosome
      mutate(tot = cumsum(as.numeric(chr_len))-as.numeric(chr_len)) %>%
      select(-chr_len) %>%
      # Add this info to the initial dataset
      left_join(jnt, ., by=c("chr" = "chr")) %>%
      # Add a cumulative position of each SNP
      arrange(chr, ps) %>%
      mutate(BPcum = ps+tot)

    axisdf <- result %>% group_by(chr) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

    threshold <- 0.05/nrow(full_gwas)

    # mark the central snp...
    centralsnp <- result[which(result$rs == site), BPcum]

    # Manhattan plot
    ggplot(result, aes(x=BPcum, y=-log10(p_lrt))) +
        # Show all points
        geom_point(aes(color=R2), size=1.3) +
        geom_hline(aes(yintercept=-log10(threshold)), color = "firebrick1", linetype="dashed", alpha=0.7) +
        #geom_vline(aes(xintercept=centralsnp), color="dodgerblue", linetype="dashed", alpha=0.6) +
        # custom X axis:
        scale_x_continuous(label = c((centralsnp-500000), centralsnp, (centralsnp+500000)), breaks= c((centralsnp-500000), centralsnp, (centralsnp+500000))) +
        scale_y_continuous(expand = c(0, 0.5)) +     # remove space between plot area and x axis
        # Custom theme:
        theme_classic() +
        theme(panel.border = element_blank(),
              panel.grid.major.x = element_blank(),
              panel.grid.minor.x = element_blank(),
              text=element_text(size=16),
              plot.title = element_text(size=12),
              plot.subtitle = element_text(size=10)) +
        xlab(gsub("chr(\\d)H", "Chromosome \\1", axisdf$chr)) +
        ylab(expression(-log[10](italic(p))))
}

lst <- sites
test <- lapply(lst, addPlot)

pdf("LD_decay_manhattans.pdf")
marrangeGrob(grobs=test, nrow=2, ncol=1)
dev.off()