#!/usr/bin/env Rscript

#SBATCH --job-name=gwas
#SBATCH -o /rhome/jmarz001/bigdata/Ag-Competition/scripts/gwas/7b_plot_AFS.stdout
#SBATCH --ntasks=1
#SBATCH --mem=60gb
#SBATCH -t 02:00:00
#SBATCH -p koeniglab

library(pacman)
p_load(tidyverse, data.table, gridExtra, ggsci, ggpubr, Cairo)

setwd("/rhome/jmarz001/bigdata/Ag-Competition/results/gwas")
source("../../scripts/CUSTOM_FNS.R")

neutral_sampled <- fread("neutral_sites_sampled.tsv")
pop_freqs <- fread("gwas_sites_pop_freq_binned.tsv")
pop_freqs_neutral <- fread("neutral_sites_pop_freq_binned.tsv")

## make binned/count data into a fraction of the whole
# gwas sites
pop_freqs[,2:ncol(pop_freqs)] <- pop_freqs %>% reframe(across(-c(bins), \(x) x/sum(x)))
colSums(pop_freqs[,2:ncol(pop_freqs)])
# neutral sites sampled w same starting freq as gwas sites
neutral_sampled[,2:ncol(neutral_sampled)] <- neutral_sampled %>% reframe(across(-c(bins), \(x) x/sum(x)))
colSums(neutral_sampled[,2:ncol(neutral_sampled)])
# all neutral sites 
pop_freqs_neutral[,2:ncol(pop_freqs_neutral)] <- pop_freqs_neutral %>% reframe(across(-c(bins), \(x) x/sum(x)))
colSums(pop_freqs_neutral[,2:ncol(pop_freqs_neutral)])


#### plot allele frequency spectrum for

### all gwas sites
gwas_sites <- 
pop_freqs %>% 
  pivot_longer(cols=c(-bins), names_to="group") %>%
  ggplot(aes(x=bins, y=value, fill=group)) +
        geom_bar(stat='identity', position='dodge', width=0.7) +
        theme_classic(base_size=20) +
        theme(axis.text.x = element_text(angle=45, hjust=1)) +
        labs(x="Frequency Ranges", y="Proportion of Sites", title="Site Frequency Spectrum", subtitle="GWAS sites") +
        scale_fill_discrete(name = "Generation",
        labels = c("Parents", "F18", "F28", "F50", "F58"))
ggsave("gwas_sites_afs.png", gwas_sites, height=9, width=15, units="in")

### gwas sites for parents, F18, and F58 generations
pop_freqs2 <- pop_freqs %>% select(c(bins, F0_bins, F18_bins, F58_bins))
gwas_sites2 <- 
pop_freqs2 %>% 
  pivot_longer(cols=c(-bins), names_to="group") %>%
  ggplot(aes(x=bins, y=value, fill=group)) +
        geom_bar(stat='identity', position='dodge', width=0.7) +
        theme_classic(base_size=20) +
        theme(axis.text.x = element_text(angle=45, hjust=1)) +
        labs(x="Frequency Ranges", y="Proportion of Sites", title="Site Frequency Spectrum", subtitle="GWAS sites") +
        scale_fill_discrete(name = "Generation",
        labels = c("Parents", "F18", "F58"))
ggsave("gwas_sites_limitedgens_afs.png", gwas_sites2, height=9, width=15, units="in")



### genome-wide AF change over time
genome_sites <- 
pop_freqs_neutral %>% 
  pivot_longer(cols=c(-bins), names_to="group") %>%
  ggplot(aes(x=bins, y=value, fill=group)) +
        geom_bar(stat='identity', position='dodge', width=0.7) +
        theme_classic(base_size=20) +
        theme(axis.text.x = element_text(angle=45, hjust=1)) +
        labs(x="Frequency Ranges", y="Proportion of Sites", title="Site Frequency Spectrum", subtitle="all genome sites") +
        scale_fill_discrete(name = "Generation",
        labels = c("Parents", "F18", "F28", "F50", "F58"))
ggsave("genome_sites_afs.png", genome_sites, height=9, width=15, units="in")


### comparing gwas & genome sites
##  F0-gwas / F0-neutral sites
##  F18-gwas / F18-neutral sites
##  F58-gwas / F58-neutral sites

# join gwas & neutral data together
pop_freqs <- pop_freqs %>% select(-c(F28_bins, F50_bins))
pop_freqs_neutral <- pop_freqs_neutral %>% select(-c(F28_bins_neutral, F50_bins_neutral))

df <- full_join(pop_freqs, pop_freqs_neutral)
combo <- df %>% pivot_longer(cols=c(-bins), names_to="group")
combo$group <- as.factor(combo$group)
combo$group <- factor(combo$group,  # Change ordering manually
levels = c("F0_bins", "F18_bins", "F18_bins_neutral", "F58_bins", "F58_bins_neutral"))


#### bar plots comparing groups w their neutral match
f0 <- combo %>% 
    filter(group=="F0_bins" | group=="F0_bins_neutral") %>%
    ggplot(aes(x=bins, y=value, fill=group)) +
      geom_bar(stat='identity', position='dodge', width=0.7) +
      theme_classic(base_size=20) +
      theme(axis.text.x = element_text(angle=45, hjust=1)) +
      labs(x="Frequency Ranges", y="Proportion of Sites", title="F0 Site Frequency Spectrum", subtitle="Comparing GWAS and neutral sites allele frequencies") +
      scale_fill_discrete(name = "Site Groups", labels = c("GWAS", "Neutral"))

f1 <- combo %>% 
    filter(group=="F18_bins" | group=="F18_bins_neutral") %>%
    ggplot(aes(x=bins, y=value, fill=group)) +
      geom_bar(stat='identity', position='dodge', width=0.7) +
      theme_classic(base_size=20) +
      theme(axis.text.x = element_text(angle=45, hjust=1)) +
      labs(x="Frequency Ranges", y="Proportion of Sites", title="F18 Site Frequency Spectrum") +
      scale_fill_discrete(name = "Site Groups", labels = c("GWAS", "Neutral"))

f2 <- combo %>% 
    filter(group=="F58_bins" | group=="F58_bins_neutral") %>%
    ggplot(aes(x=bins, y=value, fill=group)) +
      geom_bar(stat='identity', position='dodge', width=0.7) +
      theme_classic(base_size=20) +
      theme(axis.text.x = element_text(angle=45, hjust=1)) +
      labs(x="Frequency Ranges", y="Proportion of Sites", title="F58 Site Frequency Spectrum") +
      scale_fill_discrete(name = "Site Groups", labels = c("GWAS", "Neutral"))

png("P_F18_F58_gwas_neutral_freq_comparisons.png", width=10, height=21, units="in", res=300)
ggarrange(f0, f1, f2, ncol=1)
dev.off()


## all gwas and neutral generations (P, F0, F18)
combo_afs <- combo %>%
  ggplot(aes(x=bins, y=value, fill=group)) +
        geom_bar(stat='identity', position='dodge', width=0.7) +
        theme_classic(base_size=20) +
        theme(axis.text.x = element_text(angle=45, hjust=1)) +
        labs(x="Frequency Ranges", y="Proportion of Sites", title="Site Frequency Spectrum") +
        scale_fill_discrete(name = "Generation & Site Association",
        labels = c("Parents - GWAS", "F18 - GWAS", "F18 - Neutral", "F58 - GWAS", "F58 - Neutral"))
ggsave("gwas_and_genome_sites_afs.png", combo_afs, height=9, width=15, units="in")
