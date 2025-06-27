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

#neu_f0 <- fread("neutral_sites_sampled_F0.tsv")
#neu_f18 <- fread("neutral_sites_sampled_F18.tsv")
pop_freqs <- fread("gwas_sites_pop_freq_binned.tsv")
pop_freqs_neutral <- fread("neutral_sites_pop_freq_binned.tsv")

##################################################
## keep this section as notes on data format...
# sites associated w traits
#tmp <- fread("../../../IPK_Analysis/results/GWAS/POP_AF/NeutralTraits_sites.tsv") %>% select(-c(chr, rs, ps))
#tmp2 <- fread("../../../IPK_Analysis/results/GWAS/POP_AF/AssocTraits_neutral_sites.tsv")
#assoc_trait_sites$delta <- assoc_trait_sites$F58_AF - assoc_trait_sites$F18_AF
#assoc_trait_sites$group <- "trait associated"
##################################################


## make binned/count data into a fraction of the whole
# gwas sites
pop_freqs[,2:4] <- pop_freqs %>% reframe(across(c(F0_bins, F18_bins, F58_bins), \(x) x/sum(x)))

# all neutral sites 
pop_freqs_neutral[,2:3] <- pop_freqs_neutral %>% reframe(across(-bins, \(x) x/sum(x)))




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
ggsave("gwas_sites_afs.png", gwas_sites)

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
ggsave("gwas_sites_limitedgens_afs.png", gwas_sites2)



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
        labels = c("Parents", "F18", "F58"))
ggsave("genome_sites_afs.png", genome_sites)


### comparing gwas & genome sites
##  F0-gwas / F0-neutral sites
##  F18-gwas / F18-neutral sites
##  F58-gwas / F58-neutral sites

# join gwas & neutral data together
df <- full_join(pop_freqs, pop_freqs_neutral)
combo <- df %>% pivot_longer(cols=c(-bins), names_to="group")
combo$group <- as.factor(combo$group)
combo$group <- factor(combo$group,  # Change ordering manually
levels = c("F0_bins", "F18_bins", "F58_bins", "F0_bins_neutral", "F18_bins_neutral", "F58_bins_neutral"))

combo_afs <- combo %>%
  ggplot(aes(x=bins, y=value, fill=group)) +
        geom_bar(stat='identity', position='dodge', width=0.7) +
        theme_classic(base_size=20) +
        theme(axis.text.x = element_text(angle=45, hjust=1)) +
        labs(x="Frequency Ranges", y="Proportion of Sites", title="Site Frequency Spectrum")#, subtitle="all genome sites") +
        scale_fill_discrete(name = "Generation",
        labels = c("Parents", "F18", "F58"))
ggsave("gwas_and_genome_sites_afs.png", combo_afs)







#### bar plots comparing groups w their neutral match
t_sig18 <- df4 %>% filter(group=="sig18" | group=="signeu18")
t_t18 <- df4 %>% filter(group=="trait18" | group=="trneu18")
t_sig58 <- df4 %>% filter(group=="sig58" | group=="signeu58")
t_t58 <- df4 %>% filter(group=="trait58" | group=="trneu58")


f1 <- t_sig18 %>%
    ggplot(aes(x=bins, y=value, fill=group)) +
      geom_bar(stat='identity', position='dodge', width=0.7) +
      theme_classic(base_size=20) +
      theme(axis.text.x = element_text(angle=45, hjust=1)) +
      labs(x="Frequency Ranges", y="Proportion of Sites", title="F18 Site Frequency Spectrum") +
      scale_fill_discrete(name = "Site Groups", labels = c("Significant Traits", "Neutral"))


f2 <- t_t18 %>%
    ggplot(aes(x=bins, y=value, fill=group)) +
      geom_bar(stat='identity', position='dodge', width=0.7) +
      theme_classic(base_size=20) +
      theme(axis.text.x = element_text(angle=45, hjust=1)) +
      labs(x="Frequency Ranges", y="Proportion of Sites", title="F18 Site Frequency Spectrum") +
      scale_fill_discrete(name = "Site Groups", labels = c("Trait Associated", "Neutral"))


f3 <- t_sig58 %>%
    ggplot(aes(x=bins, y=value, fill=group)) +
      geom_bar(stat='identity', position='dodge', width=0.7) +
      theme_classic(base_size=20) +
      theme(axis.text.x = element_text(angle=45, hjust=1)) +
      labs(x="Frequency Ranges", y="Proportion of Sites", title="F58 Site Frequency Spectrum") +
      scale_fill_discrete(name = "Site Groups", labels = c("Significant Traits", "Neutral"))

f4 <- t_t58 %>%
    ggplot(aes(x=bins, y=value, fill=group)) +
      geom_bar(stat='identity', position='dodge', width=0.7) +
      theme_classic(base_size=20) +
      theme(axis.text.x = element_text(angle=45, hjust=1)) +
      labs(x="Frequency Ranges", y="Proportion of Sites", title="F58 Site Frequency Spectrum") +
      scale_fill_discrete(name = "Site Groups", labels = c("Trait Associated", "Neutral"))

png("SFS_2GroupComparisons.png", width=20, height=12, units="in", res=300)
ggarrange(f1, f2, f3, f4)
dev.off()


