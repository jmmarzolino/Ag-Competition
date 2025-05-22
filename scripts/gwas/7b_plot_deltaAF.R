#!/usr/bin/env Rscript

#SBATCH --job-name=gwas
#SBATCH -o /rhome/jmarz001/bigdata/Ag-Competition/scripts/gwas/7b_plot_deltaAF.stdout
#SBATCH --ntasks=1
#SBATCH --mem=80gb
#SBATCH -t 10:00:00
#SBATCH -p koeniglab

library(pacman)
p_load(tidyverse, data.table, gridExtra, ggsci, ggpubr, Cairo)

setwd("/rhome/jmarz001/bigdata/Ag-Competition/results/gwas/POP_AF")
source("../../../scripts/standardize_trait_text_function.R")

# sites associated w sig traits = group 1
sig <- fread("SigTraits_sites.tsv") %>% select(-c(chr, rs, ps))
sig$delta_AF <- sig$F58_AF - sig$F18_AF
sig_neu <- fread("SigTraits_neutral_sites.tsv")

# sites associated w neutral traits = group 2
traits <- fread("NeutralTraits_sites.tsv") %>% select(-c(chr, rs, ps))
traits$delta_AF <- traits$F58_AF - traits$F18_AF
traits_neu <- fread("AssocTraits_neutral_sites.tsv")



x <- tibble(freq_change = sig$delta_AF, group="sig_trait_sites", neutral=F)
y <- tibble(freq_change = sig_neu$delta_AF_neutral, group="sig_trait_sites", neutral=T)
w <- tibble(freq_change = traits$delta_AF, group="trait_assoc_sites", neutral=F)
z <- tibble(freq_change = traits_neu$delta_AF_neutral, group="trait_assoc_sites", neutral=T)
siggs <- bind_rows(x, y)
figgs <- bind_rows(w, z)

all <- bind_rows(siggs, figgs)

all$neu2 <- gsub("TRUE", "_neutral", all$neutral)
all$neu2 <- gsub("FALSE", "", all$neu2)

all$combogroup <- paste(all$group, all$neu2, sep="")
all$group2 <- factor(all$group, levels=c("sig_trait"))


# what is the mean, sd, and distribution of the three groups? that'll impact which statistic you use to determine if the 3 groups are from the same distribution or not
all %>% group_by(group,neutral) %>% summarise('average change' = mean(freq_change), 'median change' = median(freq_change), 'standard dev'=sd(freq_change), 'max'=max(freq_change), 'min'=min(freq_change))





### do samples come from the same or different distributions?
# trait-SNP deltaAF vs neutral-SNP deltaAF
# use absolute value of delta AF so the data should then be parametric instead of bimodal?


# check for homogenous variance of samples (but I think it's false...)
boxplot(sig$delta_AF, sig_neu$delta_AF_neutral)
var(sig$delta_AF) /var(sig_neu$delta_AF_neutral) #if ratio of var is less than 4, assume equal var
var(abs(sig$delta_AF)) /var(abs(sig_neu$delta_AF_neutral))
# both variance comparisons are huge, so groups don't have equal var and the welch t-test should be used
var(traits$delta_AF)/var(traits_neu$delta_AF_neutral)
var(abs(traits$delta_AF))/var(abs(traits_neu$delta_AF_neutral))
# variance ratios are huge for both groups, so go ahead w welch test

# just add one more verification of variance differences
#var.test(x, y, alternative = “two.sided”)
var.test(sig$delta_AF, sig_neu$delta_AF_neutral, alternative = "two.sided")
var.test(abs(sig$delta_AF), abs(sig_neu$delta_AF_neutral), alternative = "greater")

var.test(traits$delta_AF, traits_neu$delta_AF_neutral, alternative = "two.sided")
var.test(abs(traits$delta_AF), abs(traits_neu$delta_AF_neutral), alternative = "greater") #ie. testing if var of deltaAF GREATER THAN deltaAF_neutral

# variance between groups is definitely NOT equal
library(moments)
skewness(traits$delta_AF)
skewness(abs(traits$delta_AF))

# meaning use the welch t-test
# welch t-test
#t.test(g1, g2, alternative = c(“two.sided”, “less”, “greater”))
# t-test
#t.test(g1, g2, var.equal=T)

#t.test(sig$delta_AF, sig_neu$delta_AF_neutral, var.equal=F, alternative="two.sided")
#t.test(abs(sig$delta_AF), abs(sig_neu$delta_AF_neutral), var.equal=F, alternative="two.sided")





#tests don't seem right... use kruskal test?
# non parametric, one treatment variable?, two or more categories
# for this question, it's one variable (snp source of delta AF) w 2 categories (associated sites vs neutral)
# I think
kruskal.test(list(sig$delta_AF, sig_neu$delta_AF_neutral))
print("kruskal-wallis test p-value for: significant traits vs neutral sites")
print(kruskal.test(list(sig$delta_AF, sig_neu$delta_AF_neutral))$p.value)
kruskal.test(list(abs(sig$delta_AF), abs(sig_neu$delta_AF_neutral)))


kruskal.test(list(traits$delta_AF, traits_neu$delta_AF_neutral))
print("kruskal-wallis test p-value for: -- traits vs neutral sites")
print(kruskal.test(list(traits$delta_AF, traits_neu$delta_AF_neutral))$p.value)
kruskal.test(list(abs(traits$delta_AF), abs(traits_neu$delta_AF_neutral)))





# plot distribution of direction and magnitude of allele freq change over generations
siggs$neutral <- factor(siggs$neutral, levels=c("TRUE", "FALSE"))
figgs$neutral <- factor(figgs$neutral, levels=c("TRUE", "FALSE"))

## plot distribution of allele freq changes
g <- ggplot(siggs) +
        geom_histogram(aes(freq_change, fill=neutral), binwidth=0.05, alpha=0.7) +
        xlim(c(-1,1)) +
        labs(x="frequency change", title="Allele Frequency Change", subtitle="Sites Associated with Significantly Changed Traits") +
        scale_fill_discrete(name="", labels=c("trait associated sites", "neutral sample sites")) +
        theme_bw(base_size=18) +
        theme(legend.position = "top")

h <- ggplot(figgs) +
        geom_histogram(aes(freq_change, fill=neutral), binwidth=0.05, alpha=0.7) +
        xlim(c(-1,1)) +
        labs(x="frequency change", title="Allele Frequency Change", subtitle="Sites Associated with Traits") +
        scale_fill_discrete(name="", labels=c("trait associated sites", "neutral sample sites")) +
        theme_bw(base_size=18) +
        theme(legend.position = "top")

png("deltaAF_hist_AssocSitesVNeutral.png", width=16, height=12, units="in", res=300)
ggarrange(g, h)
dev.off()

## plot distribution & median, IQR, range of AF change
### comparing selected sites w their neutral site match
p <- ggplot(siggs, aes(x=neutral, y=freq_change)) +
        geom_violin() +
        geom_boxplot(width=0.1, alpha=0.8) +
        labs(title="Allele Frequency Change", subtitle="Sites Associated with Significantly Changed Traits", y="Frequency Change") +
        #scale_fill_discrete(name="", labels=c("trait associated sites", "neutral sample sites")) +
        scale_x_discrete(name="", labels=c("trait associated sites", "neutral sample sites")) +
        theme_bw(base_size=18) +
        theme(axis.text.x = element_text(angle=45, hjust=1))

q <-  ggplot(figgs, aes(x=neutral, y=freq_change)) +
        geom_violin() +
        geom_boxplot(width=0.1, alpha=0.8) +
        labs(title="Allele Frequency Change", subtitle="Sites Associated with Traits", y="Frequency Change") +
        #scale_fill_discrete(name="", labels=c("trait associated sites", "neutral sample sites")) +
        scale_x_discrete(name="", labels=c("trait associated sites", "neutral sample sites")) +
        theme_bw(base_size=18) +
        theme(axis.text.x = element_text(angle=45, hjust=1))
png("deltaAF_boxplot_AssocSitesVNeutral.png", width=14, height=10, units="in", res=300)
ggarrange(p, q)
dev.off()


p <- ggplot(all, aes(y=freq_change, x=combogroup)) +
      geom_violin(aes(fill=combogroup)) +
      guides(fill='none') +
      geom_boxplot(width=0.1, alpha=0.7) +
      theme_bw(base_size=16) +
      theme(axis.text.x = element_text(angle=45, hjust=1)) +
      #scale_x_discrete(label = labels) +
      scale_x_discrete(name = "Site Groups",
      labels = c("Significant Traits", "Significant-Traits Neutral", "Trait Associated", "Trait-Associated Neutral")) +
      labs(title="Allele Frequency Change", y="Frequency Difference", x="")

png("deltaAF_boxplot_allgroups.png", width=6, height=8, units="in", res=300)
p
dev.off()
