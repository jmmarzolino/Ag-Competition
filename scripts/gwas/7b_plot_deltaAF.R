#!/usr/bin/env Rscript
#SBATCH --job-name=gwas
#SBATCH -o /rhome/jmarz001/bigdata/Ag-Competition/scripts/gwas/7b_plot_deltaAF.stdout
#SBATCH --ntasks=1
#SBATCH --mem=60gb
#SBATCH -t 02:00:00
#SBATCH -p koeniglab

library(pacman)
p_load(tidyverse, data.table, gridExtra, ggsci, ggpubr, Cairo)

setwd("/rhome/jmarz001/bigdata/Ag-Competition/results/gwas")
source("../../scripts/CUSTOM_FNS.R")


## siggs sites, siggs generations, siggsele frequencies
## read in siggs sites' siggsele frequencies
df <- fread("siggs_sites_siggsele_freqs.tsv")
## should already be filtered to the right columns
df$snp <- paste0(df$CHR, "_", df$BP)

## polarize siggs sites to the BENEFICIAL siggsele (F18 to F58 rather than starting w F0...)
# which is the one which increases freq over gens
# do that for siggs sites
df$beneficial_siggsele <- "A1"

x <- round(df$F58_AF - df$F18_AF, 4)
dfa1 <- df[which(x >= 0), ]
dfa2 <- df[which(x < 0), ]
# sites w beneficial siggsele 2...
dfa2$beneficial_siggsele <- "A2"
dfa2[, 5:9] <- 1 - (dfa2[, 5:9])

# stitch data sets back together
df2 <- bind_rows(dfa1, dfa2)

## then pull sites associated w traits
## and filter them to their own dataset
# load associated sites list
siggs_sig <- fread("siggs_gwas_sig_sites.tsv")
siggs_sig$rs <- paste0(siggs_sig$chr, "_", siggs_sig$ps)

trait_assoc_sites <- df2 %>% filter(snp %in% siggs_sig$rs)
neutral_sites <- df2 %>% filter(!(snp %in% siggs_sig$rs))

# calculate change in siggsele freq between 0-18 and 18-58 for later...
trait_assoc_sites$F0_F18_AF_delta <- trait_assoc_sites$F18_AF - trait_assoc_sites$F0_AF 
trait_assoc_sites$F18_F58_AF_delta <- trait_assoc_sites$F58_AF - trait_assoc_sites$F18_AF


neutral_sampled <- fread("neutral_sites_sampled.tsv")


x <- tibble(freq_change = sig$delta_AF, group="sig_trait_sites", neutral=F)
y <- tibble(freq_change = sig_neu$delta_AF_neutral, group="sig_trait_sites", neutral=T)

siggs <- bind_rows(x, y)

siggs$neu2 <- gsub("TRUE", "_neutral", siggs$neutral)
siggs$neu2 <- gsub("FALSE", "", siggs$neu2)

siggs$combogroup <- paste(siggs$group, siggs$neu2, sep="")
siggs$group2 <- factor(siggs$group, levels=c("sig_trait"))


# what is the mean, sd, and distribution of the three groups? that'll impact which statistic you use to determine if the 3 groups are from the same distribution or not
siggs %>% group_by(group,neutral) %>% summarise('average change' = mean(freq_change), 'median change' = median(freq_change), 'standard dev'=sd(freq_change), 'max'=max(freq_change), 'min'=min(freq_change))





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
print("kruskal-wsiggsis test p-value for: significant traits vs neutral sites")
print(kruskal.test(list(sig$delta_AF, sig_neu$delta_AF_neutral))$p.value)
kruskal.test(list(abs(sig$delta_AF), abs(sig_neu$delta_AF_neutral)))


kruskal.test(list(traits$delta_AF, traits_neu$delta_AF_neutral))
print("kruskal-wsiggsis test p-value for: -- traits vs neutral sites")
print(kruskal.test(list(traits$delta_AF, traits_neu$delta_AF_neutral))$p.value)
kruskal.test(list(abs(traits$delta_AF), abs(traits_neu$delta_AF_neutral)))





# plot distribution of direction and magnitude of siggsele freq change over generations
siggs$neutral <- factor(siggs$neutral, levels=c("TRUE", "FALSE"))

## plot distribution of siggsele freq changes
g <- ggplot(siggs) +
        geom_histogram(aes(freq_change, fill=neutral), binwidth=0.05, alpha=0.7) +
        xlim(c(-1,1)) +
        labs(x="frequency change", title="Allele Frequency Change", subtitle="Sites Associated with Significantly Changed Traits") +
        scale_fill_discrete(name="", labels=c("trait associated sites", "neutral sample sites")) +
        theme_bw(base_size=18) +
        theme(legend.position = "top")


png("deltaAF_hist_AssocSitesVNeutral.png", width=16, height=12, units="in", res=300)
g
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

png("deltaAF_boxplot_AssocSitesVNeutral.png", width=14, height=10, units="in", res=300)
p
dev.off()


p <- ggplot(siggs, aes(y=freq_change, x=combogroup)) +
      geom_violin(aes(fill=combogroup)) +
      guides(fill='none') +
      geom_boxplot(width=0.1, alpha=0.7) +
      theme_bw(base_size=16) +
      theme(axis.text.x = element_text(angle=45, hjust=1)) +
      #scale_x_discrete(label = labels) +
      scale_x_discrete(name = "Site Groups",
      labels = c("Significant Traits", "Significant-Traits Neutral", "Trait Associated", "Trait-Associated Neutral")) +
      labs(title="Allele Frequency Change", y="Frequency Difference", x="")

png("deltaAF_boxplot_siggsgroups.png", width=6, height=8, units="in", res=300)
p
dev.off()
