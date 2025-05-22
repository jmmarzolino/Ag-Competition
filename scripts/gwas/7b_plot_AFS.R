#!/usr/bin/env Rscript

#SBATCH --job-name=gwas
#SBATCH -o /rhome/jmarz001/bigdata/Ag-Competition/scripts/gwas/7b_plot_AFS.stdout
#SBATCH --ntasks=1
#SBATCH --mem=80gb
#SBATCH -t 10:00:00
#SBATCH -p koeniglab

library(pacman)
p_load(tidyverse, data.table, gridExtra, ggsci, ggpubr, Cairo)

setwd("/rhome/jmarz001/bigdata/Ag-Competition/results/gwas/POP_AF")
source("../../../scripts/standardize_trait_text_function.R")

# sites associated w sig traits = group 1
sig_trait_sites <- fread("SigTraits_sites.tsv") %>% select(-c(chr, rs, ps))
#sig_trait_sites$delta <- sig_trait_sites$F58_AF - sig_trait_sites$F18_AF
#sig_trait_sites$group <- "significant trait associated"


# sites associated w neutral traits = group 2
assoc_trait_sites <- fread("NeutralTraits_sites.tsv") %>% select(-c(chr, rs, ps))
#assoc_trait_sites$delta <- assoc_trait_sites$F58_AF - assoc_trait_sites$F18_AF
#assoc_trait_sites$group <- "trait associated"



neu1 <- fread("SigTraits_neutral_sites.tsv")
neu2 <- fread("AssocTraits_neutral_sites.tsv")
#neu1$group <- "neutral 1"
#neu2$group <- "neutral 2"

##   join data
#neutrals <- full_join(neu1, neu2)
#df <- full_join(assoc_trait_sites, sig_trait_sites)
#df2 <- full_join(df, neutrals, by=c('F18_AF'='F18_AF_neutral', 'F58_AF'='F58_AF_neutral', 'delta'='delta_AF_neutral', 'group'))




# allele frequency bins for groups 1 & 2
#### Site Frequency Binning
# define other necessary vars, the breaks & ranges
bins = c("0.01-0.05", "0.05-0.1", '0.1-0.15', '0.15-0.2', '0.2-0.25', '0.25-0.3', '0.3-0.35', '0.35-0.4', '0.4-0.45', '0.45-0.5', '0.5-0.55', '0.55-0.6', '0.6-0.65', '0.65-0.7', '0.7-0.75', '0.75-0.8', '0.8-0.85', '0.85-0.9', '0.9-0.95', '0.95-1') #, 'unknown')

# binning function
binning <- function(x) {
    bins = c(0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0)

    a = 0
    b = 0
    c = 0
    d = 0
    e = 0
    f = 0
    g = 0
    h = 0
    i = 0
    j = 0
    k = 0
    l = 0
    m = 0
    n = 0
    o = 0
    p = 0
    q = 0
    r = 0
    s = 0
    t = 0
    u = 0
    zero = 0

    for(row in 1:nrow(x)) {
        # take in the frequency, put it in a bin
        Line = x[row, 1]

        #if(Line == bins[1]) {
        #  zero = zero + 1
        #} else if(Line >= bins[1] & Line < bins[2]) {
        if(Line >= bins[1] & Line < bins[2]) {
          a = a +1
        } else if(Line >= bins[2] & Line < bins[3]) {
          b = b +1
        } else if(Line >= bins[3] & Line < bins[4]) {
          c = c +1
        } else if(Line >= bins[4] & Line < bins[5]) {
          d = d +1
        } else if(Line >= bins[5] & Line < bins[6]) {
          e = e +1
        } else if(Line >= bins[6] & Line < bins[7]) {
          f = f +1
        } else if(Line >= bins[7] & Line < bins[8]) {
          g = g +1
        } else if(Line >= bins[8] & Line < bins[9]) {
          h = h +1
        } else if(Line >= bins[9] & Line < bins[10]) {
          i = i +1
        } else if(Line >= bins[10] & Line < bins[11]) {
          j = j +1
        } else if(Line >= bins[11] & Line < bins[12]) {
          k = k +1
        } else if(Line >= bins[12] & Line < bins[13]) {
          l = l +1
        } else if(Line >= bins[13] & Line < bins[14]) {
          m = m +1
        } else if(Line >= bins[14] & Line < bins[15]) {
          n = n +1
        } else if(Line >= bins[15] & Line < bins[16]) {
          o = o +1
        } else if(Line >= bins[16] & Line < bins[17]) {
          p = p +1
        } else if(Line >= bins[17] & Line < bins[18]) {
          q = q +1
        } else if(Line >= bins[18] & Line < bins[19]) {
          r = r +1
        } else if(Line >= bins[19] & Line < bins[20]) {
          s = s +1
        } else if(Line >= bins[20] & Line <= bins[21]) {
          t = t +1
        } else {
          u = u+1
        }
  }
    # for frequency divide bins by lines
    hist <- c(a, b, c, d, e, f, g, h, i, j, k, l, m, n, o, p, q, r, s, t) #, u) #, zero)
    hist <- hist/dim(x)[1]
    return(hist)
}


# apply plotting function and make list of those plots
a1 <- binning(assoc_trait_sites[,1])
a2 <- binning(assoc_trait_sites[,2])

s1 <- binning(sig_trait_sites[,1])
s2 <- binning(sig_trait_sites[,2])

n11 <- binning(neu1[,1])
n12 <- binning(neu1[,2])
n21 <- binning(neu2[,1])
n22 <- binning(neu2[,2])


# join data
batch <- tibble(bins, 'sig18'=a1, 'sig58'=a2, 'trait18'=s1, 'trait58'=s2, 'signeu18'=n11, 'signeu58'=n12, 'trneu18'=n21, 'trneu58'=n22)


##   format data for plotting
df4 <- batch %>% pivot_longer(cols=ends_with("8"), names_to='group')
df4$group <- as.factor(df4$group)
df4$group <- factor(df4$group,  # Change ordering manually
levels = c("sig18", "signeu18", "trait18", "trneu18", "sig58", "signeu58", "trait58", "trneu58"))

##   plot allele frequency spectrum
# bar plot comparing site freqs of neutral and signif sites
# custom colors
#custom_5col <- c("#C4961A", "#D16103", "#C3D7A4", "#52854C", "#00AFBB")
png("SFS.png", width=20, height=8, units="in", res=300)
ggplot(df4, aes(x=bins, y=value, fill=group)) +
      geom_bar(stat='identity', position='dodge', width=0.7) +
      theme_classic(base_size=20) +
      theme(axis.text.x = element_text(angle=45, hjust=1)) +
      labs(x="Frequency Ranges", y="Proportion of Sites", title="Site Frequency Spectrum") +
      scale_fill_discrete(name = "Site Groups",
      labels = c("F18 Significant Traits", "F18 Signif.-Traits Neutral", "F18 Trait Associated", "F18 Trait-Associated Neutral", "F58 Significant Traits", "F58 Signif.-Traits Neutral", "F58 Trait Associated", "F58 Trait-Associated Neutral"))
dev.off()


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



### bar plots w all groups for each generation
t_18 <- df4[grep("18", df4$group), ]
#t_18$group <- factor(t_18$group, levels = c("sig18", "signeu18", "trait18", "trneu18"))

t_58 <- df4[grep("58", df4$group), ]
#t_58$group <- factor(t_58$group, levels = c("sig58", "signeu58", "trait58", "trneu58"))

g1 <- t_18 %>%
    ggplot(aes(x=bins, y=value, fill=group)) +
      geom_bar(stat='identity', position='dodge', width=0.7) +
      theme_classic(base_size=20) +
      theme(axis.text.x = element_text(angle=45, hjust=1)) +
      labs(x="Frequency Ranges", y="Proportion of Sites", title="F18 Site Frequency Spectrum") +
      scale_fill_discrete(name = "Site Groups", labels = c("Significant Trait Associated", "Significant-Trait Neutral", "Trait Associated", "Trait-Associated Neutral"))

g2 <- t_58 %>%
    ggplot(aes(x=bins, y=value, fill=group)) +
      geom_bar(stat='identity', position='dodge', width=0.7) +
      theme_classic(base_size=20) +
      theme(axis.text.x = element_text(angle=45, hjust=1)) +
      labs(x="Frequency Ranges", y="Proportion of Sites", title="F58 Site Frequency Spectrum") +
      scale_fill_discrete(name = "Site Groups", labels = c("Significant Trait Associated", "Significant-Trait Neutral", "Trait Associated", "Trait-Associated Neutral"))


png("SFS_2GenerationComparisons.png", width=30, height=8, units="in", res=300)
ggarrange(g1, g2)
dev.off()


#+ #scale_fill_manual(values = col10)
