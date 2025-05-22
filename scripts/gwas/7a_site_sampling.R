#!/usr/bin/env Rscript

#SBATCH --job-name=gwas
#SBATCH -o /rhome/jmarz001/bigdata/Ag-Competition/scripts/gwas/7a_site_sampling.stdout
#SBATCH --ntasks=1
#SBATCH --mem=80gb
#SBATCH -t 02:00:00
#SBATCH -p koeniglab

library(pacman)
p_load(tidyverse, data.table, gridExtra, ggsci, ggpubr, Cairo)

setwd("/rhome/jmarz001/bigdata/Ag-Competition/results/gwas/POP_AF")
#source("../../../scripts/standardize_trait_text_function.R")
df <- fread("MAF_BeneficialAllele_TraitAssoc.tsv")

# list of significantly changed traits
sig <- fread("../../002A_change_over_time_derived_traits_significant.tsv")
st <- unlist(sig$traits)

# divide sites into 3 groups
# 1 sites associated w significantly changed trait
s1 <- df %>% filter(trait_name %in% st) %>% select(c('chr', 'rs', 'ps', 'F18_AF', 'F58_AF', 'trait_name'))
#%>% select(-c('F0A1', 'F0A2', 'F18A1', 'F18A2', 'F58A1', 'F58A2', 'chi_p', 'minor_allele', 'F0_AF', 'DELTA_P18', 'DELTA_1858'))
dim(s1)
# number of unique sites
unique(s1$rs) %>% length
# list of unique sites' frequencies for binnings
s1u <- unique(s1[,1:5])
write_delim(s1u, "SigTraits_sites.tsv")

# 2 sites associated w any other traits
s2 <- df %>% filter(!is.na(trait_name)) %>% filter(!(trait_name %in% st))  %>% select(c('chr', 'rs', 'ps', 'F18_AF', 'F58_AF', 'trait_name'))
dim(s2)
# number of unique sites
unique(s2$rs) %>% length
# list of unique sites' frequencies for binnings
s2u <- unique(s2[,1:5])
write_delim(s2u, "NeutralTraits_sites.tsv")

# 3 all other sites
s3 <- df %>% filter(is.na(trait_name))  %>% select(c('chr', 'rs', 'ps', 'F18_AF', 'F58_AF', 'trait_name'))
dim(s3)
# number of unique sites
unique(s3$rs) %>% length
s3u <- unique(s3[,1:5])



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
    #hist <- hist/dim(x)[1]
    return(hist)
}


# apply plotting function and make list of those plots
b14 <- binning(s1[,4])
b15 <- binning(s1[,5])
b24 <- binning(s2[,4])
b25 <- binning(s2[,5])

# number of samples to pull to match starting site freqs
b14u <- binning(s1u[,4])
b24u <- binning(s2u[,4])


b34 <- binning(s3[,4])
b35 <- binning(s3[,5])

# make frequency bin df
pop_freqs <- tibble('bins'=bins, 'F18_sig'=b14, 'F58_sig'=b15, 'F18_trait'=b24, 'F58_trait'=b25, 'F18_neu'=b34, 'F58_neu'=b35)
# determines sample number to each AF bin


# determine size of samples to match
s1num <- sum(b14u)
s2num <- sum(b24u)
#s1 ~80
#s2 ~100




# bin s3 sites for easy sampling?
# reuse binning function, but instead of counting sies for each frequencies bin,
# store the chr-position
# sites in neutral site list that belong to each bin
s3s <- s3 %>% select(c(rs, F18_AF))
s3 <- s3 %>% select(c(rs, F18_AF, F58_AF))
s3$delta <- s3$F58_AF - s3$F18_AF


sort_site_by_freq <- function(x) {
    bins = c(0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0)

    a = c()
    b = c()
    c = c()
    d = c()
    e = c()
    f = c()
    g = c()
    h = c()
    i = c()
    j = c()
    k = c()
    l = c()
    m = c()
    n = c()
    o = c()
    p = c()
    q = c()
    r = c()
    s = c()
    t = c()
    u = c()

    for(row in 1:nrow(x)) {
        # take in the frequency
        freq = unlist(x[row, 2])
        # when the bin matches, store the site in that bin
        site = unlist(x[row, 1])

        if(freq >= bins[1] & freq < bins[2]) {
          a[length(a)+1] = site
        } else if(freq >= bins[2] & freq < bins[3]) {
          b[length(b)+1] = site
        } else if(freq >= bins[3] & freq < bins[4]) {
          c[length(c)+1] = site
        } else if(freq >= bins[4] & freq < bins[5]) {
          d[length(d)+1] = site
        } else if(freq >= bins[5] & freq < bins[6]) {
          e[length(e)+1] = site
        } else if(freq >= bins[6] & freq < bins[7]) {
          f[length(f)+1] = site
        } else if(freq >= bins[7] & freq < bins[8]) {
          g[length(g)+1] = site
        } else if(freq >= bins[8] & freq < bins[9]) {
          h[length(h)+1] = site
        } else if(freq >= bins[9] & freq < bins[10]) {
          i[length(i)+1] = site
        } else if(freq >= bins[10] & freq < bins[11]) {
          j[length(j)+1] = site
        } else if(freq >= bins[11] & freq < bins[12]) {
          k[length(k)+1] = site
        } else if(freq >= bins[12] & freq < bins[13]) {
          l[length(l)+1] = site
        } else if(freq >= bins[13] & freq < bins[14]) {
          m[length(m)+1] = site
        } else if(freq >= bins[14] & freq < bins[15]) {
          n[length(n)+1] = site
        } else if(freq >= bins[15] & freq < bins[16]) {
          o[length(o)+1] = site
        } else if(freq >= bins[16] & freq < bins[17]) {
          p[length(p)+1] = site
        } else if(freq >= bins[17] & freq < bins[18]) {
          q[length(q)+1] = site
        } else if(freq >= bins[18] & freq < bins[19]) {
          r[length(r)+1] = site
        } else if(freq >= bins[19] & freq < bins[20]) {
          s[length(s)+1] = site
        } else if(freq >= bins[20] & freq <= bins[21]) {
          t[length(t)+1] = site
        } else {
          u[length(u)+1] = site
        }
  }
    freq_sites <- list(a, b, c, d, e, f, g, h, i, j, k, l, m, n, o, p, q, r, s, t)
    if(length(u) > 0){
      print("ERROR: uncategorized site frequencies")
    }
    return(freq_sites)
}

sfq <- sort_site_by_freq(s3s)
# are any of the freq categories empty?
for(i in 1:length(sfq)){print(length(sfq[[i]]))}



neutral_site_AF <- function(sample_sites=sfq, sample_freqs= s3, freq_sample_numbers=b14u) {
        change_sample <- c()
        early_sample <- c()
        late_sample <- c()

        # sample loop 1000 times
        for(i in 1:1000){
          # run through each frequency bin in turn
          F18_AF_list <- c()
          F58_AF_list <- c()
          new_delta_list <- c()

          # for each frequency bin
            ## each frequency bin is a list of site IDs
          for(j in 1:length(sample_sites)){
            bbin <- sample_freqs[[1]]
            sample_num <- freq_sample_numbers[j]
            # sample the sites from each bin
            fx <- sample(bbin, sample_num, replace=F)

            # pull allele frequencies at the sites
            daf <- sample_freqs %>% filter(rs %in% fx) %>% select(delta) %>% unlist
            af1 <- sample_freqs %>% filter(rs %in% fx) %>% select(F18_AF) %>% unlist
            af2 <- sample_freqs %>% filter(rs %in% fx) %>% select(F58_AF) %>% unlist

            # and save it for averaging over all freq bins
            new_delta_list <- c(new_delta_list, daf)
            F18_AF_list <- c(F18_AF_list, af1)
            F58_AF_list <- c(F58_AF_list, af2)
          }

          # average change in AF over sites
          # add sample average to your neutral site AF list
          change_sample <- c(change_sample, mean(new_delta_list))
          early_sample <- c(early_sample, mean(F18_AF_list))
          late_sample <- c(late_sample, mean(F58_AF_list))
        }

        xyz <- tibble("F18_AF_neutral"=early_sample, "F58_AF_neutral"=late_sample, "delta_AF_neutral"=change_sample)
        return(xyz)
      }



safe <- neutral_site_AF(sample_sites=sfq, sample_freqs= s3, freq_sample_numbers=b14u)
saff <- neutral_site_AF(sample_sites=sfq, sample_freqs= s3, freq_sample_numbers=b24u)


write_delim(safe, "SigTraits_neutral_sites.tsv", "\t")
write_delim(saff, "AssocTraits_neutral_sites.tsv", "\t")
# the average change in allele frequency
# for sites with the same starting allele frequency
# randomly sample 1000 times to simulate the random probability paths of frequency changes
