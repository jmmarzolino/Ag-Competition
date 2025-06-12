#!/usr/bin/env Rscript

#SBATCH --job-name=gwas
#SBATCH -o /rhome/jmarz001/bigdata/Ag-Competition/scripts/gwas/7a_site_sampling.stdout
#SBATCH --ntasks=1
#SBATCH --mem=80gb
#SBATCH -t 02:00:00
#SBATCH -p koeniglab

library(pacman)
p_load(tidyverse, data.table, gridExtra, ggsci, ggpubr, Cairo)

setwd("/rhome/jmarz001/bigdata/Ag-Competition/results/gwas")
source("../../scripts/CUSTOM_FNS.R")


## all sites, all generations, allele frequencies
## read in all sites' allele frequencies
df <- fread("all_sites_allele_freqs.tsv")
## should already be filtered to the right columns
df$snp <- paste0(df$CHR, "_", df$BP)


## polarize all sites to the BENEFICIAL allele (F18 to F58 rather than starting w F0...)
# which is the one which increases freq over gens
# do that for all sites
df$beneficial_allele <- "A1"

x <- round(df$F58_AF - df$F18_AF, 4)
dfa1 <- df[which(x >= 0), ]
dfa2 <- df[which(x < 0), ]
# sites w beneficial allele 2...
dfa2$beneficial_allele <- "A2"
dfa2[, 5:9] <- 1 - (dfa2[, 5:9])

# stitch data sets back together
df2 <- bind_rows(dfa1, dfa2)

## then pull sites associated w traits
## and filter them to their own dataset
# load associated sites list
all_sig <- fread("all_gwas_sig_sites.tsv")
all_sig$rs <- paste0(all_sig$chr, "_", all_sig$ps)

trait_assoc_sites <- df2 %>% filter(snp %in% all_sig$rs)
neutral_sites <- df2 %>% filter(!(snp %in% all_sig$rs))


## find starting allele freqs for all trait associated sites
#### Site Frequency Binning
# define other necessary vars, the breaks & ranges
bins = c("0.00-0.05", "0.05-0.1", '0.1-0.15', '0.15-0.2', '0.2-0.25', '0.25-0.3', '0.3-0.35', '0.35-0.4', '0.4-0.45', '0.45-0.5', '0.5-0.55', '0.55-0.6', '0.6-0.65', '0.65-0.7', '0.7-0.75', '0.75-0.8', '0.8-0.85', '0.85-0.9', '0.9-0.95', '0.95-1') #, 'unknown')
#levels(cut(0:1, 20))
#levels(cut(0:1, 10))
bigbins = c("0.0-0.1", '0.1-0.2', '0.2-0.3', '0.3-0.4', '0.4-0.5', '0.5-0.6', '0.6-0.7', '0.7-0.8', '0.8-0.9', '0.9-1', 'unknown')


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

# divide site list into allele freq bins
## one list starting w F0 and one w F18

# apply plotting function and make list of those plots
# allele freq bins for gwas sites, F0 and F18
bin_F0 <- binning(trait_assoc_sites[,5])
bin_F18 <- binning(trait_assoc_sites[,6])

# number of samples to pull to match starting site freqs
sample_n <- nrow(trait_assoc_sites)



# make frequency bin df
pop_freqs <- tibble('bins'=bins, 'F0_bins'=bin_F0, 'F18_bins'=bin_F18)
# determines sample number to each AF bin
fwrite(pop_freqs, "gwas_sites_pop_freq_binned.tsv")



## with remaining (neutral) sites of genome, 
## divide those sites into bins based on their starting allele freq
neutral_sites
bin_F0_neutral <- binning(neutral_sites[,5])
bin_F18_neutral <- binning(neutral_sites[,6])
pop_freqs_neutral <- tibble('bins'=bins, 'F0_bins_neutral'=bin_F0_neutral, 'F18_bins_neutral'=bin_F18_neutral)

fwrite(pop_freqs_neutral, "neutral_sites_pop_freq_binned.tsv")

full_join(pop_freqs, pop_freqs_neutral)

### bin all neutral sites based on their starting frequency
### this will make sampling based on starting frequency easy 
# reuse binning function, but instead of counting sites for each frequencies bin,
# store the chr-position sites in neutral site list that belong to each bin

sort_site_by_freq <- function(x=netral_sites, gen_col="F18_AF") {
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
        freq = unlist(x[row, ..gen_col])
        # when the bin matches, store the site in that bin
        site = unlist(x[row, snp])

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

sfq18 <- sort_site_by_freq(neutral_sites, "F18_AF")
sfq0 <- sort_site_by_freq(neutral_sites, "F0_AF")
# are any of the freq categories empty?
for(i in 1:length(sfq0)){
  print(length(sfq0[[i]])) }
for(i in 1:length(sfq18)){
  print(length(sfq18[[i]])) }


# randomly sample a matching number of starting allele freq sites
# record the allele freq change over gens for each of those sites
# record median & absolute val median
# repeat 1000 times


#sample sites = table of genome sites binned by frequency
# list of genome-site lists binned based on allele frequency
# or, a list of neutral sites for a given starting allele freq 
#bins
#str(sfq0)
#str(sfq18)
# sample_freqs = ??? 's3' full data?
#neutral_sites$delta_F0F18 <- neutral_sites$F18_AF - neutral_sites$F0_AF
#neutral_sites$delta_F18F58 <- neutral_sites$F58_AF - neutral_sites$F18_AF
# freq_sample_numbers = 
#pop_freqs
#sfq0 & F0_bins 
#sfq18 & F18_bins 

neutral_site_AF <- function(sample_sites=sfq0, sample_freqs= neutral_sites, freq_sample_numbers=pop_freqs$F0_bins, EARLY=0, LATE=58) {
        delta_sample <- c()
        delta_abs_sample <- c()
        early_sample <- c()
        late_sample <- c()

        
        # define 'early' and 'late' generation cols
        earlygen <- paste0("F", EARLY, "_AF")
        lategen <- paste0("F", LATE, "_AF")


        # sample loop 1000 times
        for(i in 1:1000){
          # run through each frequency bin in turn
            early_AF_list <- c()
            late_AF_list <- c()
            delta_list <- c()

          # for each frequency bin
            ## each frequency bin is a list of site IDs
          for(j in 1:length(sample_sites)){
            sites <- sample_sites[[j]] # select the list of IDs
            sample_num <- freq_sample_numbers[j] # select the number of IDs to sample from the list
            fx <- sample(sites, sample_num, replace=F) # sample the sites from each bin

            # pull allele frequencies at the sites
            row <- sample_freqs %>% filter(snp %in% fx)

            earlygenaf <- row %>% select(all_of(earlygen)) %>% unlist
            lategenaf <- row %>% select(all_of(lategen)) %>% unlist
            dddelta <- lategenaf - earlygenaf

            # and save it for averaging over all freq bins
            early_AF_list <- c(early_AF_list, earlygenaf)
            late_AF_list <- c(late_AF_list, lategenaf)
            delta_list <- c(delta_list, dddelta)
          }

          # find the median change in AF over sites
          # add sampled median to your neutral-site AF change list
          delta_sample <- c(delta_sample, median(delta_list))
          delta_abs_sample <- c(delta_abs_sample, median(abs(delta_list)))
          early_sample <- c(early_sample, median(early_AF_list))
          late_sample <- c(late_sample, median(late_AF_list))
        }

        xyz <- tibble("early_AF_neutral"=early_sample, "late_AF_neutral"=late_sample, "delta_AF_neutral"=delta_sample, "delta_AF_neutral_absolute"=delta_abs_sample)
        return(xyz)
      }



random_sample_neutral_F0 <- neutral_site_AF(sample_sites=sfq0, sample_freqs= neutral_sites, freq_sample_numbers=pop_freqs$F0_bins, EARLY=0, LATE=58)
random_sample_neutral_F18 <- neutral_site_AF(sample_sites=sfq18, sample_freqs= neutral_sites, freq_sample_numbers=pop_freqs$F18_bins, EARLY=18, LATE=58)


write_delim(random_sample_neutral_F0, "neutral_sites_sampled_F0.tsv", "\t")
write_delim(random_sample_neutral_F18, "neutral_sites_sampled_F18.tsv", "\t")
# the average change in allele frequency
# for sites with the same starting allele frequency
# randomly sample 1000 times to simulate the random probability paths of frequency changes



# plot allele freq & allele freq chnages as bar plot/allele freq spectrum
# next script