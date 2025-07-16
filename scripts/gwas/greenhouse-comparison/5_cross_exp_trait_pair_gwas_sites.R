#!/usr/bin/env Rscript
#SBATCH -o /rhome/jmarz001/bigdata/Ag-Competition/scripts/gwas/greenhouse-comparison/5_cross_exp_trait_pair_gwas_sites.stdout
#SBATCH --ntasks=1
#SBATCH --mem=64gb
#SBATCH -t 01:30:00
#SBATCH -p koeniglab

library(pacman)
p_load(tidyverse, data.table, gridExtra, ggsci, Cairo, qqman, CMplot)
options(stringsAsFactors = F)

setwd("/rhome/jmarz001/bigdata/Ag-Competition/results/gwas/CCII_greenhouse_exp_gwas")
source("../../../scripts/CUSTOM_FNS.R")

## it's only a few plots so you don't need to do it so systematically but can just bang it out in repetetive code blocks per trait pair
# yeah for lm and lmm

# trait pairs
## FT - Flowering_days_2017
## FT - Flowering_2018_Median
## TOTAL_MASS - Total_mass
## SEED_WEIGHT_100 - Mass_100
## FECUNDITY - Seed_Estimate

# manually compile matched lists of traits & files
matched_exp_ft <- tibble(
          "experiment"=c("field", "field", "greenhouse", "greenhouse", "greenhouse", "greenhouse") ,
          "trait"=c("FT", "FT", "Flowering_days_2017", "Flowering_days_2017", "Flowering_2018_Median", "Flowering_2018_Median") , 
          "gwas_method"=c("lm", "lmm", "lm", "lmm", "lm", "lmm") ,
          "file"=c("../ASSOC_6.assoc.txt", "../ASSOC_6_lmm.assoc.txt", "ASSOC_12.assoc.txt", "ASSOC_12_lmm.assoc.txt", "ASSOC_13.assoc.txt", "ASSOC_13_lmm.assoc.txt")
          )

matched_exp_fec <- tibble(
          "experiment"=c("field", "field", "greenhouse", "greenhouse") ,
          "trait"=c("FECUNDITY", "FECUNDITY", "Seed_Estimate", "Seed_Estimate") , 
          "gwas_method"=c("lm", "lmm", "lm", "lmm") ,
          "file"=c("../ASSOC_11.assoc.txt", "../ASSOC_11_lmm.assoc.txt", "ASSOC_11.assoc.txt", "ASSOC_11_lmm.assoc.txt")
          )

matched_exp_sw <- tibble(
          "experiment"=c("field", "field", "greenhouse", "greenhouse") ,
          "trait"=c("SEED_WEIGHT_100", "SEED_WEIGHT_100", "Mass_100", "Mass_100") , 
          "gwas_method"=c("lm", "lmm", "lm", "lmm") ,
          "file"=c("../ASSOC_8.assoc.txt", "../ASSOC_8_lmm.assoc.txt", "ASSOC_9.assoc.txt", "ASSOC_9_lmm.assoc.txt")
          )


matched_exp_tm <- tibble(
          "experiment"=c("field", "field", "greenhouse", "greenhouse") ,
          "trait"=c("TOTAL_MASS", "TOTAL_MASS", "Total_mass", "Total_mass") , 
          "gwas_method"=c("lm", "lmm", "lm", "lmm") ,
          "file"=c("../ASSOC_7.assoc.txt", "../ASSOC_7_lmm.assoc.txt", "ASSOC_10.assoc.txt", "ASSOC_10_lmm.assoc.txt")
          )

#combo_traits <- bind_rows(matched_exp_ft, matched_exp_fec, matched_exp_sw, matched_exp_tm)



# read in pairs of tables
# join them together by position
# can filter high p-val rows & have a table of common gwas sites across diff tests & write out


### flowering time overlap
x1 <- fread(matched_exp_ft[[1,4]]) %>% select(c("chr", "ps", "af", "beta", "p_lrt"))
x2 <- fread(matched_exp_ft[[2,4]]) %>% select(c("chr", "ps", "af", "beta", "p_lrt"))
x3 <- fread(matched_exp_ft[[3,4]]) %>% select(c("chr", "ps", "af", "beta", "p_lrt"))
x4 <- fread(matched_exp_ft[[4,4]]) %>% select(c("chr", "ps", "af", "beta", "p_lrt"))
x5 <- fread(matched_exp_ft[[5,4]]) %>% select(c("chr", "ps", "af", "beta", "p_lrt"))
x6 <- fread(matched_exp_ft[[6,4]]) %>% select(c("chr", "ps", "af", "beta", "p_lrt"))

## join datasets, keep the higher p-val for lm and lmm 
# field flowering time
y1 <- full_join(x1, x2, by=c("chr", "ps", "af"), suffix=c("_lm", "_lmm"))
y1$p <- y1$p_lrt_lm
y1[which(y1$p_lrt_lm < y1$p_lrt_lmm), 8] <- y1[which(y1$p_lrt_lm < y1$p_lrt_lmm), p_lrt_lmm]
y1 <- y1 %>% select(c(chr, ps, af, p))
y1$group <- "field"
# greenhouse flowering time 2017
y2 <- full_join(x3, x4, by=c("chr", "ps", "af"), suffix=c("_lm", "_lmm"))
y2$p <- y2$p_lrt_lm
y2[which(y2$p_lrt_lm < y2$p_lrt_lmm), 8] <- y2[which(y2$p_lrt_lm < y2$p_lrt_lmm), p_lrt_lmm]
y2 <- y2 %>% select(c(chr, ps, af, p))
y2$group <- "greenhouse_2017"
# greenhouse flowering time 2018
y3 <- full_join(x5, x6, by=c("chr", "ps", "af"), suffix=c("_lm", "_lmm"))
y3$p <- y3$p_lrt_lm
y3[which(y3$p_lrt_lm < y3$p_lrt_lmm), 8] <- y3[which(y3$p_lrt_lm < y3$p_lrt_lmm), p_lrt_lmm]
y3 <- y3 %>% select(c(chr, ps, af, p))
y3$group <- "greenhouse_2018"

# now join the three experiments
z <- full_join(y1, y2)
z <- full_join(z, y3)

# filter to top 5%-ish??? of sites across experiments?
threshhold <- quantile(z$p, 0.0005)
z_filt <- z[which(z$p <= threshhold), ]
OutName <- paste0("cross_exp_gwas_suggestive_peaks_", "flowering_time", ".tsv")
fwrite(z_filt, OutName)

# chr col needs to be numeric
z$chr <- as.numeric(gsub("chr(\\d)H", "\\1", z$chr))
z$snp <- paste0(z$chr, "_", z$ps)
tmp <- z %>% select(c(snp, chr, ps, p, group))
tmp2 <- tmp %>% pivot_wider(names_from="group", values_from=p)

# plot type 1 - CMplot
CMplot(tmp2, type="p", plot.type="m", LOG10=TRUE, 
        threshold=-log10(bon), threshold.col="black", threshold.lty=1, col=c("grey60","#4197d8"), 
        signal.cex=1.2, signal.col="red", verbose=TRUE, multracks=TRUE, file="png", dpi=300, file.output=TRUE)
        # Plots are stored in: /bigdata/koeniglab/jmarz001/Ag-Competition/results/gwas/CCII_greenhouse_exp_gwas

# plot type 2 - all-in-one vertical manhattan
png("combo_ft_manhattan.png")
manhattan(tmp, chr="chr", bp="ps", p="p", snp="snp", genomewideline=-log10(bon), suggestiveline=-log10(bon)-1)
dev.off()


# plot type 3 - all-in-one overlaping manhattan
tmp$BP <- as.numeric(tmp$ps)
threshold <- 0.05/length(y1$p)
result <- tmp %>%
  # Compute chromosome size
  group_by(chr) %>%
  summarise(chr_len = max(BP)) %>%
  # Calculate cumulative position of each chromosome
  mutate(tot = cumsum(as.numeric(chr_len))-as.numeric(chr_len)) %>%
  select(-chr_len) %>%
  # Add this info to the initial dataset
  left_join(tmp, ., by=c("chr" = "chr")) %>%
  # Add a cumulative position of each SNP
  arrange(chr, BP) %>%
  mutate(BPcum = BP+tot)

#result <- result %>% filter(-log10(p_lrt)>2)
axisdf <- result %>% group_by(chr) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

png("cross_exp_manhattan_overlap_ft.png")
# Manhattan plot
ggplot(result, aes(x=BPcum, y=-log10(p), group=group)) +
    # Show all points
    geom_point(aes(color=group, shape=as.factor(chr)), alpha=0.5, size=1.3) +
    geom_hline(aes(yintercept=-log10(threshold)), color = "firebrick1", linetype="dashed", alpha=0.7) +
    scale_shape_manual(values=rep(c(16, 18), 22)) + 
    scale_color_manual(values = c("dark green", "dodgerblue4", "deepskyblue")) +
    # custom X axis:
    scale_x_continuous(label = axisdf$chr, breaks= axisdf$center) +
    scale_y_continuous(expand = c(0, 0.5)) +   # remove space between plot area and x axis
    # Custom theme:
    theme_classic() +
    theme(#legend.position="none",
          panel.border = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          text=element_text(size=16),
          plot.title = element_text(size=12),
          plot.subtitle = element_text(size=10)) +
    xlab("Chromosome") +
    ylab(expression(-log[10](italic(p)))) +
    ggtitle("Flowering Time GWAS Results across Experiments")
dev.off()






### fecundity
x1 <- fread(matched_exp_fec[[1,4]]) %>% select(c("chr", "ps", "af", "beta", "p_lrt"))
x2 <- fread(matched_exp_fec[[2,4]]) %>% select(c("chr", "ps", "af", "beta", "p_lrt"))
x3 <- fread(matched_exp_fec[[3,4]]) %>% select(c("chr", "ps", "af", "beta", "p_lrt"))
x4 <- fread(matched_exp_fec[[4,4]]) %>% select(c("chr", "ps", "af", "beta", "p_lrt"))

## join datasets, keep the higher p-val for lm and lmm 
# field flowering time
y1 <- full_join(x1, x2, by=c("chr", "ps", "af"), suffix=c("_lm", "_lmm"))
y1$p <- y1$p_lrt_lm
y1[which(y1$p_lrt_lm < y1$p_lrt_lmm), 8] <- y1[which(y1$p_lrt_lm < y1$p_lrt_lmm), p_lrt_lmm]
y1 <- y1 %>% select(c(chr, ps, af, p))
y1$group <- "field"
# greenhouse flowering time 2017
y2 <- full_join(x3, x4, by=c("chr", "ps", "af"), suffix=c("_lm", "_lmm"))
y2$p <- y2$p_lrt_lm
y2[which(y2$p_lrt_lm < y2$p_lrt_lmm), 8] <- y2[which(y2$p_lrt_lm < y2$p_lrt_lmm), p_lrt_lmm]
y2 <- y2 %>% select(c(chr, ps, af, p))
y2$group <- "greenhouse"

# now join the three experiments
z <- full_join(y1, y2)

# filter to top 5%-ish??? of sites across experiments?
threshhold <- quantile(z$p, 0.0005)
z_filt <- z[which(z$p <= threshhold), ]
OutName <- paste0("cross_exp_gwas_suggestive_peaks_", "fecundity", ".tsv")
fwrite(z_filt, OutName)

# chr col needs to be numeric
z$chr <- as.numeric(gsub("chr(\\d)H", "\\1", z$chr))
z$snp <- paste0(z$chr, "_", z$ps)
tmp <- z %>% select(c(snp, chr, ps, p, group))

# plot type 3 - all-in-one overlaping manhattan
tmp$BP <- as.numeric(tmp$ps)
threshold <- 0.05/length(y1$p)
result <- tmp %>%
  # Compute chromosome size
  group_by(chr) %>%
  summarise(chr_len = max(BP)) %>%
  # Calculate cumulative position of each chromosome
  mutate(tot = cumsum(as.numeric(chr_len))-as.numeric(chr_len)) %>%
  select(-chr_len) %>%
  # Add this info to the initial dataset
  left_join(tmp, ., by=c("chr" = "chr")) %>%
  # Add a cumulative position of each SNP
  arrange(chr, BP) %>%
  mutate(BPcum = BP+tot)

#result <- result %>% filter(-log10(p_lrt)>2)
axisdf <- result %>% group_by(chr) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

png("cross_exp_manhattan_overlap_fec.png")
# Manhattan plot
ggplot(result, aes(x=BPcum, y=-log10(p), group=group)) +
    # Show all points
    geom_point(aes(color=group, shape=as.factor(chr)), alpha=0.5, size=1.3) +
    geom_hline(aes(yintercept=-log10(threshold)), color = "firebrick1", linetype="dashed", alpha=0.7) +
    scale_shape_manual(values=rep(c(16, 18), 22)) + 
    scale_color_manual(values = c("dodgerblue4", "deepskyblue")) +
    # custom X axis:
    scale_x_continuous(label = axisdf$chr, breaks= axisdf$center) +
    scale_y_continuous(expand = c(0, 0.5)) +   # remove space between plot area and x axis
    # Custom theme:
    theme_classic() +
    theme(#legend.position="none",
          panel.border = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          text=element_text(size=16),
          plot.title = element_text(size=12),
          plot.subtitle = element_text(size=10)) +
    xlab("Chromosome") +
    ylab(expression(-log[10](italic(p)))) +
    ggtitle("Fecundity GWAS Results across Experiments")
dev.off()



### total mass
x1 <- fread(matched_exp_tm[[1,4]]) %>% select(c("chr", "ps", "af", "beta", "p_lrt"))
x2 <- fread(matched_exp_tm[[2,4]]) %>% select(c("chr", "ps", "af", "beta", "p_lrt"))
x3 <- fread(matched_exp_tm[[3,4]]) %>% select(c("chr", "ps", "af", "beta", "p_lrt"))
x4 <- fread(matched_exp_tm[[4,4]]) %>% select(c("chr", "ps", "af", "beta", "p_lrt"))

## join datasets, keep the higher p-val for lm and lmm 
# field flowering time
y1 <- full_join(x1, x2, by=c("chr", "ps", "af"), suffix=c("_lm", "_lmm"))
y1$p <- y1$p_lrt_lm
y1[which(y1$p_lrt_lm < y1$p_lrt_lmm), 8] <- y1[which(y1$p_lrt_lm < y1$p_lrt_lmm), p_lrt_lmm]
y1 <- y1 %>% select(c(chr, ps, af, p))
y1$group <- "field"
# greenhouse flowering time 2017
y2 <- full_join(x3, x4, by=c("chr", "ps", "af"), suffix=c("_lm", "_lmm"))
y2$p <- y2$p_lrt_lm
y2[which(y2$p_lrt_lm < y2$p_lrt_lmm), 8] <- y2[which(y2$p_lrt_lm < y2$p_lrt_lmm), p_lrt_lmm]
y2 <- y2 %>% select(c(chr, ps, af, p))
y2$group <- "greenhouse"

# now join the three experiments
z <- full_join(y1, y2)

# filter to top 5%-ish??? of sites across experiments?
threshhold <- quantile(z$p, 0.0005)
z_filt <- z[which(z$p <= threshhold), ]
OutName <- paste0("cross_exp_gwas_suggestive_peaks_", "total_mass", ".tsv")
fwrite(z_filt, OutName)

# chr col needs to be numeric
z$chr <- as.numeric(gsub("chr(\\d)H", "\\1", z$chr))
z$snp <- paste0(z$chr, "_", z$ps)
tmp <- z %>% select(c(snp, chr, ps, p, group))

# plot type 3 - all-in-one overlaping manhattan
tmp$BP <- as.numeric(tmp$ps)
threshold <- 0.05/length(y1$p)
result <- tmp %>%
  # Compute chromosome size
  group_by(chr) %>%
  summarise(chr_len = max(BP)) %>%
  # Calculate cumulative position of each chromosome
  mutate(tot = cumsum(as.numeric(chr_len))-as.numeric(chr_len)) %>%
  select(-chr_len) %>%
  # Add this info to the initial dataset
  left_join(tmp, ., by=c("chr" = "chr")) %>%
  # Add a cumulative position of each SNP
  arrange(chr, BP) %>%
  mutate(BPcum = BP+tot)

#result <- result %>% filter(-log10(p_lrt)>2)
axisdf <- result %>% group_by(chr) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

png("cross_exp_manhattan_overlap_total_mass.png")
# Manhattan plot
ggplot(result, aes(x=BPcum, y=-log10(p), group=group)) +
    # Show all points
    geom_point(aes(color=group, shape=as.factor(chr)), alpha=0.5, size=1.3) +
    geom_hline(aes(yintercept=-log10(threshold)), color = "firebrick1", linetype="dashed", alpha=0.7) +
    scale_shape_manual(values=rep(c(16, 18), 22)) + 
    scale_color_manual(values = c("dodgerblue4", "deepskyblue")) +
    # custom X axis:
    scale_x_continuous(label = axisdf$chr, breaks= axisdf$center) +
    scale_y_continuous(expand = c(0, 0.5)) +   # remove space between plot area and x axis
    # Custom theme:
    theme_classic() +
    theme(#legend.position="none",
          panel.border = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          text=element_text(size=16),
          plot.title = element_text(size=12),
          plot.subtitle = element_text(size=10)) +
    xlab("Chromosome") +
    ylab(expression(-log[10](italic(p)))) +
    ggtitle("Total Mass GWAS Results across Experiments")
dev.off()








### 100 seed weight
x1 <- fread(matched_exp_sw[[1,4]]) %>% select(c("chr", "ps", "af", "beta", "p_lrt"))
x2 <- fread(matched_exp_sw[[2,4]]) %>% select(c("chr", "ps", "af", "beta", "p_lrt"))
x3 <- fread(matched_exp_sw[[3,4]]) %>% select(c("chr", "ps", "af", "beta", "p_lrt"))
x4 <- fread(matched_exp_sw[[4,4]]) %>% select(c("chr", "ps", "af", "beta", "p_lrt"))

## join datasets, keep the higher p-val for lm and lmm 
# field flowering time
y1 <- full_join(x1, x2, by=c("chr", "ps", "af"), suffix=c("_lm", "_lmm"))
y1$p <- y1$p_lrt_lm
y1[which(y1$p_lrt_lm < y1$p_lrt_lmm), 8] <- y1[which(y1$p_lrt_lm < y1$p_lrt_lmm), p_lrt_lmm]
y1 <- y1 %>% select(c(chr, ps, af, p))
y1$group <- "field"
# greenhouse flowering time 2017
y2 <- full_join(x3, x4, by=c("chr", "ps", "af"), suffix=c("_lm", "_lmm"))
y2$p <- y2$p_lrt_lm
y2[which(y2$p_lrt_lm < y2$p_lrt_lmm), 8] <- y2[which(y2$p_lrt_lm < y2$p_lrt_lmm), p_lrt_lmm]
y2 <- y2 %>% select(c(chr, ps, af, p))
y2$group <- "greenhouse"

# now join the three experiments
z <- full_join(y1, y2)

# filter to top 5%-ish??? of sites across experiments?
threshhold <- quantile(z$p, 0.0005)
z_filt <- z[which(z$p <= threshhold), ]
OutName <- paste0("cross_exp_gwas_suggestive_peaks_", "100_seed_weight", ".tsv")
fwrite(z_filt, OutName)

# chr col needs to be numeric
z$chr <- as.numeric(gsub("chr(\\d)H", "\\1", z$chr))
z$snp <- paste0(z$chr, "_", z$ps)
tmp <- z %>% select(c(snp, chr, ps, p, group))

# plot type 3 - all-in-one overlaping manhattan
tmp$BP <- as.numeric(tmp$ps)
threshold <- 0.05/length(y1$p)
result <- tmp %>%
  # Compute chromosome size
  group_by(chr) %>%
  summarise(chr_len = max(BP)) %>%
  # Calculate cumulative position of each chromosome
  mutate(tot = cumsum(as.numeric(chr_len))-as.numeric(chr_len)) %>%
  select(-chr_len) %>%
  # Add this info to the initial dataset
  left_join(tmp, ., by=c("chr" = "chr")) %>%
  # Add a cumulative position of each SNP
  arrange(chr, BP) %>%
  mutate(BPcum = BP+tot)

#result <- result %>% filter(-log10(p_lrt)>2)
axisdf <- result %>% group_by(chr) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

png("cross_exp_manhattan_overlap_100sw.png")
# Manhattan plot
ggplot(result, aes(x=BPcum, y=-log10(p), group=group)) +
    # Show all points
    geom_point(aes(color=group, shape=as.factor(chr)), alpha=0.5, size=1.3) +
    geom_hline(aes(yintercept=-log10(threshold)), color = "firebrick1", linetype="dashed", alpha=0.7) +
    scale_shape_manual(values=rep(c(16, 18), 22)) + 
    scale_color_manual(values = c("dodgerblue4", "deepskyblue")) +
    # custom X axis:
    scale_x_continuous(label = axisdf$chr, breaks= axisdf$center) +
    scale_y_continuous(expand = c(0, 0.5)) +   # remove space between plot area and x axis
    # Custom theme:
    theme_classic() +
    theme(#legend.position="none",
          panel.border = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          text=element_text(size=16),
          plot.title = element_text(size=12),
          plot.subtitle = element_text(size=10)) +
    xlab("Chromosome") +
    ylab(expression(-log[10](italic(p)))) +
    ggtitle("100-Seed Weight GWAS Results across Experiments")
dev.off()
