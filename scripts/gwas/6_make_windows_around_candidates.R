#!/usr/bin/env Rscript

#SBATCH --job-name=manhattan
#SBATCH -o /rhome/jmarz001/bigdata/Ag-Competition/scripts/gwas/6_make_windows_around_candidates.stdout
#SBATCH --ntasks=1
#SBATCH --mem=20gb
#SBATCH -t 00:30:00
#SBATCH -p short

setwd("/bigdata/koeniglab/jmarz001/Ag-Competition/results/gwas")
library(pacman)
p_load(tidyverse, data.table, ggsci, Cairo, qqman)

# read in file
ft <- fread("ASSOC_6_lmm.assoc.clumped")
sw <- fread("ASSOC_7_lmm.assoc.clumped")
fec <- fread("ASSOC_8_lmm.assoc.clumped")




## make lead snps list from lmm gwas & clumped file
## set windows around those snps (100 kb)
# - make list of top 100-1000 sites + LD prune them
## calculate LD
## plot LD decay around lead snps

# plot chr 4 peak and vrn h2
# zoom in view of chr 4 & 5 regions
# snp effect of lead snps
# FT allele freq over time, using only main snps


# sort sites by chr/position
ft2 <- ft[,1:11]
ft2 <- ft2[order(ft2$CHR, ft2$BP),]
table(ft2$CHR)
ft2$CHR <- as.numeric(gsub("chr(\\d)H", "\\1", ft2$CHR))


ft0 <- fread("ASSOC_6_lmm.assoc.txt")
ft0$chr <- as.numeric(gsub("chr(\\d)H", "\\1", ft0$chr))

#pdf("ft_clump_snps_highlighted_manhattan.png")
#manhattan(ft0, chr="chr", bp="ps", snp="rs", p="p_lrt", highlight=ft2$SNP)
#dev.off()


for(i in 1:5){
  df <- fread(paste0("plink_clump_lmm_r0.", i, ".clumped"))
  highs <- df$SNP

  png(paste0("clump_test_r", i, ".png"))

  manhattan(ft0, chr="chr", bp="ps", snp="rs", p="p_lrt", highlight=highs)

  dev.off()
}



for(i in 1:6){
  df <- fread(paste0("plink_clump_lmm_", i, "kb.clumped"))
  highs <- df$SNP

  png(paste0("clump_test_", i, ".png"))

  manhattan(ft0, chr="chr", bp="ps", snp="rs", p="p_lrt", highlight=highs)

  dev.off()
}





## 
manhattan(subset(ft0, chr==5), chr="chr", bp="ps", snp="rs", p="p_lrt", highlight=ft2$SNP)
manhattan(subset(ft0, chr==5), chr="chr", bp="ps", snp="rs", p="p_lrt", highlight=ft2$SNP, xlim=c(450000000, 550000000))

manhattan(subset(ft0, chr==4), chr="chr", bp="ps", snp="rs", p="p_lrt", highlight=ft2$SNP, xlim=c(588000000, 592000000))

annotatePval = 0.005, annotateTop = FALSE

for(i in 1:7){
    # set data to run through each chr subset in turn
    #chr <- paste0("chr", i, "H")
    #tmp <- ft2 %>% filter(CHR==chr)

    # find lowest p-value (its already at the top)
    min(tmp$P)
}

# manhattan function requires numeric chromosome col
chr1 <- ft2[which(ft2$CHR=="chr1H"),]
chr1$CHR <- as.numeric(gsub("chr(\\d)H", "\\1", chr1$CHR))
manhattan(chr1)

manhattan(
       x,
       chr = "CHR",
       bp = "BP",
       p = "P",
       snp = "SNP",
       col = c("gray10", "gray60"),
       chrlabs = NULL,
       suggestiveline = -log10(1e-05),
       genomewideline = -log10(5e-08),
       highlight = NULL,
       logp = TRUE,
       annotatePval = NULL,
       annotateTop = TRUE,
       ...
     )



# take first peak area and make 100kb window around it
# check if other sites in list are within that region, if so, skip it as it's already included in that window
# cont. through file & output file of chr /t peak position /t window start /t window end