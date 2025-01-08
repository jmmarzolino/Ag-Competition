#!/usr/bin/env Rscript
#SBATCH --job-name=gwas
#SBATCH -o /rhome/jmarz001/bigdata/Ag-Competition/scripts/gwas/5a_format_chr_pos.stdout
#SBATCH --mem=20G
#SBATCH -t 00:10:00
#SBATCH -p short

library(pacman)
p_load(tidyverse, data.table, gridExtra, ggsci, Cairo)

options(stringsAsFactors = F)

setwd("/rhome/jmarz001/bigdata/Ag-Competition/results/gwas")

### Fill/Update variant ID column (col 2) with the combination of chr and position
# read in bim file
bim <- read_delim("all_traits.bim", "\t", col_names=F)
# make 2nd col a combo of chr and position
bim$X2 <- paste0(bim$X1,"_", bim$X4)
# write bim file out
write_delim(bim, "all_traits.bim", "\t", col_names=F)


############################################################
############################################################
### Read gwas results (ASSOC files), extract 4 columns, write out

# read list of assoc files to run over
assoc_files_list <- read_delim("trait_name_to_col_numbers.tsv")

assoc_files_list <- assoc_files_list[order(as.character(assoc_files_list$trait_num)),]
assoc_files_list$clump_order <- 1:nrow(assoc_files_list)

write_delim(assoc_files_list, "association_files_traits.txt", "\t", col_names=T)

# write out each SNP association table formatted for plink to clump snp positions

format_assoc <- function(filename){
  # import data
  df <- fread(filename)

  # make a table w the columns CHR, SNP (variant ID), BP (position on chr), P (association p-value)
  x <- tibble(CHR=df$chr, SNP=paste(df$chr, df$ps, sep="_"), BP=df$ps, P=df$p_lrt)
  # p_wald, p_lrt, p_score

  filename_out <- gsub("txt", "tmp", filename)
  write_delim(x, filename_out, "\t", col_names=T)
  #return(x)
}


# set variables
lst <- assoc_files_list$file
lapply(lst, format_assoc)
