#!/bin/bash
#SBATCH --output=/rhome/jmarz001/bigdata/Ag-Competition/scripts/gwas/5b_clump_sig_regions.stdout
#SBATCH --ntasks=1
#SBATCH --mem=20G
#SBATCH --time=1:00:00
#SBATCH -p short

cd /rhome/jmarz001/bigdata/Ag-Competition/results/gwas
module load plink/1.90b6.25

ls ASSOC_*assoc.tmp > tmp.txt #backup record of file order fed into plink
FILE_LIST=$(ls ASSOC_*assoc.tmp | grep -v / | xargs echo | sed 's/ /,/g')

plink --bfile ../all_traits --double-id --allow-no-sex --allow-extra-chr --clump $FILE_LIST
# clump has two main applications:
#    To report the top X single SNP results from a genome-wide scan in terms of a smaller number of clumps of correlated SNPs (i.e. to assess how many independent loci are associated, for example)
#    To provide a quick way to combine sets of results from two or more studies, when the studies might also be genotyped on different marker sets

#--clump loads the named PLINK-format association report(s) (text files with a header line, a column containing variant IDs, and another column containing p-values) and groups results into LD-based clumps, writing a new report to plink2.clumps[.zst]. Multiple filenames can be separated by spaces or commas.

#By default, PLINK scans these files and extracts fields with the headers SNP and P
