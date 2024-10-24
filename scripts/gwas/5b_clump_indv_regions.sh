#!/bin/bash
#SBATCH --output=/rhome/jmarz001/bigdata/Ag-Competition/scripts/gwas/5b_clump_indv_regions.stdout
#SBATCH --ntasks=1
#SBATCH --mem=20G
#SBATCH --time=2:00:00
#SBATCH -p short

cd /rhome/jmarz001/bigdata/Ag-Competition/results/gwas
module load plink/1.90b6.25

ASSOC_FILE=$(cut -f3 association_files_traits.txt | tail -n+2 | head -n $SLURM_ARRAY_TASK_ID | tail -n1)
OUT=$(basename $ASSOC_FILE | cut -d. -f1-2)
#T_NUM=$(expr $SLURM_ARRAY_TASK_ID + 5)
ASSOC_TMP=$(echo ""$OUT".tmp")

#OUT=$(basename "$ASSOC_FILE" | cut -d. -f1)
plink --bfile all_traits --double-id --allow-no-sex --allow-extra-chr --clump $ASSOC_TMP --out $OUT
# clump has two main applications:
#    To report the top X single SNP results from a genome-wide scan in terms of a smaller number of clumps of correlated SNPs (i.e. to assess how many independent loci are associated, for example)
#    To provide a quick way to combine sets of results from two or more studies, when the studies might also be genotyped on different marker sets

#--clump loads the named PLINK-format association report(s) (text files with a header line, a column containing variant IDs, and another column containing p-values) and groups results into LD-based clumps, writing a new report to plink2.clumps[.zst]. Multiple filenames can be separated by spaces or commas.

 #By default, PLINK scans these files and extracts fields with the headers SNP and P

 #/rhome/dkoenig/bigdata/BARLEY_CCII_PROGENY_RAD/DATA/OUTPUT/TMP/Barley_Morex_V2_TE_annotation.gff
 #format.gt.pos
