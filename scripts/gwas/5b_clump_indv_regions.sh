#!/bin/bash
#SBATCH --output=/rhome/jmarz001/bigdata/Ag-Competition/scripts/gwas/5b_clump_indv_regions.stdout
#SBATCH --ntasks=1
#SBATCH --mem=20G
#SBATCH --time=2:00:00
#SBATCH -p short

cd /rhome/jmarz001/bigdata/Ag-Competition/results/gwas
module load plink/1.90b6.25

ASSOC_FILE=$(cut -f1 all_assoc_files.txt | head -n $SLURM_ARRAY_TASK_ID | tail -n1)
OUT=$(basename $ASSOC_FILE | cut -d. -f1-2)
ASSOC_TMP=$(echo ""$OUT".tmp")

plink --bfile all_traits --double-id --allow-no-sex --allow-extra-chr --clump $ASSOC_TMP --out $OUT




FILE_LIST_LMM=ASSOC_6_lmm.assoc.tmp
# test different window sizes
# clump based on list of top p-value sites
plink --bfile all_traits --double-id --allow-no-sex --allow-extra-chr --clump $FILE_LIST_LMM --clump-kb 4500 --clump-r2 0.1 --out plink_clump_lmm_1

plink --bfile all_traits --double-id --allow-no-sex --allow-extra-chr --clump $FILE_LIST_LMM --clump-kb 5000 --clump-r2 0.1 --out plink_clump_lmm_2

plink --bfile all_traits --double-id --allow-no-sex --allow-extra-chr --clump $FILE_LIST_LMM --clump-kb 5500 --clump-r2 0.1 --out plink_clump_lmm_3

plink --bfile all_traits --double-id --allow-no-sex --allow-extra-chr --clump $FILE_LIST_LMM --clump-kb 6000 --clump-r2 0.1 --out plink_clump_lmm_4

plink --bfile all_traits --double-id --allow-no-sex --allow-extra-chr --clump $FILE_LIST_LMM --clump-kb 7000 --clump-r2 0.1 --out plink_clump_lmm_5

plink --bfile all_traits --double-id --allow-no-sex --allow-extra-chr --clump $FILE_LIST_LMM --clump-kb 8000 --clump-r2 0.1 --out plink_clump_lmm_6

plink --bfile all_traits --double-id --allow-no-sex --allow-extra-chr --clump $FILE_LIST_LMM --clump-kb 9000 --clump-r2 0.1 --out plink_clump_lmm_7

plink --bfile all_traits --double-id --allow-no-sex --allow-extra-chr --clump $FILE_LIST_LMM --clump-kb 10000 --clump-r2 0.1 --out plink_clump_lmm_8
# 10 Mb
plink --bfile all_traits --double-id --allow-no-sex --allow-extra-chr --clump $FILE_LIST_LMM --clump-kb 11000 --clump-r2 0.1 --out plink_clump_lmm_9

plink --bfile all_traits --double-id --allow-no-sex --allow-extra-chr --clump $FILE_LIST_LMM --clump-kb 12000 --clump-r2 0.1 --out plink_clump_lmm_10

plink --bfile all_traits --double-id --allow-no-sex --allow-extra-chr --clump $FILE_LIST_LMM --clump-kb 13000 --clump-r2 0.1 --out plink_clump_lmm_11

plink --bfile all_traits --double-id --allow-no-sex --allow-extra-chr --clump $FILE_LIST_LMM --clump-kb 15000 --clump-r2 0.1 --out plink_clump_lmm_12



plink --bfile all_traits --double-id --allow-no-sex --allow-extra-chr --clump $FILE_LIST_LMM --clump-kb 20000 --clump-r2 0.1 --out plink_clump_lmm_13

wc -l plink_clump_lmm_*clumped

plink --bfile all_traits --double-id --allow-no-sex --allow-extra-chr --clump $FILE_LIST_LMM --clump-kb 35000 --clump-r2 0.1 --out plink_clump_lmm_14
