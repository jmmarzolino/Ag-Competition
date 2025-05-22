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
