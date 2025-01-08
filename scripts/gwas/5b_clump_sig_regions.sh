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

plink --bfile all_traits --double-id --allow-no-sex --allow-extra-chr --clump $FILE_LIST
