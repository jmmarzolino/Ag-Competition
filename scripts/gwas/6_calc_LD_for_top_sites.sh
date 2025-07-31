#!/bin/bash
#SBATCH --job-name=gwas
#SBATCH --output=/rhome/jmarz001/bigdata/Ag-Competition/scripts/gwas/6_calc_LD_for_top_sites.stdout
#SBATCH --ntasks=10
#SBATCH --mem=100G
#SBATCH --time=02:00:00
#SBATCH -p koeniglab

cd /rhome/jmarz001/bigdata/Ag-Competition/results/gwas

module load plink/1.90b6.25

plink --allow-extra-chr \
--allow-no-sex \
--double-id \
--maf 0.01 \
--set-missing-var-ids @:# \
--vcf combined_filt.vcf.gz \
--ld-snp-list top_sites.txt \
--r2 gz --ld-window-r2 0 --ld-window 1000 --ld-window-kb 5000 \
--out LD_10kbwin

#--set-missing-var-ids @_# 
# @==chr "_" #==position