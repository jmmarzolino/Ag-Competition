#!/bin/bash
#SBATCH --job-name=gwas
#SBATCH --output=/rhome/jmarz001/bigdata/Ag-Competition/scripts/gwas/6_calc_LD_for_top_sites.stdout
#SBATCH --ntasks=10
#SBATCH --mem=100G
#SBATCH --time=10:00:00
#SBATCH -p koeniglab

cd /rhome/jmarz001/bigdata/Ag-Competition/results/gwas

module load plink/1.90b6.25

plink --allow-extra-chr \
--allow-no-sex \
--double-id \
--maf 0.01 \
--set-missing-var-ids @:# \
--vcf combined_filt.vcf.gz \
--r2 gz --ld-window-r2 0 --ld-window 100 --ld-window-kb 100 \
--out LD_10kbwin


#--ld-snp chr1H:14499156
#--ld-snp-list top_sites.tsv \
#topsiteLD
#--set-missing-var-ids @_# 
# @==chr "_" #==position