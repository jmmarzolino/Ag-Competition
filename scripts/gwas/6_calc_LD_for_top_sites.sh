#!/bin/bash
#SBATCH --job-name=gwas
#SBATCH --output=/rhome/jmarz001/bigdata/Ag-Competition/scripts/gwas/6_calc_LD_for_top_sites.stdout
#SBATCH --ntasks=14
#SBATCH --mem=100G
#SBATCH -p koeniglab

cd /rhome/jmarz001/bigdata/Ag-Competition/results/gwas

module load plink/1.90b6.25

plink --allow-extra-chr \
--allow-no-sex \
--double-id \
--maf 0.01 \
#--keep AgComp_genotypes.tsv \
--set-missing-var-ids @:#$1,$2 \
--vcf combined_filt.vcf.gz \
--r2 gz --ld-window-r2 0 --ld-window 10000 \
--out LD_10kbwin
