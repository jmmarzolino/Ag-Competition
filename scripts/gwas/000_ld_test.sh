#!/bin/bash
#SBATCH --job-name=gwas
#SBATCH --output=/rhome/jmarz001/bigdata/Ag-Competition/scripts/gwas/000_ld_test.stdout
#SBATCH --ntasks=14
#SBATCH --mem=100G
#SBATCH -p koeniglab

cd /rhome/jmarz001/bigdata/Ag-Competition/results/gwas

module load plink/1.90b6.25

CHR="chr"$SLURM_ARRAY_TASK_ID"H"
OUT="all_traits_chr"$SLURM_ARRAY_TASK_ID

plink --allow-extra-chr \
--allow-no-sex \
--double-id \
--maf 0.01 \
--keep AgComp_genotypes.tsv \
--set-missing-var-ids @:#$1,$2 \
--vcf combined_filt.vcf.gz \
--r2 gz --ld-window 10 --ld-window-kb 10000 --ld-window-r2 0 \
--chr $CHR \
--out $OUT

#--indiv-sort f AgComp_genotypes.tsv \
#--make-bed \
# --distance triange
