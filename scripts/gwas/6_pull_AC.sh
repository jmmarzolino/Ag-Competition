#!/bin/bash
#SBATCH --job-name=GWAS
#SBATCH -o /rhome/jmarz001/bigdata/Ag-Competition/scripts/gwas/6_pull_AC.stdout
#SBATCH --ntasks=1
#SBATCH --mem=50gb
#SBATCH -t 10:00:00
#SBATCH -p koeniglab


# calculate generations' minor allele frequencies
cd /rhome/jmarz001/bigdata/Ag-Competition/results/gwas
#module load vcftools/0.1.16-18
module load bcftools/1.19

VCF='imputed_filter.recode.vcf'
#VCF2='/rhome/jmarz001/bigdata/IPK_Analysis/results/GWAS/PARENTS.vcf.gz'
#vcftools --gzvcf $VCF2 --remove-indels --not-chr chrUn --positions CALLED_POS.txt --recode --recode-INFO-all --out PARENT_SITEFILT
#VCF2='/rhome/jmarz001/bigdata/IPK_Analysis/results/GWAS/PARENT_SITEFILT.recode.vcf'

### ALLELE FREQUENCIES
# bcftools needs an index
#bgzip $VCF
VCF='imputed_filter.recode.vcf.gz'
#bcftools index $VCF

# bcftools query
# Extracts fields from VCF or BCF files and outputs them in user-defined format.
## list generations' genotype IDs
#bcftools query -l $VCF2 > PARENTS.gt.names
bcftools query -l $VCF | grep '^1_' > F18.gt.names
bcftools query -l $VCF | grep '^2_' > F28.gt.names
bcftools query -l $VCF | grep '^3_' > F50.gt.names
bcftools query -l $VCF | grep '^7_' > F58.gt.names

### calculate MAF for generations
# list chr-position data with ref and alt alelles for both site lists
#bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' $VCF2 > pos_data.txt
## column order is REF count / ALT count
# list chr-position data for gwas sites
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' $VCF > prog_sites.txt


# script counts occurance of alt and ref alleles per line/site
PY_SCRIPT=$'../../scripts/gwas/calc_afs_from_vcf.py'
# extract allele counts for both site lists
#bcftools query -f '[\t%GT]\n' $VCF2 | python $PY_SCRIPT > AFS1.txt

bcftools query -S F18.gt.names -f '[\t%GT]\n' $VCF | python $PY_SCRIPT > AFS18.txt
bcftools query -S F28.gt.names -f '[\t%GT]\n' $VCF | python $PY_SCRIPT > AFS28.txt
bcftools query -S F50.gt.names -f '[\t%GT]\n' $VCF | python $PY_SCRIPT > AFS50.txt
bcftools query -S F58.gt.names -f '[\t%GT]\n' $VCF | python $PY_SCRIPT > AFS58.txt

#paste pos_data.txt AFS1.txt AFS18.txt AFS28.txt AFS50.txt AFS58.txt > COMBINED_AFS.txt
paste prog_sites.txt AFS18.txt AFS28.txt AFS50.txt AFS58.txt > COMBINED_AFS.txt

#rm F* AFS* all_sites.txt pos_data.txt