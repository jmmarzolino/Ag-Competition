#!/bin/bash
#SBATCH --output=/rhome/jmarz001/bigdata/Ag-Competition/scripts/gwas/greenhouse-comparison/1_gwas_GH.stdout
#SBATCH --ntasks=1
#SBATCH --mem=20G
#SBATCH --time=2:00:00
#SBATCH -p short

cd /rhome/jmarz001/bigdata/Ag-Competition/results/gwas/CCII_greenhouse_exp_gwas
#############################################################################
#source files location

#### set up running gwas w greenhouse experiment 
## mkdir rhome/jmarz001/bigdata/Ag-Competition/results/gwas/CCII_greenhouse_exp_gwas

#cp /rhome/jlandis/bigdata/RADSeq/GEMMA/CCII_Raw_phenotype_data.txt rhome/jmarz001/bigdata/Ag-Competition/results/gwas/CCII_greenhouse_exp_gwas

#### this lower bit is probably not necessary
### potential code for parental vcf files & sites
#cp rhome/dkoenig/bigdata/BARLEY_CCII_PARENTS_RAD/DATA/OUTPUT/VARCALLS/FILTERED.vcf.gz rhome/jmarz001/bigdata/Ag-Competition/results/gwas/PARENTS.vcf.gz
#TMPVCF=PARENTS.vcf.gz
#bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' $TMPVCF > pos_data2.txt
#/rhome/jmarz001/shared/MANUSCRIPT_FILES/2019_LANDIS/CCIISELECT/FINAL_PARENTAL_RAD.vcf.gz
#FINAL_PROGENY_RAD.vcf: The unimputed progeny calls
#############################################################################

#module load vcftools/0.1.16-18
module load bcftools/1.19

# make list of genotypes in both experiments
sub /rhome/jmarz001/bigdata/Ag-Competition/scripts/gwas/greenhouse-comparison/2_joining_exp_phenotypes.stdout

### filter vcf file to indvs used in greenhouse experiment
bcftools view

#CCII_GH_trait_file_nums.tsv  
#CCII_Raw_phenotype_data.txt  
#common_exp.fam  

VCF="../imputed_filter.recode.vcf"
# format indv list for plink
cat exp_common_genos | awk '{$1=$1}{print $1" "$1}' > exp_common_genos2

# create plink format files from vcf (bed, bim, fam)
module load plink/1.90b6.25

plink --allow-extra-chr \
--allow-no-sex \
--double-id \
--maf 0.01 \
--keep exp_common_genos2 \
--indiv-sort f exp_common_genos2 \
--make-bed \
--out all_traits \
--set-missing-var-ids @:#$1,$2 \
--vcf $VCF

### format fam file with phenotypes from greenhouse experiment

#### add phenotypes to .fam file
sbatch /rhome/jmarz001/bigdata/Ag-Competition/scripts/gwas/greenhouse-comparison/2_joining_exp_phenotypes.R

# make pca covar file
plink --allow-extra-chr \
--allow-no-sex \
--double-id \
--keep AgComp_genotypes.tsv \
--vcf imputed_filter.recode.vcf \
--pca 10 \
--out all_traits_pca

# cut only eigenvec values from file; ensure proper formatting of value-tab-value (awk reconstitutes all fields)
cut -d" " -f3- all_traits_pca.eigenvec | awk '{OFS="\t"};{$1=$1}{print 1"\t"$0}' > pca.txt
# awk prints col of "1" \t eigenvecs
# command removes genotype IDs from first two cols, replaces w one col of 1 gemma can use (1 in first col indicates non-missing/included data I believe)


# Relatedness Matrix
module load gemma/0.98.5
/rhome/jmarz001/software/gemma0.98.5 -bfile all_traits -gk 1 -outdir ../output -o related_matrix 
#-miss 1 -notsnp


# GWAS
# set number of univariate gwas to run based on
# the number of traits in trait-gwas number file
ARRAY_LIM=$(tail -n +2 trait_name_to_col_numbers.tsv | wc -l | cut -d\  -f1)
# genotype-phenotype association
sbatch --array=1-$ARRAY_LIM%10 /rhome/jmarz001/bigdata/Ag-Competition/scripts/gwas/3_univariate_association_array.sh

# and gwas with all the traits
#sbatch /rhome/jmarz001/bigdata/Ag-Competition/scripts/gwas/3_multivariate_association.sh 
### proglem w correlated traits, pca? kinship matrix?

# plot results
sbatch /rhome/jmarz001/bigdata/Ag-Competition/scripts/gwas/4_manhattan_plot.R
