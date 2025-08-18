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

# make list of genotypes in both experiments
sub /rhome/jmarz001/bigdata/Ag-Competition/scripts/gwas/greenhouse-comparison/2a_list_genotypes.R

#CCII_GH_trait_file_nums.tsv  
#CCII_Raw_phenotype_data.txt  
#common_exp.fam  

# create plink format files from vcf (bed, bim, fam)
module load plink/1.90b6.25
VCF="../combined_filt.vcf.gz"

plink --allow-extra-chr \
--allow-no-sex \
--double-id \
--maf 0.01 \
--keep GH_genos \
--indiv-sort f GH_genos \
--make-bed \
--out gh_only \
--set-missing-var-ids @:# \
--vcf $VCF

plink --allow-extra-chr \
--allow-no-sex \
--double-id \
--maf 0.01 \
--keep exp_common_genos \
--indiv-sort f exp_common_genos \
--make-bed \
--out exp_common \
--set-missing-var-ids @:# \
--vcf $VCF

### format fam file with phenotypes from greenhouse experiment

#### add phenotypes to .fam file
sbatch /rhome/jmarz001/bigdata/Ag-Competition/scripts/gwas/greenhouse-comparison/2b_joining_exp_phenotypes.R
sbatch /rhome/jmarz001/bigdata/Ag-Competition/scripts/gwas/greenhouse-comparison/2b_join_GH_phenos.R

# make pca covar file
plink --allow-extra-chr \
--allow-no-sex \
--double-id \
--keep exp_common_genos \
--vcf $VCF \
--pca 10 \
--out exp_common_pca

plink --allow-extra-chr \
--allow-no-sex \
--double-id \
--keep GH_genos \
--vcf $VCF \
--pca 10 \
--out gh_only_pca

# cut only eigenvec values from file; ensure proper formatting of value-tab-value (awk reconstitutes all fields)
cut -d" " -f3- exp_common_pca.eigenvec | awk '{OFS="\t"};{$1=$1}{print 1"\t"$0}' > exp_common_pca.txt
cut -d" " -f3- gh_only_pca.eigenvec | awk '{OFS="\t"};{$1=$1}{print 1"\t"$0}' > gh_only_pca.txt
# awk prints col of "1" \t eigenvecs
# command removes genotype IDs from first two cols, replaces w one col of 1 gemma can use (1 in first col indicates non-missing/included data I believe)


# Relatedness Matrix
module load gemma/0.98.5
/rhome/jmarz001/software/gemma0.98.5 -bfile exp_common -gk 1 -outdir output -o related_matrix 
#-miss 1 -notsnp
/rhome/jmarz001/software/gemma0.98.5 -bfile gh_only -gk 1 -outdir output -o related_matrix_gh_only


# GWAS
# set number of univariate gwas to run based on
# the number of traits in trait-gwas number file
ARRAY_LIM=$(tail -n +2 CCII_GH_trait_file_nums.tsv | wc -l | cut -d\  -f1)
# genotype-phenotype association
sbatch --array=1-$ARRAY_LIM%10 /rhome/jmarz001/bigdata/Ag-Competition/scripts/gwas/greenhouse-comparison/3_univariate_association_array.sh

# and gwas with all the traits
#sbatch /rhome/jmarz001/bigdata/Ag-Competition/scripts/gwas/greenhouse-comparison/3_multivariate_association.sh 
### proglem w correlated traits, pca? kinship matrix?

# plot results
sbatch /rhome/jmarz001/bigdata/Ag-Competition/scripts/gwas/greenhouse-comparison/4_manhattan_plot.R
sbatch /rhome/jmarz001/bigdata/Ag-Competition/scripts/gwas/greenhouse-comparison/4_manhattan_plot_cross_exp_plots.R
##sbatch /rhome/jmarz001/bigdata/Ag-Competition/scripts/gwas/greenhouse-comparison/5_cross_exp_trait_pair_gwas_sites.R
