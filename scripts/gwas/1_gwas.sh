#!/bin/bash
#SBATCH --output=/rhome/jmarz001/bigdata/Ag-Competition/scripts/gwas/1_gwas.stdout
#SBATCH --ntasks=1
#SBATCH --mem=20G
#SBATCH --time=2:00:00
#SBATCH -p short

#cd /rhome/jmarz001/bigdata/Ag-Competition/results
# mkdir gwas
cd /rhome/jmarz001/bigdata/Ag-Competition/results/gwas
#############################################################################
#source files location, alignment by chromosome then stitched back into one vcf. see documentation in BARLEY_CCII_PROGENY_RAD directories
#cp /rhome/dkoenig/bigdata/BARLEY_CCII_PROGENY_RAD/DATA/OUTPUT/STITCH/PROGENY.vcf.gz PROGENY.vcf.gz
#cp /rhome/dkoenig/bigdata/BARLEY_CCII_PROGENY_RAD/DATA/INPUT/SAMPLES.txt SAMPLES.txt
#cp /rhome/dkoenig/bigdata/BARLEY_CCII_PROGENY_RAD/DATA/OUTPUT/STITCH/FULL_FILTER.vcf FULL_FILTER.vcf#parents and progeny lines
#############################################################################

module load vcftools/0.1.16-18
# remove indels for gemma
vcftools --gzvcf PROGENY.vcf.gz --remove-indels --not-chr chrUn --recode --recode-INFO-all --out AG

vcftools --vcf AG.recode.vcf --extract-FORMAT-info GT
head -n1 out.GT.FORMAT | cut -f3- | sed 's/\t/\n/g' | awk '{$1=$1}{print $1" "$1}' > progeny_geno_pheno_list


# create plink format files from vcf (bed, bim, fam)
module load plink/1.90b6.25
plink --vcf AG.recode.vcf --double-id --allow-no-sex --allow-extra-chr --keep progeny_geno_pheno_list -indiv-sort file progeny_geno_pheno_list --make-bed --out all_traits
mv all_traits.log all_traits_bed.log


#### add phenotypes to .fam file
sbatch /rhome/jmarz001/bigdata/Ag-Competition/scripts/gwas/2_prep_phenotypes.R
# and filter vcf file to individuals in phenotype data
cut -d\  -f1 all_traits.fam | awk '{$1=$1}{print $1" "$1}' > common_progeny_geno_pheno_list


plink --vcf AG.recode.vcf --double-id --allow-no-sex --allow-extra-chr --keep common_progeny_geno_pheno_list --pca 20 --out all_traits
# --distance square
# cut only eigenvec values from file; ensure proper formatting of value-tab-value (awk reconstitutes all fields)
cut -d" " -f3-23 all_traits.eigenvec | awk '{OFS="\t"};{$1=$1}{print 1"\t"$0}' > pca.txt


# Relatedness Matrix
module load gemma/0.98.5
/rhome/jmarz001/software/gemma0.98.5 -bfile all_traits -gk 1 -o related_matrix 
#-miss 1 -notsnp


# set number of univariate gwas to run based on
# the number of traits in trait-gwas number file
ARRAY_LIM=$(tail -n +2 trait_name_to_col_numbers.tsv | wc -l | cut -d\  -f1)
# genotype-phenotype association
sbatch --array=1-$ARRAY_LIM%10 /rhome/jmarz001/bigdata/Ag-Competition/scripts/gwas/3_univariate_association_array.sh


sbatch /rhome/jmarz001/bigdata/Ag-Competition/scripts/gwas/3_multivariate_association.sh














#### plot results
sbatch /rhome/jmarz001/bigdata/Ag-Competition/scripts/gwas/5_plot_AllDerivedTraits_GWAS.R

####
# format all_traits.bim
# change 2nd col "." to chr#_pos##
sbatch /rhome/jmarz001/bigdata/Ag-Competition/scripts/gwas/6a_format_chr_pos.R
sbatch /rhome/jmarz001/bigdata/Ag-Competition/scripts/gwas/6b_clump_sig_regions.sh
sbatch --array=1-$ARRAY_LIM /rhome/jmarz001/bigdata/Ag-Competition/scripts/gwas/6b_clump_indv_regions.sh
sbatch /rhome/jmarz001/bigdata/Ag-Competition/scripts/gwas/6c_clump_man.R

####
sbatch /rhome/jmarz001/bigdata/Ag-Competition/scripts/gwas/7_count_sig_regions.R
sbatch /rhome/jmarz001/bigdata/Ag-Competition/scripts/gwas/7_sig_sites_over_traits.R


## plot those gwas results with vertically aligned manhattan plots for each category
sbatch /rhome/jmarz001/bigdata/Ag-Competition/scripts/gwas/7_plot_trait_categories.R
sbatch /rhome/jmarz001/bigdata/Ag-Competition/scripts/gwas/7_plot_trait_categories_sigtraits.R
## plot common gwas sites with upset plots - within and across categories
sbatch /rhome/jmarz001/bigdata/Ag-Competition/scripts/gwas/7_upset.R


# list sites in vcf file
#TAILNUM=$(($(grep -c "##" PROGENY.vcf) + 1))
#cut -f1-2 PROGENY.vcf | tail -n+${TAILNUM} > positions.txt
# do that cleaner w tool~
module load bcftools/1.19
bcftools query -f '%CHROM\t%POS\n' PROGENY.vcf > CALLED_POS.txt

# pull allele counts from vcf
sbatch /rhome/jmarz001/bigdata/Ag-Competition/scripts/gwas/8a_pull_AC.sh

# filter list of vcf sites to common, segregating sites
sbatch /rhome/jmarz001/bigdata/Ag-Competition/scripts/gwas/8b_IPK_segregating_sites_AF.R
# investigate genome sites w significant chang in allele counts
sbatch /rhome/jmarz001/bigdata/Ag-Competition/scripts/gwas/8b_sig_allele_change_sites.R


# plot allele frequency spectra for all sites & for sites identified as significant from gwas
## match starting allele frequencies for sites associated w traits & randomly sample them to create neutral comparison sets
sbatch /rhome/jmarz001/bigdata/Ag-Competition/scripts/gwas/8c_site_sampling.R

sbatch /rhome/jmarz001/bigdata/Ag-Competition/scripts/gwas/8d_plot_AFS.R
sbatch /rhome/jmarz001/bigdata/Ag-Competition/scripts/gwas/8d_plot_deltaAF.R
