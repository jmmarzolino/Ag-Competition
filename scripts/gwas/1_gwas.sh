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
vcftools --gzvcf imputed.vcf.gz --remove-indels --not-chr chrUn --recode --recode-INFO-all --out AG

# list genotypes in raw vcf as basis for plink phenotype file
vcftools --vcf AG.recode.vcf --extract-FORMAT-info GT

# take list of genotypes from first line
# cut cols 3+ (excludes CHR & POS cols)
# sed: search for \t, replace w \n, globally. changes tab-delim col list of genotypes to one col of genos in rows
# awk: copy first col, print col \s col. replicates genotype list into two matching columns for plink
head -n1 out.GT.FORMAT | cut -f3- | sed 's/\t/\n/g' | awk '{$1=$1}{print $1" "$1}' > progeny_geno_pheno_list

# create plink format files from vcf (bed, bim, fam)
module load plink/1.90b6.25
plink --vcf AG.recode.vcf --double-id --allow-no-sex --allow-extra-chr --keep progeny_geno_pheno_list -indiv-sort file progeny_geno_pheno_list --make-bed --out all_traits
mv all_traits.log all_traits_bed.log


#### add phenotypes to .fam file
sbatch /rhome/jmarz001/bigdata/Ag-Competition/scripts/gwas/2_phenotypes.R
# and filter vcf file to individuals in phenotype data
cut -d\  -f1 all_traits.fam | awk '{$1=$1}{print $1" "$1}' > common_progeny_geno_pheno_list


# filter vcf to variant sites in retained progeny
# re-create plink files w filtered progeny list
vcftools --gzvcf imputed.vcf.gz --remove-indels --not-chr chrUn --recode --recode-INFO-all --keep common_progeny_geno_pheno_list --maf 0.002 --out AG
# minor allele freq set to 1/(2*208 seq'd indvs) = 0.0024
# so maf filter will only remove freqs of 0

#Output nucleotide diversity at a list of positions
vcftools --vcf MAF_filt.recode.vcf --freq

# create plink format files from vcf (bed, bim, fam)
plink --vcf AG.recode.vcf --double-id --allow-no-sex --allow-extra-chr --keep common_progeny_geno_pheno_list -indiv-sort file common_progeny_geno_pheno_list --make-bed --out all_traits
mv all_traits.log all_traits_bed2.log

# add phenotypes to .fam file
sbatch /rhome/jmarz001/bigdata/Ag-Competition/scripts/gwas/2_phenotypes.R

# make pca covar file
plink --vcf AG.recode.vcf --double-id --allow-no-sex --allow-extra-chr --keep common_progeny_geno_pheno_list --pca 5 --out all_traits
mv all_traits.log all_traits_pca.log
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
sbatch /rhome/jmarz001/bigdata/Ag-Competition/scripts/gwas/3_multivariate_association.sh

# plot results
sbatch /rhome/jmarz001/bigdata/Ag-Competition/scripts/gwas/4_manhattan_plot.R


## clump gwas results to potentially find sig sites or narrow identified regions

# format all_traits.bim for clumping gwas regions
# change 2nd col "." to chr#_pos##
sbatch /rhome/jmarz001/bigdata/Ag-Competition/scripts/gwas/5a_format_chr_pos.R

# clump gwas results
sbatch /rhome/jmarz001/bigdata/Ag-Competition/scripts/gwas/5b_clump_sig_regions.sh
sbatch --array=1-$ARRAY_LIM /rhome/jmarz001/bigdata/Ag-Competition/scripts/gwas/5b_clump_indv_regions.sh

# plot results
sbatch /rhome/jmarz001/bigdata/Ag-Competition/scripts/gwas/5c_clump_man.R














####
sbatch /rhome/jmarz001/bigdata/Ag-Competition/scripts/gwas/7_count_sig_regions.R
sbatch /rhome/jmarz001/bigdata/Ag-Competition/scripts/gwas/7_sig_sites_over_traits.R


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
