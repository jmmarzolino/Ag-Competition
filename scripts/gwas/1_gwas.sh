#!/bin/bash
#SBATCH --output=/rhome/jmarz001/bigdata/Ag-Competition/scripts/gwas/1_gwas.stdout
#SBATCH --ntasks=1
#SBATCH --mem=20G
#SBATCH --time=2:00:00
#SBATCH -p short

cd /rhome/jmarz001/bigdata/Ag-Competition/results/gwas
module load bcftools/1.19

######
cp ~/shared/for_JILL/combined.vcf.gz ~/bigdata/Ag-Competition/results/gwas
# list sites in vcf file
bcftools query -f '%CHROM\t%POS\n' PROGENY.vcf > CALLED_POS.txt

# remove indels for gemma & filter to variant sites
#tabix -C -p vcf imputed.vcf.gz
#bcftools view imputed.vcf.gz --targets ^chrUn --exclude-types indels --exclude "MAF>0.1" -o imputed_filter.vcf 
#--regions "^"<expression>, leading carot changes inclusion to exclusion

bcftools view combinded.vcf.gz --targets ^chrUn --exclude-types indels --exclude "MAF>0.1" -o combined_filt.vcf.gz 

# Unlike -r, targets can be prefixed with "^"
#bgzip imputed_filter.vcf
#tabix -C -p vcf imputed_filter.vcf.gz
# bcftools recognizes MAF and claculates it on the fly
#frequency of minor alleles (MAF=MAC/AN)
#exclude sites for which EXPRESSION is true.
# MAF expression is better than min-af, because it does actually look for *minor* allele freqs
#--write-index # use for bgzf compressed files

# list genotypes in raw vcf as basis for plink phenotype file
#vcftools --vcf imputed_filter.recode.vcf --extract-FORMAT-info GT
# take list of genotypes from first line
# cut cols 3+ (excludes CHR & POS cols)
# sed: search for \t, replace w \n, globally. changes tab-delim col list of genotypes to one col of genos in rows
# awk: copy first col, print col \s col. replicates genotype list into two matching columns for plink
#head -n1 out.GT.FORMAT | cut -f3- | sed 's/\t/\n/g' | awk '{$1=$1}{print $1" "$1}' > progeny_geno_pheno_list

#
#bcftools query -l imputed_filter.vcf.gz > imputed_filter.gt.names
#cat imputed_filter.gt.names | awk '{$1=$1}{print $1" "$1}' > progeny_geno_pheno_list

bcftools query -f '[\t%GT]\n' combined_filt.vcf.gz | sed -e s:"0/0":0:g -e s:"0/1":1:g -e s:"1/1":2:g -e s:"\./\.":NA:g > combined_filter.gt

# create plink format files from vcf (bed, bim, fam)
module load plink/1.90b6.25
plink --allow-extra-chr \
--allow-no-sex \
--double-id \
--maf 0.01 \
--keep AgComp_genotypes.tsv \
--indiv-sort f AgComp_genotypes.tsv \
--make-bed \
--out all_traits \
--set-missing-var-ids @:#$1,$2 \
--vcf combined_filt.vcf.gz


#### add phenotypes to .fam file
sbatch /rhome/jmarz001/bigdata/Ag-Competition/scripts/gwas/2_phenotypes.R

# make pca covar file
plink --allow-extra-chr \
--allow-no-sex \
--double-id \
--keep AgComp_genotypes.tsv \
--vcf combined_filt.vcf.gz \
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


## clump gwas results to potentially find sig sites or narrow identified regions

# format all_traits.bim for clumping gwas regions
# change 2nd col "." to chr#_pos##
sbatch /rhome/jmarz001/bigdata/Ag-Competition/scripts/gwas/5a_format_chr_pos.R

# combine all sig gwas sites across all traits into one table for paper



# clump gwas results
ls ASSOC_*.assoc.txt > all_assoc_files.txt
sbatch /rhome/jmarz001/bigdata/Ag-Competition/scripts/gwas/5b_clump_sig_regions.sh
sbatch /rhome/jmarz001/bigdata/Ag-Competition/scripts/gwas/5c_clumped_sites_manhattan.R

ARRAY_LIM=$(wc -l all_assoc_files.txt | cut -d\  -f1)
sbatch --array=1-$ARRAY_LIM /rhome/jmarz001/bigdata/Ag-Competition/scripts/gwas/5b_clump_indv_regions.sh

# plot results
sbatch /rhome/jmarz001/bigdata/Ag-Competition/scripts/gwas/5c_clumped_indv_sites_manhattan.R


# make a list of significant snps, top 5%, 1%, and 0.1% of snps from each gwas association file
sbatch /rhome/jmarz001/bigdata/Ag-Competition/scripts/gwas/5d_extract_leading_snps.R

# extract allele counts for all sites from progeny sequencing
sbatch /rhome/jmarz001/bigdata/Ag-Competition/scripts/gwas/6_pull_AC.sh

# calculate allele frequncy for each generation based on allele counts; filter non-segregating sites
sbatch /rhome/jmarz001/bigdata/Ag-Competition/scripts/gwas/6b_segregating_sites_AF.R

# test site allele counts for significant differences between generations F18 and F58; calculate allele frequency change for those sites; plot allele frequency spectra & distribution of change in allele frequency; list which allele (A1, A2) is beneficial
sbatch /rhome/jmarz001/bigdata/Ag-Competition/scripts/gwas/6b_sig_allele_change_sites.R

# 
sbatch /rhome/jmarz001/bigdata/Ag-Competition/scripts/gwas/6c_plot_top_sites_allele_effect.R





############
# submit controlling script for greenhouse experiment gwas
sbatch /rhome/jmarz001/bigdata/Ag-Competition/scripts/gwas/greenhouse-comparison/1_gwas_GH.sh










####
sbatch /rhome/jmarz001/bigdata/Ag-Competition/scripts/gwas/7_count_sig_regions.R
sbatch /rhome/jmarz001/bigdata/Ag-Competition/scripts/gwas/7_sig_sites_over_traits.R



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
