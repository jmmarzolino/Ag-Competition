#!/bin/bash
#SBATCH --output=/rhome/jmarz001/bigdata/Ag-Competition/scripts/gwas/1_gwas.stdout
#SBATCH --ntasks=1
#SBATCH --mem=20G
#SBATCH --time=2:00:00
#SBATCH -p short

cd /rhome/jmarz001/bigdata/Ag-Competition/results/gwas
module load bcftools/1.19

######
cp ~/shared/for_JILL/combinded.vcf.gz ~/bigdata/Ag-Competition/results/gwas/combined.vcf.gz
# list sites in vcf file
bcftools query -f '%CHROM\t%POS\n' combined.vcf.gz > CALLED_POS.txt

# compare number of CC II samples/lines in each vcf
# bcftools query -l PROGENY.vcf.gz | wc -l 
# bcftools query -l imputed.vcf.gz | wc -l 
# bcftools query -l imputed_filter.vcf.gz | wc -l 
# bcftools query -l AG.recode.vcf | wc -l 
# bcftools query -l imputed_filter.recode.vcf.gz | wc -l 
# bcftools query -l PARENTS.vcf.gz | wc -l 
# bcftools query -l combinded.vcf.gz | wc -l 
# bcftools query -l combined_filt.vcf.gz | wc -l 

# remove indels for gemma & filter to variant sites
#tabix -C -p vcf imputed.vcf.gz
#--regions "^"<expression>, leading carot changes inclusion to exclusion


## vcf file already filtered ---
###bcftools_filterCommand=filter -G 3; Date=Thu Dec 14 08:51:17 2023
##bcftools_viewCommand=view -i 'TYPE="snp" & COUNT(GT="het")<4 & MIN(AD[*:0]+AD[*:1])>3 & DP>100 & DP<1000 & QUAL>500 & N_ALT=1'; Date=Thu Dec 14 08:51:17 2023

bcftools view combined.vcf.gz --targets ^chrUn --exclude-types indels -o combined_filt.vcf.gz 
tabix -C -p vcf combined_filt.vcf.gz

#Filter calls for biallelic, no missing calls, minimum coverage 2, less than 6 heterozygotes,
#max depth 500, QUAL more than 30, snps only, at least one homozygous alternate call of GQ>5

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


# pull the genotype names
bcftools query -l combined_filt.vcf.gz  > combined_filt.gt.names
# cut and copy genotype names into two columns (plink format) 
# and remove parental lines (first 28 entries)
tail -n+29 combined_filt.gt.names | awk '{$1=$1}{print $1" "$1}' > progeny_geno_pheno_list

# pull ea. position and all genotypes
# replace genotype indicators (0/0, 0/1, 1/1, ./.) w 0/1/2/NA
bcftools query -f '%CHROM\t%POS[\t%GT]\n' combined_filt.vcf.gz | sed -e s:"0|0":0:g -e s:"0|1":1:g -e s:"1|1":2:g -e s:"\.|\.":NA:g > combined_filter.gt

# check that the number of columns (genotypes) match the number of individuals
wc -l combined_filt.gt.names
head -n1 combined_filter.gt | awk '{print NF}'
# 2 extra fields in combined_filter.gt b/c of chr/pos cols


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


# combine all sig gwas sites across all traits into one table for paper
# make a list of significant snps, top 5%, 1%, and 0.1% of snps from each gwas association file
sbatch /rhome/jmarz001/bigdata/Ag-Competition/scripts/gwas/5_extract_SIGandTOP_sites_across_traits.R



## clump gwas results to potentially find sig sites or narrow identified regions

# format all_traits.bim for clumping gwas regions
# change 2nd col "." to chr#_pos##
sbatch /rhome/jmarz001/bigdata/Ag-Competition/scripts/gwas/5a_format_chr_pos.R


# clump gwas results
ls ASSOC_*.assoc.txt > all_assoc_files.txt
sbatch /rhome/jmarz001/bigdata/Ag-Competition/scripts/gwas/5b_clump_sig_regions.sh
sbatch /rhome/jmarz001/bigdata/Ag-Competition/scripts/gwas/5c_clumped_sites_manhattan.R

ARRAY_LIM=$(wc -l all_assoc_files.txt | cut -d\  -f1)
sbatch --array=1-$ARRAY_LIM /rhome/jmarz001/bigdata/Ag-Competition/scripts/gwas/5b_clump_indv_regions.sh

# plot results
sbatch /rhome/jmarz001/bigdata/Ag-Competition/scripts/gwas/5c_clumped_indv_sites_manhattan.R


## calculate LD
## plot LD decay around lead snps

# plot chr 4 peak and vrn h2
# zoom in view of chr 4 & 5 regions
# snp effect of lead snps

## FT allele freq over time, using only main snps

# pull the full gwas info (position, p-val, beta) for site identified in clumping
sbatch /rhome/jmarz001/bigdata/Ag-Competition/scripts/gwas/6_format_clumped_sites_for_AFchange.R



## calculate LD for all sites
sbatch /rhome/jmarz001/bigdata/Ag-Competition/scripts/gwas/6_calc_LD_for_top_sites.sh

# filter the LD output to r2 vals of top sites
# extract the SNP list from all 3 clumped files and copy them into a new, one col list/file
awk '{print $3}' ASSOC_6_lmm.assoc.clumped | head -n-2 > top_sites.txt
awk '{print $3}' ASSOC_7_lmm.assoc.clumped | tail -n+2 | head -n-2 >> top_sites.txt
awk '{print $3}' ASSOC_8_lmm.assoc.clumped | tail -n+2 | head -n-2 >> top_sites.txt

## plot LD decay around peak snps
sbatch /rhome/jmarz001/bigdata/Ag-Competition/scripts/gwas/000_LD_plot.R





### ALLELE FREQUENCIES
# extract allele counts for all sites from progeny sequencing
sbatch /rhome/jmarz001/bigdata/Ag-Competition/scripts/gwas/6_pull_AC.sh

sbatch /rhome/jmarz001/bigdata/Ag-Competition/scripts/gwas/6b_polarize_sites_to_gwas_beta.R

### test genome sites for significant changes in allele frequency between generations (0-18, 18-58)
# actually need to troubleshoot it for 0-18 period...
# and plot distributions of allele change btwn generations

# test site allele counts for significant differences between generations F18 and F58; calculate allele frequency change for those sites; plot allele frequency spectra & distribution of change in allele frequency; list which allele (A1, A2) is beneficial
sbatch /rhome/jmarz001/bigdata/Ag-Competition/scripts/gwas/6b_test_for_sig_allele_change_sites.R

## filter the genotype file before reading into R
## by limiting the number of sites
# start by formatting your list of positions for bcftools
awk '{print $1"\t"$3}' gwas_top_sites.tsv | tail -n+2 > top_sites.tsv


bcftools view --regions-file top_sites.tsv combined_filt.vcf.gz -o top_sites.vcf
bcftools query -f '%CHROM\t%POS[\t%GT]\n' top_sites.vcf | sed -e s:"0|0":0:g -e s:"0|1":1:g -e s:"1|0":1:g -e s:"1|1":2:g -e s:"\.|\.":NA:g > top_sites.gt

sbatch /rhome/jmarz001/bigdata/Ag-Competition/scripts/gwas/6c_plot_top_sites_allele_effect.R




############

# plot allele frequency spectra for all sites & for sites identified as significant from gwas
## match starting allele frequencies for sites associated w traits & randomly sample them to create neutral comparison sets
sbatch /rhome/jmarz001/bigdata/Ag-Competition/scripts/gwas/7a_site_sampling.R

sbatch /rhome/jmarz001/bigdata/Ag-Competition/scripts/gwas/7b_plot_AFS.R
sbatch /rhome/jmarz001/bigdata/Ag-Competition/scripts/gwas/7b_plot_deltaAF.R






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
