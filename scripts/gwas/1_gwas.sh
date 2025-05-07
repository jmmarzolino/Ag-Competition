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
#/rhome/dkoenig/bigdata/BARLEY_CCII_PROGENY_RAD/SCRIPTS
#############################################################################

module load vcftools/0.1.16-18
module load bcftools/1.19




# remove indels for gemma & filter to variant sites
vcftools --gzvcf imputed.vcf.gz --remove-indels --not-chr chrUn --maf 0.1 --recode --recode-INFO-all --out imputed_filter
# maf filter to remove non-varying sites
#After filtering, kept 155221 out of a possible 425098 Sites
tabix -C -p vcf imputed.vcf.gz
bcftools view imputed.vcf.gz --exclude-types indels --exclude "MAF>0.1" -o output
--min-af FLOAT
# bcftools recognizes MAF and claculates it on the fly
#frequency of minor alleles (MAF=MAC/AN)
bcftools view imputed.vcf.gz --exclude "MAF>0.1" -o output_MAF
#exclude sites for which EXPRESSION is true.
# MAF expression is better than min-af, because it does actually look for *minor* allele freqs
# min-af just looks at allele count / allele number
bcftools view imputed.vcf.gz --exclude-types indels -o output_indels
bcftools view imputed.vcf.gz --min-af 0.1 -o output_min-af


--write-index


# bcftools view -S <(awk '{if ($5<.6) print $1}' NUM_MISS.imiss) -i 'F_MISSING<.4' FINAL_PROGENY_RAD.vcf.gz > FINAL_FILTER_RAD.vcf
# /rhome/dkoenig/bigdata/LOCAL_SOFTWARE/htslib/htslib/bgzip FINAL_FILTER_RAD.vcf
# /rhome/dkoenig/bigdata/LOCAL_SOFTWARE/htslib/htslib/tabix -C -p vcf FINAL_FILTER_RAD.vcf.gz
# bcftools query -f '%CHROM\t%POS\n' FINAL_FILTER_RAD.vcf.gz > RAD_SITES.txt
# /rhome/dkoenig/bigdata/LOCAL_SOFTWARE/htslib/htslib/tabix -h -R RAD_SITES.txt ~/bigdata/BARLEY_CCII_PARENTS_RAD/DATA/OUTPUT/VARCALLS/FILTERED.vcf.gz > PARENTS.vcf
# /rhome/dkoenig/bigdata/LOCAL_SOFTWARE/htslib/htslib/bgzip PARENTS.vcf
# /rhome/dkoenig/bigdata/LOCAL_SOFTWARE/htslib/htslib/tabix -C -p vcf PARENTS.vcf.gz
# bcftools merge PARENTS.vcf.gz FINAL_FILTER_RAD.vcf.gz > PAR_PROG_FILTER_RAD.vcf
# /rhome/dkoenig/bigdata/LOCAL_SOFTWARE/htslib/htslib/bgzip PAR_PROG_FILTER_RAD.vcf
# /rhome/dkoenig/bigdata/LOCAL_SOFTWARE/htslib/htslib/tabix -C -p vcf PAR_PROG_FILTER_RAD.vcf.gz



# list genotypes in raw vcf as basis for plink phenotype file
vcftools --vcf imputed_filter.recode.vcf --extract-FORMAT-info GT
vcftools --gzvcf imputed_filter.recode.vcf.gz --remove-indels --not-chr chrUn --maf 0.1 --kept-sites
vcftools --vcf FINAL_PROGENY_RAD_INDFILT.vcf --out missing --missing-site
#Output nucleotide diversity at a list of positions
#vcftools --vcf imputed_filter.recode.vcf --freq
# vcftools --vcf FINAL_PROGENY_RAD.vcf --missing-indv --out NUM_MISS



bcftools query -f '%CHROM\t%POS\n' FULL_FILTER.vcf > FULL_FILTER.gt.pos
bcftools query -l FULL_FILTER.vcf > FULL_FILTER.gt.names
bcftools query -f '[\t%GT]\n' FULL_FILTER.vcf | sed -e s:"0/0":0:g -e s:"0/1":1:g -e s:"1/1":2:g -e s:"\./\.":NA:g > FULL_FILTER.gt
#
bcftools view -S samp_.6.txt FINAL_PROGENY_RAD.vcf > FINAL_PROGENY_RAD_INDFILT.vcf

bcftools view -i 'N_MISSING<273' FINAL_PROGENY_RAD_INDFILT.vcf > FINAL_PROGENY_RAD_FILTER_0.35.vcf

bcftools query -f '%CHROM\t%POS\n' FINAL_PROGENY_RAD_FILTER_0.35.vcf > FINALPOS_0.35.vcf
sed s/"\/"/"\|"/g FILTERED.vcf | bcftools view -T FINALPOS_0.35.vcf > PARENTS_PHASED_0.35.vcf

cat <(bcftools view -h FINAL_PROGENY_RAD_FILTER_0.35.vcf)  <(bcftools view -H FINAL_PROGENY_RAD_FILTER_0.35.vcf | awk '{OFS="\t"};{$3=$1"."$2"."$4"."$5; print $0}') > PROGENY.vcf
cat <(bcftools view -h PARENTS_PHASED_0.35.vcf)  <(bcftools view -H PARENTS_PHASED_0.35.vcf | awk '{OFS="\t"};{$3=$1"."$2"."$4"."$5; print $0}') > PARENTS.vcf



mkdir ~/bigdata/BARLEY_CCII_PROGENY_RAD/DATA/OUTPUT/AFS
cd ~/bigdata/BARLEY_CCII_PROGENY_RAD/DATA/OUTPUT/AFS
#
bcftools query -l ../VARCALLS/PROGENY.vcf.gz | grep '^1_' > F18.txt
bcftools query -l ../VARCALLS/PROGENY.vcf.gz | grep '^2_' > F28.txt
bcftools query -l ../VARCALLS/PROGENY.vcf.gz | grep '^3_' > F50.txt
bcftools query -l ../VARCALLS/PROGENY.vcf.gz | grep '^7_' > F58.txt
#
bcftools query -f '%CHROM\t%POS\n' ../VARCALLS/PROGENY.vcf.gz > CALLED_POS.txt
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' ../VARCALLS/PROGENY.vcf.gz > pos_data.txt
bcftools query -R CALLED_POS.txt -f '[\t%GT]\n' ../VARCALLS/PARENTS.vcf.gz | python ../../../SCRIPTS/calc_afs_from_vcf.py > AFS1.txt
bcftools query -S F18.txt -f '[\t%GT]\n' ../VARCALLS/PROGENY.vcf.gz | python ../../../SCRIPTS/calc_afs_from_vcf.py > AFS18.txt
bcftools query -S F28.txt -f '[\t%GT]\n' ../VARCALLS/PROGENY.vcf.gz | python ../../../SCRIPTS/calc_afs_from_vcf.py > AFS28.txt
bcftools query -S F50.txt -f '[\t%GT]\n' ../VARCALLS/PROGENY.vcf.gz | python ../../../SCRIPTS/calc_afs_from_vcf.py > AFS50.txt
bcftools query -S F58.txt -f '[\t%GT]\n' ../VARCALLS/PROGENY.vcf.gz | python ../../../SCRIPTS/calc_afs_from_vcf.py > AFS58.txt
paste pos_data.txt AFS1.txt AFS18.txt AFS28.txt AFS50.txt AFS58.txt > COMBINDED_AFS.txt
rm F* AFS* CALLED_POS.txt pos_data.txtmodule load plink


bcftools query -f '%CHROM\t%POS\t%CHROM\_%POS\_%REF\_%ALT\t%REF\t%ALT\t%QUAL\t%FILTER\t\.\t\GT[\t%GT]\n' RAD_CALLS_FULL.vcf | sed s/'\/'/'|'/g > RAD_CALLS_FULL_PHASE.vcf
bcftools query -l RAD_CALLS_FULL.vcf > parent_calls.txt
#Create vcf
cd ~/bigdata/BARLEY_CCII_PROGENY_RAD/DATA/OUTPUT/COMBINED_VARIANTS/
paste ../CALL_VARIANTS/* > FULLVAR_CALLS.txt
paste ~/bigdata/BARLEY_CCII_PARENTS_RAD/DATA/OUTPUT/VARCALLS/RAD_CALLS_FULL_PHASE.vcf FULLVAR_CALLS.txt > combined_Par_Prog_calls_noheader.vcf
cat ~/bigdata/BARLEY_CCII_PARENTS_RAD/DATA/OUTPUT/VARCALLS/parent_calls.txt <(ls ../CALL_VARIANTS) | sed s/'.bam.calls'//g > sample_genos.txt

bcftools query -l filtered_Par_Prog.vcf > s_names.txt
bcftools query -f '[\t%GT]\n' filtered_Par_Prog.vcf > temp.vcf
bcftools view -S lt_50.txt -i 'COUNT(GT="het")<20 && COUNT(GT="ref")>0 && COUNT(GT="alt")>0' filtered_Par_Prog.vcf > FINAL.vcf
bcftools view -S fpar_names.txt FINAL.vcf > PARENTS.vcf
bcftools view -S fprog_names.txt FINAL.vcf > PROGENY.vcf

bcftools view -S <(awk '{if ($5<.6) print $1}' NUM_MISS.imiss) -i 'F_MISSING<.4' FINAL_PROGENY_RAD.vcf.gz > FINAL_FILTER_RAD.vcf
bcftools query -f '%CHROM\t%POS\n' FINAL_FILTER_RAD.vcf.gz > RAD_SITES.txt


mkdir /rhome/dkoenig/bigdata/BARLEY_CCII_PROGENY_RAD/DATA/OUTPUT/AFS
cd /rhome/dkoenig/bigdata/BARLEY_CCII_PROGENY_RAD/DATA/OUTPUT/AFS
#
bcftools query -l ../STITCH/PROGENY.vcf.gz | grep "^1_" > F18.txt
bcftools query -l ../STITCH/PROGENY.vcf.gz | grep "^2_" > F28.txt
bcftools query -l ../STITCH/PROGENY.vcf.gz | grep "^3_" > F50.txt
bcftools query -l ../STITCH/PROGENY.vcf.gz | grep "^7_" > F58.txt
#
vcftools --gzvcf  FILTERED.vcf.gz --counts2 --out F0
vcftools --gzvcf ../STITCH/PROGENY.vcf.gz --keep F18.txt --counts2 --out F18
vcftools --gzvcf ../STITCH/PROGENY.vcf.gz --keep F28.txt --counts2 --out F28
vcftools --gzvcf ../STITCH/PROGENY.vcf.gz --keep F50.txt --counts2 --out F50
vcftools --gzvcf ../STITCH/PROGENY.vcf.gz --keep F58.txt --counts2 --out F58





# take list of genotypes from first line
# cut cols 3+ (excludes CHR & POS cols)
# sed: search for \t, replace w \n, globally. changes tab-delim col list of genotypes to one col of genos in rows
# awk: copy first col, print col \s col. replicates genotype list into two matching columns for plink
head -n1 out.GT.FORMAT | cut -f3- | sed 's/\t/\n/g' | awk '{$1=$1}{print $1" "$1}' > progeny_geno_pheno_list

# create plink format files from vcf (bed, bim, fam)
module load plink/1.90b6.25
plink --allow-extra-chr \
--allow-no-sex \
--double-id \
--maf 0.01 \
--keep progeny_geno_pheno_list \
--indiv-sort f progeny_geno_pheno_list \
--make-bed \
--out all_traits \
--set-missing-var-ids @:#$1,$2 \
--vcf imputed_filter.recode.vcf

#~/bigdata/TMP/TEST_JILL_GWAS
# not previously included
#--maf .01 \ 
#--set-missing-var-ids @:#$1,$2 \
#--indiv-sort f riverside.fam \ #previous syntax/args were slightly diff, ie. -indv-sort file vs. --indv-sort f
#--vcf imputed_filter.vcf.gz # unclear what vcf he used & how it was filtered


#### add phenotypes to .fam file
sbatch /rhome/jmarz001/bigdata/Ag-Competition/scripts/gwas/2_phenotypes.R


# make pca covar file
plink --vcf imputed_filter.recode.vcf --double-id --allow-no-sex --allow-extra-chr --pca 10 --out all_traits_pca
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
ls ASSOC_*.assoc.txt > all_assoc_files.txt
sbatch /rhome/jmarz001/bigdata/Ag-Competition/scripts/gwas/5b_clump_sig_regions.sh
sbatch --array=1-$ARRAY_LIM /rhome/jmarz001/bigdata/Ag-Competition/scripts/gwas/5b_clump_indv_regions.sh

# plot results
sbatch /rhome/jmarz001/bigdata/Ag-Competition/scripts/gwas/5c_clumped_indv_sites_manhattan.R





## mkdir CCII_greenhouse_exp_gwas
cp /rhome/jlandis/bigdata/RADSeq/GEMMA/output/MixedModel_*assoc.txt CCII_greenhouse_exp_gwas/


cp rhome/dkoenig/bigdata/BARLEY_CCII_PARENTS_RAD/DATA/OUTPUT/VARCALLS/FILTERED.vcf.gz rhome/jmarz001/bigdata/Ag-Competition/results/gwas/PARENTS.vcf.gz

TMPVCF=PARENTS.vcf.gz
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' $TMPVCF > pos_data2.txt

chr1H   145277
chr1H   145277

#final parental rad
chr1H   145289
chr1H   159194

/rhome/jmarz001/shared/MANUSCRIPT_FILES/2019_LANDIS/CCIISELECT/FINAL_PARENTAL_RAD.vcf.gz
#FINAL_PROGENY_RAD.vcf: The unimputed progeny calls

RAD_CALLS_FILT.vcf
~/bigdata/BARLEY_CCII_PARENTS_RAD/DATA/INPUT/RECALL_FILT_NAMES.bed





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
