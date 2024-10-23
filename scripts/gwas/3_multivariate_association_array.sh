#!/bin/bash
#SBATCH --job-name=gwas
#SBATCH --output=/rhome/jmarz001/bigdata/IPK_Analysis/scripts/GWAS/3_multivariate_association_array.stdout
#SBATCH --ntasks=2
#SBATCH --mem=80gb
#SBATCH -p koeniglab

# GEMMA 0.98
cd /rhome/jmarz001/bigdata/IPK_Analysis/results/GWAS

# define variables
GENO=all_traits
KINSHIP=output/related_matrix.cXX.txt
PCS=pca.txt
OUT=derived_traits_ASSOC
COL=$(expr $SLURM_ARRAY_TASK_ID + 5)

PHENO_NAME="ASSOC_${COL}"
#TOTAL_COLS=$(head -n1 all_traits.fam | awk '{print NF}')
#echo $TOTAL_COLS

### Univariate Linear Mixed Model
# -km 1: PLINK binary ped as relatedness matrix; 1 is default
# -k [filename] specifies the relatedness matrix file name
#/rhome/jmarz001/software/gemma0.98.5 -bfile "$GENO" -k "$KINSHIP" -c "$PCS" -lmm 4 -miss 1 -notsnp -outdir $OUT -n $SLURM_ARRAY_TASK_ID -o $PHENO_NAME -debug

/rhome/jmarz001/software/gemma0.98.5 -bfile "$GENO" -k "$KINSHIP" -c "$PCS" -lmm 4 -outdir $OUT -n $SLURM_ARRAY_TASK_ID -o $PHENO_NAME -debug
## if no covar use lm instead of lmm
#/rhome/jmarz001/software/gemma0.98.5 -bfile FT -k $RMATRIX -notsnp -lmm 4 -miss 1 -outdir "ASSOC" -n 1 -o ASSOC_FT -debug
#One can specify a different column as the phenotype column by using “-n [num]”, where ”-n 1” uses the original sixth column as phenotypes, and “-n 2” uses the seventh column, and so on and so forth.

# -c [filename] specify input covariates file name (optional); an intercept term is needed in the covariates file

### Multivariate Linear Mixed model
#-n [num1] [num2] [num3]
#-n option is employed to specify which phenotypes in the phenotype file are used for association tests.
#it is highly recommended to restrict the number of phenotypes to be small, say, less than ten.
#when a small proportion of phenotypes are partially missing, one can impute these missing values before association tests
#-predict
