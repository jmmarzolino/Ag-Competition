#!/bin/bash
#SBATCH --job-name=gwas
#SBATCH --output=/rhome/jmarz001/bigdata/Ag-Competition/scripts/gwas/3_multivariate_association.stdout
#SBATCH --ntasks=4
#SBATCH --mem=100G
#SBATCH -p koeniglab

# GEMMA 0.98
cd /rhome/jmarz001/bigdata/Ag-Competition/results/gwas

# define variables
GENO=all_traits
KINSHIP=../output/related_matrix.cXX.txt
PCS=pca.txt
PHENO_NAME="ASSOC_MULI"

### Univariate Linear Mixed Model
# -km 1: PLINK binary ped as relatedness matrix; 1 is default
# -k [file] specifies the relatedness matrix file name
# if no covar use lm instead of lmm
# -c [file] specify input covariates file name (optional); an intercept term is needed in the covariates file

#/rhome/jmarz001/software/gemma0.98.5 -bfile "$GENO" -k "$KINSHIP" -c "$PCS" -lmm 4 -miss 1 -notsnp -outdir $OUT -n $SLURM_ARRAY_TASK_ID -o $PHENO_NAME -debug


### Multivariate Linear Mixed model
#-n [num1] [num2] [num3]
#-n option is employed to specify which phenotypes in the phenotype file are used for association tests

# restrict the number of phenotypes to be < ten
#when a small proportion of phenotypes are partially missing, one can impute these missing values before association tests w/ -predict

/rhome/jmarz001/software/gemma0.98.5 -bfile "$GENO" -k "$KINSHIP" -lmm 4 -n 1 2 3 4 5 6 -outdir . -o $PHENO_NAME -debug

#-c "$PCS"