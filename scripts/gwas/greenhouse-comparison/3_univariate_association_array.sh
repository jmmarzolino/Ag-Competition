#!/bin/bash
#SBATCH --job-name=gwas
#SBATCH --output=/rhome/jmarz001/bigdata/Ag-Competition/scripts/gwas/greenhouse-comparison/3_univariate_association_array.stdout
#SBATCH --ntasks=2
#SBATCH --mem=80gb
#SBATCH --time=2:00:00
#SBATCH -p koeniglab

# GEMMA 0.98
cd /rhome/jmarz001/bigdata/Ag-Competition/results/gwas/CCII_greenhouse_exp_gwas

# define variables
GENO=gh
KINSHIP=output/related_matrix.cXX.txt
PCS=pca.txt
COL=$(expr $SLURM_ARRAY_TASK_ID + 5)
PHENO_NAME="ASSOC_${COL}"

### Univariate Linear Mixed Model
# -km 1: PLINK binary ped as relatedness matrix; 1 is default
# -k [filename] specifies the relatedness matrix file name
#/rhome/jmarz001/software/gemma0.98.5 -bfile "$GENO" -k "$KINSHIP" -c "$PCS" -lmm 4 -miss 1 -notsnp -outdir $OUT -n $SLURM_ARRAY_TASK_ID -o $PHENO_NAME -debug

## if no covar use lm instead of lmm
#One can specify a different column as the phenotype column by using “-n [num]”, where ”-n 1” uses the original sixth column as phenotypes, and “-n 2” uses the seventh column, and so on and so forth.

# univariate linear MIXED model
/rhome/jmarz001/software/gemma0.98.5 -bfile "$GENO" -k "$KINSHIP" -c "$PCS" -lmm 4 -n $SLURM_ARRAY_TASK_ID -outdir . -o ${PHENO_NAME}_lmm -debug


## adding another quick gwas with only greenhouse results
PHENO_NAME_GHONLY="ASSOC_GH_${COL}_lmm"

/rhome/jmarz001/software/gemma0.98.5 -bfile gh_only -k output/related_matrix_gh_only.cXX.txt -c gh_only_pca.txt -lmm 4 -n $SLURM_ARRAY_TASK_ID -outdir . -o $PHENO_NAME -debug
