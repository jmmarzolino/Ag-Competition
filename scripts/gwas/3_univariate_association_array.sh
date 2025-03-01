#!/bin/bash
#SBATCH --job-name=gwas
#SBATCH --output=/rhome/jmarz001/bigdata/Ag-Competition/scripts/gwas/3_univariate_association_array.stdout
#SBATCH --ntasks=2
#SBATCH --mem=80gb
#SBATCH -p koeniglab

# GEMMA 0.98
cd /rhome/jmarz001/bigdata/Ag-Competition/results/gwas

# define variables
GENO=all_traits
KINSHIP=../output/related_matrix.cXX.txt
PCS=pca.txt
COL=$(expr $SLURM_ARRAY_TASK_ID + 5)
PHENO_NAME="ASSOC_${COL}"

### Univariate Linear Mixed Model
# -km 1: PLINK binary ped as relatedness matrix; 1 is default
# -k [filename] specifies the relatedness matrix file name
#/rhome/jmarz001/software/gemma0.98.5 -bfile "$GENO" -k "$KINSHIP" -c "$PCS" -lmm 4 -miss 1 -notsnp -outdir $OUT -n $SLURM_ARRAY_TASK_ID -o $PHENO_NAME -debug

## if no covar use lm instead of lmm
#One can specify a different column as the phenotype column by using “-n [num]”, where ”-n 1” uses the original sixth column as phenotypes, and “-n 2” uses the seventh column, and so on and so forth.

# test univariate linear model & univariate linear MIXED model
/rhome/jmarz001/software/gemma0.98.5 -bfile "$GENO" -c "$PCS" -lm 4 -n $SLURM_ARRAY_TASK_ID -outdir . -o $PHENO_NAME -debug

#/rhome/jmarz001/software/gemma0.98.5 -bfile "$GENO" -k "$KINSHIP" -c "$PCS" -lmm 4 -n $SLURM_ARRAY_TASK_ID -outdir . -o ${PHENO_NAME}_lmm -debug
