This readme concerns the files produced for a gwas using the software gemma

File directory
    /rhome/jmarz001/bigdata/Ag-Competition/results/gwas

All of the GWAS result files have the following structure
    ASSOC_\\d+{_lmm}.assoc.<file suffix>
More specifically, all files start with "ASSOC_", followed by a digit (starting at 6), some files have the string "_lmm" indicating the model used for that gwas, ".assoc.", and a file suffix

so all GWAS result files can be found by running this command   
    ls ASSOC*

Trait/file numbers
    the file 'trait_name_to_col_numbers.tsv' connects the tested trait name with its corresponding trait number and file names.
    All gwas result numbers start with 6 (referencing the format of the phenotype (.fam) file, where 6 is the number of the first testable phenotype in the file)
    The trait - number connection is reiterated here for convenience
    FT 6 (Flowering Time)
    TOTAL_MASS 7
    SEED_WEIGHT_100 8
    MASS_PER_PLANT 9
    Germination 10
    FECUNDITY 11

File suffixes
*.log
    contains details about each sepecific gemma run
*.nosex
    list of line IDs without sex information (all lines in this experiment)
*.txt
    raw association text file produced by gemma
*.tmp
    filtered & formatted version of *.txt file, retains only chromosome, snp id, base pair, and p-value; temporary imput file for producing clumped gwas results
*.clumped
    result of running the clump command with plink software; file contains snps identified as the top-snp identified as part of a gwas peak

"_lmm"
all traits were run with lm and lmm models to compare results, so lmm versions should exist for all gwas files

The previous holds true for all gwas analysis files also contained in the sub-directory 'CCII_greenhouse_exp_gwas'. The file with trait names, numbers, and file names is 'CCII_GH_trait_file_nums.tsv'.