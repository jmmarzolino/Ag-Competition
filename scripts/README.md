rhome/jmarz001/bigdata/Ag-Competition/scripts/
---

Data_Manipulation

Manipulating Data

Dist...Quality Check






# 000_row_and_col_boxplots.R
Experimental replicates were not correctly randomized, so this script attempts to check for a pattern of differences between field 'replicates' based on planting location. It plots boxplots for both rows (field short rows) and columns (field beds) per experimental year (2022 or 2023) for measured phenotypes (100 seed weight, Flowering time).

#000_single_v_mixed_PerSeedWeight_correlation.R
Measures correlation of plot replicates for 2 phenotypes for both mixed- and single-genotype conditions.

#001_compare_ft_over_years.R
Compares flowering time between experimental years with anova tests for significant difference in average flowering time per year.
Tests explanatory variables besides year as a source of variance in flowering time.
Checks correlation coefficient across experimental conditions per year.
Produces boxplots and density plots comparing flowering time over 1) Experimental Condition  2) CCII generations  3) CCII generation and experimental condition  4) experimental replicates across years.

#001_format_2022.R & #001_format_2023.R
Load and filter raw data of genotype, field position, field-phenotyped flowering time. Calculate average flowering time for generation, experimental condition, and replicate, separately. Calculates average flowering time for for parental genotypes in field year 2023. Explores variables as explantory sources of variance in flowering time.


############################################
#Ag_Comp_Graphs_2.R
Needs to be adapted for data on biocluster instead of local computer.
Creates plots comparing seed weight phenotypes (fitness) between CC II progeny genotypes and Atlas parent.

#001_format_plot_floweringtime.R
Is this script redundant with others?
#competition1.R
Probably redundant with other scripts.

#Field 2023.R
Script acts as a quality control for the data used in the Ag_Comp Project. Contains correlation graphs for replicates in relation to Total seed weight, Fecundity, 100 seed weight single, and 100 seed weight mixed. Also shows histograms to compare the shift in averages over generations for Total seed weight and 100 seed weight.
