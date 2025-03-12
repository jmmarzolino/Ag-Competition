# rhome/jmarz001/bigdata/Ag-Competition/scripts/

---

### CUSTOM_FNS.R

Defines custom functions written for and used throughout scripts in this project.  

- "add_generation" adds a generation column to a data frame based on the UCRKL number in the genotype column
- "tidy_text_substitution" finds common abbreviations or formats used in coding environment & replaces them with a human-readable text string, for cleaner plotting.


### 001_format_raw_data.R

Load and filter raw data of genotype, field position, field-phenotyped flowering time.

### 002_replicate_correlations.R

Checks for replicability of experiment replicates.  

- Calculate pearson's correlation value between replicates 1 & 2 within the two experimental years.  
- Generates dot plots of trait values for rep 1 vs. rep 2 & plots residuals for the model rep 2 ~ rep 1.  


### 002_year_correlations.R

Calculate and plot trait values for experiment years (2022 vs 2023).  

### 002b_derived_traits_and_QC.R

Filter outliers from raw data based on threshholds of (median +/- 2*IQR). Also quantifies amount of filtered data.  

- Plots trait distributions as histograms with threshhold lines, for data before and after outlier filtering.  

Uses collected phenotypes (total mass, number of plants, 100-seed weight) to calculate derived phenotypes (mass per plant, seed count, fecundity).  

**TMP: identifies & filters replicates with large differences??**

### 002c_derived_data_replicate_correlations.R

Calculates & plots correlation between experiment years & within-year replicates for the new derived traits (ie. including mass per plant, seed count, fecundity).  

### 003a_blups_heritability.R

Imports filtered derived trait data & fits a mixed model for each trait. Extracts BLUPs for each genotype to remove the effects of experimental year and random factors from the genetically attributible trait value. Also estimates trait's heritability value based on model fit. Combines data for genotype blups & trait heritability values, respectively, and writes out as data tables. 

### 003b_change_over_gens.R

Calculate each traits average and variance per generation, for filtered, derived & blup traits.  
Test for normal distributions within each generation with Shapiro-Wilkes test & test for equal variance between generations with Levene test.  
Then test for significant differences in generations means with ANOVA or Kruskal-Wallis test. Follow up significant test results with posthoc Dunn test to investiage if significant change in mean happened between parent and F18 generations or F18 and F58 generations.  

### 003b_response_selection.R

Calculate response (total trait change between two generations average, divided by the number of generations between the two periods; aka average amount trait changes over generations) & selection estimate (response / heritability) for all traits based off of genotype BLUPs.  

- Calculates response in normal/default units & standard deviations.
- Response calculated for two time periods: between parents and F18 generation, and between F18 and F58 generations
- Plots response between generations, scaled for comparison between traits, & selection for each trait

### 004_comparing_haplotypes.R

Purpose is to re-sample phenotype data to match its representation in source-generation population. Does so by sampling greater or fewer numbers of phenotypes attributible to a given genetic background, based on how common that genetic background was at a given time point.  
To do so, it uses a file connecting sequenced genotypes with their corresponding haplotype number. This file is based on individual sequencing & pooled-population sequencing, based on random sampling of seed for each tested generations.  

- Sums the frequency of each haplotype number in a given generation & divides by the total number of samples to estimate the fraction of the population's individuals with the same genetic background (haplotype). 
- Averages trait values across haplotype to represent the genetic background with a common value. 
- Samples haplotype-phenotypes based on how common a given haplotype was in each generation of the experiment (18, 28, 50, and 58). 
- Calculates population's weighted or 'haplotype-population adjusted' trait averages. 

### 004b_change_over_gens_HapRepPop.R

Uses the same basis as '004_comparing_haplotypes.R' to sample phenotypes representative of a given generation (haplotype-represented populations).  

- plots histogram of trait frequencies after re-sampling adjustment
- finds average and variance for each generation
- repeats tests in '003b_change_over_gens.R' script for normal trait distributions, equal variance between generations, significant changes in mean between generations, and the posthoc results for specific inter-generational average difference



## GWAS/

### 1_gwas.sh

Master GWAS script, controls & submits all other scripts defined in folder.  

- calls & filters vcf file
- creates plink files for analysis
- submits scripts in order, with appropriate array numbers

### 2_phenotypes.R

### 3_multivariate_association.sh

### 3_univariate_association_array.sh

### 4_manhattan_plot.R

### 5a_format_chr_pos.R

### 5b_clump_indv_regions.sh

### 5b_clump_sig_regions.sh

### 5c_clumped_indv_sites_manhattan.R

### 5c_clumped_sites_manhattan.R



# plotting/

### 2_blup_correlations.R

### 2_raw_value_vs_blups_correlation.R

### 2_trait_distributions_GenCmnValLine.R

### 2_traits_distributions.R

### 3_hapadjusted_genavg_over_gens.R

### 3_hapadjusted_traits_over_gens_scatter_and_line_and_box.R

### 3_trait_avg_and_var_over_gens.R

### 3_traits_over_gens_scatter_and_line.R


---
---
---


# TMP


### 000_row_and_col_boxplots.R

Experimental replicates were not correctly randomized, so this script attempts to check for a pattern of differences between field "replicates" based on planting location. It plots boxplots for both rows (field short rows) and columns (field beds) per experimental year (2022 or 2023) for measured phenotypes (100 seed weight, Flowering time).

### 000_single_v_mixed_PerSeedWeight_correlation.R

Measures correlation of plot replicates for 2 phenotypes for both mixed- and single-genotype conditions.

### 001_compare_ft_over_years.R

Compares flowering time between experimental years with anova tests for significant difference in average flowering time per year.
Tests explanatory variables besides year as a source of variance in flowering time.
Checks correlation coefficient across experimental conditions per year.
Produces boxplots and density plots comparing flowering time over 1) Experimental Condition  2) CCII generations  3) CCII generation and experimental condition  4) experimental replicates across years.

### Ag_Comp_Graphs_2.R

Needs to be adapted for data on biocluster instead of local computer.
Creates plots comparing seed weight phenotypes (fitness) between CC II progeny genotypes and Atlas parent.

### Field 2023.R

Script acts as a quality control for the data used in the Ag_Comp Project. Contains correlation graphs for replicates in relation to Total seed weight, Fecundity, 100 seed weight single, and 100 seed weight mixed. Also shows histograms to compare the shift in averages over generations for Total seed weight and 100 seed weight.
