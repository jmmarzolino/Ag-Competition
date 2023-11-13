library(dplyr)
library(googlesheets4)


ID_and_Geno <- read_delim("~/Documents/GitHub/Ag-Competition/Fig1.tsv")
Missing_ID <- slice(df, 999:1029)
Only_Missing_ID <- select(Missing_ID, "PLOT_ID", "Genotypes")
                          