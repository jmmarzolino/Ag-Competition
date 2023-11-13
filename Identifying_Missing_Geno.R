library(tidyverse)
library(readr)
library(dplyr)
library(ggpubr)
library(ggplot2)
library(tidyr)
library(googlesheets4)


ID_and_Geno <- read_delim("~/Documents/GitHub/Ag-Competition/Fig1.tsv")
Missing_data <- read_sheet("https://docs.google.com/spreadsheets/d/1NJKf-l3UnXrBJIoOADxhlGTuIfnYzr8yZpKUpqbrSJ4/edit#gid=2141534954")
Missing_data$Genotypes <- as.character(Missing_data$Genotypes)
Missing_ID <- slice(df, 999:1029)
Only_Missing_ID <- select(Missing_ID, "PLOT_ID", "Genotypes", "Condition", "replicate", "Bed_2022", "Row_2022")
Full_Data <- inner_join(Only_Missing_ID, Missing_data, by = c("Genotypes", "Condition", "replicate"))


                          