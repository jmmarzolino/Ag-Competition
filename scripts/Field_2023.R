library(tidyverse)
library(readr)
library(dplyr)
library(ggpubr)
library(ggplot2)

df <- read_delim("~/Desktop/FT_DAYS_2022-2023.xlsx - FIELD_LAYOUT.tsv",
                 delim = "\t", escape_double = F)
df1 <- read_delim("~/Desktop/Seed Weights - Field 2023.tsv",
                  delim = "\t", escape_double = F)

dfs <- full_join(df, df1, by = ("PLOT_ID"))



