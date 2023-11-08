library(tidyverse)
library(readr)
library(dplyr)
library(ggpubr)
library(ggplot2)
library(googlesheets4)
library(tidyr)

df <- read_sheet('https://docs.google.com/spreadsheets/d/1mxSidDcodD7-Iju9JJZ_YoejhjEqt6JDad3LMGOh61s/edit#gid=1749035245')
df1 <- read_sheet('https://docs.google.com/spreadsheets/d/1Rb1oN4yeqcQDFKfezf3uGLEOlp-s8x1iAIpAO6cqT4E/edit#gid=682364943')
df2 <- read_sheet('https://docs.google.com/spreadsheets/d/1pKOlthCtyF-T8bbU_96xfoUDcQbbVSnP3jAk6iP3egc/edit#gid=326137907')

###Filtering Necessary columns
df <- select(df, 1:4)
df1 <- select(df1, 1:2)
dfs <- full_join(df, df1, by = ("PLOT_ID"))
dfs$PLOT_ID <- as.numeric(dfs$PLOT_ID)
dfs <- subset(dfs, dfs$PLOT_ID <= 1036, )

remove(df1)
remove(df)

df <- full_join(dfs, df2, by = ("PLOT_ID"))
df <- select(df, !c(Date, ROW, `albino count (not included in germination / survival since they don't survive)`, PLANT_ID, Plot_Survival))
df$Generation <- gsub("^1_.*", 18, df$Genotypes)
df$Generation <- gsub("^2_.*", 28, df$Generation)
df$Generation <- gsub("^3_.*", 50, df$Generation)
df$Generation <- gsub("^7_.*", 58, df$Generation)
df$Generation <- gsub("^*.*_.*", 0, df$Generation)
df$Fecundity <- df$`Brown Bag Weight`/(df$`100 seed weight`/100)
df$Fitness <- df$Fecundity * df$Plot_Germination

remove(dfs)
remove(df2)

#### Check for missing data

col1 <- as.vector(df$PLOT_ID)
col1 <- as.integer(col1)
seq2 <- (9:1036)
setdiff(seq2, df$PLOT_ID)


