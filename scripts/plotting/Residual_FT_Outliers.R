#!/usr/bin/env Rscript

#SBATCH --ntasks=1
#SBATCH --mem=30G
#SBATCH --time=02:00:00
#SBATCH --output=/rhome/jmarz001/bigdata/Ag-Competition/competition1.stdout
#SBATCH -p koeniglab

library(tidyverse)
library(readr)
library(dplyr)
library(ggpubr)
library(ggplot2)
library(googlesheets4)
library(tidyr)
library(car)
library(gridExtra)
library(dunn.test)


FT_2022_ranking_mixed$Genotypes <- as.factor(FT_2022_ranking_mixed$Genotypes)
FT_2022_ranking_mixed <- FT_2022_ranking_mixed %>% filter(Genotypes != "51_5_3" & Genotypes != '2_17' & Genotypes != '2_202' & Genotypes != '53_5' & Genotypes != '2_148')

a <- ggplot(FT_2022_ranking_mixed, aes(x= Genotypes, y= Flowering_Date)) +
  geom_point() +
  scale_y_continuous(breaks = seq(85, 135, 5)) +
  labs(y = "Flowering Time",
      title = 'Residual Outliers Mixed 2022')

FT_2022_ranking_single$Genotypes <- as.factor(FT_2022_ranking_single$Genotypes)
FT_2022_ranking_single <- FT_2022_ranking_single %>% filter(Genotypes != "1_61" & Genotypes != '2_74' & Genotypes != '3_4' & Genotypes != "1_165_3" & Genotypes != "1_191" & Genotypes != '7_11')

b <- ggplot(FT_2022_ranking_single, aes(Genotypes, Flowering_Date)) +
  geom_point() +
  scale_y_continuous(breaks = seq(65, 125, 5)) +
  labs(y = "Flowering Time",
       title = "Residual Outliers Single 2022")

FT_2023_ranking_mixed$Genotypes <- as.factor(FT_2023_ranking_mixed$Genotypes)
FT_2023_ranking_mixed <- FT_2023_ranking_mixed %>% filter(Genotypes != '2_142')

c <- ggplot(FT_2023_ranking_mixed, aes(Genotypes, FT_DAYS)) +
  geom_point() +
  scale_y_continuous(breaks = seq(100, 135, 5)) +
  labs(y = "Flowering Time",
       title = "Residual Outliers Mixed 2023")

FT_2023_ranking_single$Genotypes <- as.factor(FT_2023_ranking_single$Genotypes)
FT_2023_ranking_single <- FT_2023_ranking_single %>% filter(Genotypes != '7_211' & Genotypes != '7_74' & Genotypes != '72_7_1')

d <- ggplot(FT_2023_ranking_single, aes(Genotypes, FT_DAYS)) +
  geom_point() +
  scale_y_continuous(breaks = seq(110, 140, 5)) +
  labs(y = "Flowering Time",
       title = "Residual Outliers Single 2023")

y <- arrangeGrob(a, b, c, d, top = "Residual Flowering Time Outliers", nrow = 2, ncol =2)
ggsave("scripts/plotting/Residual_FT_Outliers.png", y, width = 14, height = 10)
