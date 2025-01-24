#!/usr/bin/env Rscript
#SBATCH --mem=30G
#SBATCH --time=02:00:00
#SBATCH --output=/rhome/jmarz001/bigdata/Ag-Competition/scripts/004b_compare_w_CCII_GH.stdout
#SBATCH -p koeniglab

library(tidyverse)
library(ggpubr)
library(data.table)

setwd("/bigdata/koeniglab/jmarz001/Ag-Competition/data")
source("../scripts/CUSTOM_FNS.R")

df <- fread("")
