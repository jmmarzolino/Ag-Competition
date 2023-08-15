library(tidyverse)
library(readr)
library(dplyr)
library(ggpubr)
library(ggplot2)
library(tidyr)

df <- read_delim("~/Documents/GitHub/Ag-Competition/Fig1.tsv")
Combinedfull <- read_delim("~/Documents/GitHub/Ag-Competition/Combinedfull.tsv")

### Create new data tables for Two conditions
Mixed_df <- subset(df, df$Condition == "mixed")
Single_df <- subset(df, df$Condition == "single")

#### Restructure Single_df to average the replicates

Single_df1 <- subset(Single_df, replicate == "rep 1")
Single_df2 <- subset(Single_df, replicate == "rep 2")
colnames(Single_df2)[1:14] <- paste(colnames(Single_df2)[c(1:14)], '_2', sep = '_')
Single_df <- inner_join(Single_df1, Single_df2, by = c("Genotypes" = "Genotypes__2","Condition" = "Condition__2"))
Single_df$AvgFec <- (Single_df$Fecundity + Single_df$Fecundity__2)/2
Single_df$AvgFit <- (Single_df$Fitness + Single_df$Fitness__2)/2

#### Make a Table for the Average Atlas Seed Weight (Atlas x Atlas = 48_5)
agg_tbl <- subset(df, df$Genotypes == "48_5")
agg_tbl <- filter(agg_tbl, agg_tbl$Condition != "mixed")
agg_tbl$AvgAtlas <- (sum(agg_tbl$Fecundity)/2)/2

#### Restructure Mixed_df to average replicates
Mixed_1 <- subset(Mixed_df, replicate == "rep 1")
Mixed_2 <- subset(Mixed_df, replicate == "rep 2")
colnames(Mixed_2)[1:14] <- paste(colnames(Mixed_2)[c(1:14)], '_2', sep = '_')
Mixed_df <- inner_join(Mixed_1, Mixed_2, by = c("Genotypes" = "Genotypes__2","Condition" = "Condition__2"))

#### Making the Expected Fecundity column (number of seeds) (1/2(Avg of Atlas x Atlas line) + 1/2(Avg of Geno x Geno line))
Mixed_df$AvgGeno <- ((Mixed_df$Fecundity + Mixed_df$Fecundity__2)/2)/2
Mixed_df$ExpFecundity <- Mixed_df$AvgGeno + (sum(agg_tbl$Fecundity)/2)/2

#### Add an Expected Contribution column to Single and Mixed dataframes
Single_df$ExpContribution_s <- Single_df$`Brown Bag Weight`/10
Mixed_df$ExpContribution_m <- (Mixed_df$`Brown Bag Weight`/5)
Mixed_df$AvgFec <- (Mixed_df$Fecundity + Mixed_df$Fecundity__2)/2
Mixed_df$AvgFit <- (Mixed_df$Fitness + Mixed_df$Fitness__2)/2

#### 02a_Single_Fitness_Over_Generations.R 
#### (fitness is defined as Plot Germination * Fecundity)
Single_df$Generation <- as.numeric(Single_df$Generation)
ggplot(Single_df, aes(Generation, AvgFit)) +
  geom_jitter(alpha =.3) + 
  geom_smooth(method = lm) +
  labs(x = "Generation",
       y = "Fitness",
       title = "Fitness over Generations") 
ggsave("Sin_Fitness_Over_Generations.png")

##### 2ai_FT_over_Generations.R
ggplot(Combinedfull, aes(Generation, FT_DAYS)) +
  geom_jitter(alpha = .3) +
  geom_smooth(method = lm) +
  labs(y = "Flowering Time (Days)",
       title = "Flowering Time over Generations") +
  facet_wrap(~Condition)

ggsave("FT_over_Generations.png")

### 02b_Single_Fitness_over_Generation.R

ggplot(Single_df, aes(Generation, AvgFit)) +
geom_jitter(alpha = .5) +
  geom_smooth(method = lm) +
  labs(x = "Generation",
       y = "Average Fitness",
       title = "Average Fitness Over Generation") 
ggsave("Sin_Fitness_over_Generation.png")

### 02c_Fecundity_Distributions
ggplot(Combinedfull, aes(x = AvgFec)) +
  geom_histogram() +
  labs(x = "Fecundity",
       y = "Frequency",
       title = "Fecundity Over Generations") +
  facet_wrap(~Generation)
ggsave("Fecundity_Over_Gen_Distributions.png")

ggplot(Combinedfull, aes(x = AvgFec, group = Generation, color = Generation)) +
  geom_histogram(aes(fill = Generation), alpha = .5) +
  scale_fill_brewer(palette = "Blues")

### 02ci_100_SW_Distributions
ggplot(Combinedfull, aes(x = AvgSW_100)) +
  geom_histogram() +
  labs(x = "Fitness",
       y = "Frequency",
       title = "Fitness Over Generations") +
  facet_wrap(~Generation)
ggsave("100SW_Over_Gen_Distributions.png")



ggplot(Combinedfull, aes(x = `100 seed weight`, group = Generation, color = Generation)) +
  geom_histogram(aes(fill = Generation), alpha =.5) +
  scale_fill_brewer(palette = "Blues")

### 02di_Fitness_Compared_to_Atlas.R

ggplot(Single_df, aes(Generation, AvgFit)) +
  geom_jitter() +
  ylim(0, 30000) +
  geom_hline(yintercept = 30000, color = "red") +
  labs(y = "Relative Fitness",
       title = "Relative Fitness Over Generations") +
  facet_wrap(~Generation)

ggplot(Single_df, aes(Generation, AvgFit)) +
  geom_jitter() +
  ylim(0, 30000) +
  geom_hline(yintercept = 30000, color = "red") +
  geom_smooth(method = lm)
  labs(y = "Relative Fitness",
       title = "Relative Fitness Over Generations")
  




##### Write Delims

write_delim(Mixed_df, "Mixed_df.tsv")
write_delim(Single_df, "Single_df.tsv")
