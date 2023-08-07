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

### Join data frames by PLOT_ID and clean environment
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


#### Quality Control
df %>% 

### Create new data tables for Two conditions
Mixed_df <- subset(df, df$Condition == "mixed")
Single_df <- subset(df, df$Condition == "single")

#### Make a Table for Average Atlas Seed Weight
agg_tbl <- subset(df, df$Genotypes == "48_5")
aggMixed <- subset(agg_tbl, agg_tbl$Condition == "mixed")
aggMixed$Avg <- summarise(aggMixed, AvgAtlas = mean(`Brown Bag Weight`))
fn <- function(x) x/2
aggMixed$AvgNew <- apply(aggMixed[c("Avg")], 
                         FUN = fn,
                         MARGIN = 2)

### 
Mixed_new <- aggregate(`Brown Bag Weight` ~ Genotypes, Mixed_df, mean)
Mixed_new$Add <- apply(Mixed_new[c("Brown Bag Weight")],
                                FUN = fn,
                                MARGIN = 2)

Mixed_new <- full_join(Mixed_new, aggMixed, by = "Genotypes")



Mixed_df$Average <- 
  

Mixed_df$`Brown Bag Weight`

Mixed_df <- full_join(Mixed_df, aggMixed, by = "PLOT_ID")

##### Flowering Times over Generations
df$Generation <- as.numeric(df$Generation)
ggplot(df, aes(Generation, FT_DAYS)) +
  geom_jitter(alpha = .3) +
  geom_smooth(method = lm) +
  labs(y = "Flowering Time (Days)",
       title = "Flowering Time over Generations")

#### Mixed vs. Single Yield over Generations
ggplot(df, aes(Generation, `Brown Bag Weight`, color = Condition)) +
  geom_jitter(alpha = .3) +
  geom_smooth(method = lm) +
  labs(y = "Total Mass",
       title = " Mixed vs. Single Yield over Generations")


####
df <- full_join(df, agg_tbl, by = "PLOT_ID")


"48_5" %in% df$Genotypes

df$ExpectedWPS <- 

