#!/usr/bin/env Rscript
#SBATCH --mem=30G
#SBATCH --time=02:00:00
#SBATCH --output=/rhome/jmarz001/bigdata/Ag-Competition/scripts/002_QC.stdout
#SBATCH -p koeniglab

library(tidyverse)
library(ggpubr)
library(data.table)
library(car)
library(gridExtra)
library(dunn.test)

setwd("/bigdata/koeniglab/jmarz001/Ag-Competition")
df <- fread("data/JOINED_PHENOTYPES.tsv")


# filter any problematic outliers that remain
df_long <- df %>% select(-c(BED, ROW)) %>% pivot_longer(cols=c(FT, TOTAL_MASS, SEED_WEIGHT_100), names_to='PHENOTYPE', values_to="VALUE")


# Plot trait distributions
# Distribution of number of plants survived

df2 <- df %>% filter(Plants != 0 & !is.na(Plants))

df2 %>%
  ggplot() +
    geom_histogram(aes(Plants), fill = "#eca50b") +
    #scale_y_continuous(breaks = seq(0, 1000, 50)) +
    #scale_x_continuous(breaks = seq(0, 12, 1)) +
    geom_vline(aes(xintercept = mean(df2$Plants)), color = "#0c820c") +
    geom_vline(aes(xintercept = median(df2$Plants)), color = "#3d3de8") +
    labs(color = "Mean and Median",
        x = "number of plants",
       title = "Distribution of number of plants")
ggsave("results/plot_survival.png")#, width = 12, height = 10)




## distribution of experimental conditions
png("distribution_condition.png")
ggplot(df3, aes(x=FT, group=Condition, color=Condition, fill=Condition)) +
geom_density(alpha=0.5) +
theme_minimal()
dev.off()


## distribution of Generations
png("distribution_generation.png")
ggplot(df3, aes(x=FT, group=as.factor(Generation), color=as.factor(Generation), fill=as.factor(Generation))) +
geom_density(alpha=0.5) +
theme_minimal()
dev.off()

## distribution of Generations - single (non-mixed) plots only
png("distribution_generation_single_condition.png")
df3 %>% filter(Condition == "single") %>% ggplot(aes(x=FT, group=as.factor(Generation), color=as.factor(Generation), fill=as.factor(Generation))) +
geom_density(alpha=0.5) +
theme_minimal()
dev.off()



## distribution of Replicates
png("distribution_Replicates.png")
ggplot(df3, aes(x=FT, group=Replicate, color=Replicate, fill=Replicate)) +
geom_density(alpha=0.5) +
theme_minimal()
dev.off()











#p1 <- graph_correlation(cmp, "Total Weight (g)", "2021-2022") +
#  stat_cor(label.y = 170)
#p2 <- graph_correlation(cmp, "Relative Fecundity", "2021-2022") +
#  stat_cor(label.y = 2.5, label.x = 6.5) +
#  geom_text_repel(label = ifelse(cmp$REP_1 > 3 | cmp$REP_2 > 3,
 #                                cmp$Genotype,
  #                               ""), size = 3, hjust =1, max.overlaps = 30)



upper <- median(mix1$TOTAL_MASS, na.rm = T) + (2 * IQR(mix1$TOTAL_MASS, na.rm = T))
lower <- median(mix1$TOTAL_MASS, na.rm = T) - (2 * IQR(mix1$TOTAL_MASS, na.rm = T))
outlier_data <- mix1 %>% filter(TOTAL_MASS > upper | TOTAL_MASS < lower)








### Extrememly high 100 seed weights



ggplot(X, aes(SEED_WEIGHT_100)) + geom_histogram()
  df[which(df$SEED_WEIGHT_100 > 30),]




#df$FIT_SEED_PER_PLANT <- df$AVG_SEED_PER_PLANT/ mean(df$AVG_SEED_PER_PLANT, na.rm=TRUE)
#df$FIT_TOTAL_SEED_COUNT <- df$TOTAL_SEED_COUNT/ mean(df$TOTAL_SEED_COUNT, na.rm=TRUE)
#quantile(df$FIT_SEED_PER_PLANT, na.rm=TRUE)
#quantile(df$FIT_TOTAL_SEED_COUNT, na.rm=TRUE)
#df$POP_FIT <- df$FITNESS
#PHENO2023[which(PHENO2023$TOTAL_MASS <0),]
#range(PHENO2023$TOTAL_MASS, na.rm=TRUE)
#PHENO2023[which(PHENO2023$SEED_WEIGHT_100 <0),]
#range(PHENO2023$SEED_WEIGHT_100, na.rm=TRUE)


hist(sw$FT)

summary(sw$FT)

## identify outlier values
outlier_cutoff = quantile(test$FEC,0.75, na.rm = TRUE) + (1.5 * IQR(test$FEC, na.rm = TRUE))

# total mass
summary(sw$TOTAL_MASS)
high <- median(sw$TOTAL_MASS, na.rm=TRUE) + 2*IQR(sw$TOTAL_MASS, na.rm=TRUE)
low <- median(sw$TOTAL_MASS, na.rm=TRUE) - 2*IQR(sw$TOTAL_MASS, na.rm=TRUE)

ggplot(sw, aes(TOTAL_MASS)) + geom_histogram() + geom_vline(aes(xintercept=high)) + geom_vline(aes(xintercept=low)) + geom_vline(aes(xintercept=median(sw$TOTAL_MASS, na.rm=TRUE)), linetype="dashed") + theme_bw()


ggplot(FT_FITNESS, aes(x = FECUNDITY)) +
  geom_histogram(bins = 50) +
  geom_vline(xintercept = outlier_cutoff, color = "#F31919")








# Calculating derived phenotypes

sw$FEC <- sw$TOTAL_MASS / sw$Plants
sw$SURVIVAL <- sw$Plants / 10

# seed count based on seed weight and seed weight per 100 seeds
df$TOTAL_SEED_COUNT <- round(df$TOTAL_MASS * (100 / df$SEED_WEIGHT_100))

# estimate FECUNDITY from total and 100 seed mass
sw$FEC
sw$FITNESS <- sw$FEC * (sw$Plants/10)



# seed produced per individual
df$FEC <- df$TOTAL_SEED_COUNT/ df$Plants
df$SURVIVAL <- df$Plants / 10
df$ABS_FITNESS <- df$SURVIVAL * df$FEC

df$REL_FITNESS <- df$ABS_FITNESS / max(df$ABS_FITNESS, na.rm=TRUE)




### We have 5 genotypes where we accidentally planted 12 seeds instead of 10. For those individuals, it makes sense to adjust their survival rate relative to the 12 seeds planted

# Isolating those individuals and adjusting their survival rates

hmp <- sw %>% filter(Plants > 10)
hmp$SURVIVAL <- hmp$Plants / 12

# Adding back into original dataframe

sw <- sw %>% filter(Plants <= 10)
sw <- rbind(sw, hmp)

sw$ABS_FITNESS <- sw$SURVIVAL * sw$FEC
sw$REL_FITNESS <- sw$ABS_FITNESS / max(sw$ABS_FITNESS, na.rm=TRUE)

# Adding Centered data for the phenos of 2022

sw$FEC <- as.vector(scale(sw$FEC, center = TRUE, scale =TRUE))
sw$ABS_FITNESS <- as.vector(scale(sw$ABS_FITNESS, center = TRUE, scale = TRUE))








### Average the Replicates
sw_avg <- sw %>% group_by(Genotype, Condition, Replicate) %>% summarise(across(.cols = c(TOTAL_MASS, SEED_WEIGHT_100, FECUNDITY, FITNESS, FT), function(x) mean(x))) %>% ungroup()











## calculate measures of fitness

avg_seed_per_plant <- df %>% group_by(Generation) %>% summarise(gen_avg = mean(FIT_SEED_PER_PLANT, na.rm=TRUE))
f18 <- df %>% filter(Generation==18)

# FITNESS (single only, by year)
## population relative fitness
### each genotypes" fitness relative to whole pop in field

## generation relative fitness
### genotypes" fitness relative to avg seed of the same generation


## Atalas relative fitness

  Full_Data$FEC <- Full_Data$TOTAL_WEIGHT/(Full_Data$`100 seed weight`/100)
  Full_Data$Fitness <- Full_Data$FEC * Full_Data$Plot_Germination
  Full_Data <- na.omit(Full_Data)
  Full_Data <- Full_Data %>% filter(Genotype != "7_5")

  delete_geno <- c("396", "516", "910", "1030", "442", "444", "956", "958")
  Full_Data <- Full_Data %>% filter(!(PLOT_ID %in% delete_geno))

  #### Creates a data frame that isolates Atlas Genotype, then averages the Replicates

  Atlas_tbl <- filter(Full_Data, Full_Data$Genotype == "48_5")
  Atlas_tbl <- na.omit(Atlas_tbl)
  Atlas_tbl <- Atlas_tbl %>% mutate(Atlas_Avg_Fecundity = (sum(Atlas_tbl$FEC)/3),
                Atlas_Avg_Fitness = (sum(Atlas_tbl$Fitness)/3),
                Atlas_Avg_Total_Weight = (sum(Atlas_tbl$TOTAL_WEIGHT)/3))




  ### Importing Haplotype Data and cleaning the dataframe

  Haplo_raw <- Haplo_raw %>% select(c("Pedigree", "Haplotype", "Generation", "Family"))
  Haplo_raw$Family <- unlist(Haplo_raw$Family)
  Haplo_raw$Haplotype <- unlist(Haplo_raw$Haplotype)
  Haplo_raw <- Haplo_raw %>% mutate(Pedigree = paste0("UCRKL00000", Haplo_raw$Family))
  colnames(Haplo_raw)[which(names(Haplo_raw) == "Family")] <- "Genotype"

  Average_Haplo_rep <- full_join(Haplo_raw, Average_Haplo_rep, by = c("Genotype", "Generation"))
  Average_Haplo_rep <- Average_Haplo_rep %>% filter(TOTAL_WEIGHT != "NA")
  Average_Haplo_rep <- Average_Haplo_rep %>%  mutate(Atlas_Avg_Fec = 2363.51,
                    Atlas_Avg_Fitness = 21347.22,
                    Atlas_Avg_Total_Weight = 126.8267)
  Average_Haplo_rep$Numbers <- ifelse(Average_Haplo_rep$Condition == "mixed", 1, 0)


  ### Creating Replicate dataframe for correlation graphs

  rep1 <- Full_Data %>% filter(Replicate == "rep 1")
  rep2 <- Full_Data %>% filter(Replicate == "rep 2")
  colnames(rep2) <- paste(colnames(rep2), 2, sep = "_")
  Replicate_corr_tbl <- full_join(rep1, rep2, by = c("Condition" = "Condition_2", "Generation" = "Generation_2", "Genotype" = "Genotype_2"))
  Replicate_corr_tbl <- na.omit(Replicate_corr_tbl)

  ###### Functions to get expected fitness, fecundity, and yield per plant for both conditions
  ## Numbers col: 1 = mixed, 0 = single

  ### Fecundity

  Average_Haplo_rep$Exp_Fec_Per_Plant <- ifelse(Average_Haplo_rep$Numbers == 1,
         Exp_Fec_Mixed(Average_Haplo_rep$FEC),
         Exp_Single(Average_Haplo_rep$FEC))

  ### Fitness

  Average_Haplo_rep$Exp_Fit_Per_Plant <- ifelse(Average_Haplo_rep$Numbers == 1,
                Exp_Fit_Mixed(Average_Haplo_rep$Fitness),
                Exp_Single(Average_Haplo_rep$Fitness))

  ### Total Weight
  Exp_TW_mix <- function(x){
    TW_mix <- (x/2) + (Average_Haplo_rep$Atlas_Avg_Total_Weight/2)
    Exp_TW <- TW_mix/10
    return(Exp_TW)
  }

  Average_Haplo_rep$Exp_TW_Per_Plant <- ifelse(Average_Haplo_rep$Numbers == 1,
          Exp_TW_mix(Average_Haplo_rep$TOTAL_WEIGHT),
          Exp_Single(Average_Haplo_rep$TOTAL_WEIGHT))

  ### Adding a column for centered data
  Average_Haplo_rep <- Average_Haplo_rep %>% mutate(Centered_Fit = Fitness - mean(Fitness),
          Centered_FT = FT - mean(FT),
          Centered_Fec = Fecundity - mean(Fecundity),
          Centered_TW = TOTAL_WEIGHT - mean(TOTAL_WEIGHT, na.rm = TRUE))





  # seed count based on seed weight and seed weight per 100 seeds
  df$TOTAL_SEED_COUNT <- round(df$TOTAL_MASS * (100 / df$SEED_WEIGHT_100))


  # seed produced per individual
  df$FEC <- df$TOTAL_SEED_COUNT/ df$Plants
  df$SURVIVAL <- df$Plants / 10
  df$ABS_FITNESS <- df$SURVIVAL * df$FEC

  df$REL_FITNESS <- df$ABS_FITNESS / max(df$ABS_FITNESS, na.rm=TRUE)










  write_delim(df, "FT_FITNESS.tsv", "\t")
  write_delim(Full_Data, "Full_Data")
  write_delim(Average_Haplo_rep, "Average_Haplo_rep")












Single_2021_2022 <- PHENO_FULL_AVERAGE %>% filter(Condition == "single" & Exp_year == 2022)

### Average Total Weight Distributions
graph_tmp <- Single_2021_2022 %>% group_by(Generation, Genotype, Exp_year) %>% summarise(Avg_TW = mean(TOTAL_MASS)) %>% ungroup()
tmp <- graph_tmp %>% group_by(Generation) %>% summarise(Generation_Avg = mean(Avg_TW, na.rm = T))
graph_tmp <- full_join(graph_tmp, tmp, by = "Generation")

ggplot(graph_tmp, aes(x = Avg_TW)) +
  geom_histogram(binwidth = 4) +
  geom_vline(aes(xintercept = graph_tmp$Generation_Avg), color = "red") +
  facet_grid(~Generation) +
  scale_x_continuous(breaks = seq(0, 150, 25)) +
  labs(x = "Average Total Weight (grams)")


# Levene Test (checking for homogeneity of variance, assumption of ANOVA) - insignificant means we can assume homogeneity of variance

leveneTest(TOTAL_MASS ~ as.factor(Generation), Single_2021_2022)

### QQ-plot & Shapiro Normality for TW - (SOME GROUPS NOT NORMAL, need Kruskal Wallis and Dunn test)

h <- Single_2021_2022 %>% select(TOTAL_MASS, Generation, Genotype) %>%  arrange(Generation)
h <- h %>% mutate(TW = TOTAL_MASS)
h <- h %>% select(c("Generation", "Genotype", "TW"))
h$Genotype <- factor(h$Genotype, levels = h$Genotype[order(h$TW)])

par(mfrow = c(2, 3))
for (i in unique(h$Generation)) {
  p <- h %>% filter(Generation == i) %>% select(TW)
  qqnorm(p$TW, main = paste0("Generation ", i))
  qqline(p$TW)
}

for(i in unique(h$Generation)){
  p <- h %>% filter(Generation == i) %>% select(TW)
  s <- shapiro.test(p$TW)
  print(s)
}

### Kruskal Wallis and Dunn Tests

kruskal.test(TOTAL_MASS ~ Generation, Single_2021_2022)
dunn.test(Single_2021_2022$TOTAL_MASS, Single_2021_2022$Generation)


### Average Flowering Time Distributions

graph_tmp <- Single_2021_2022 %>% group_by(Generation, Genotype, Exp_year) %>% summarise(Avg_FT = mean(FT)) %>% ungroup()
tmp <- graph_tmp %>% group_by(Generation) %>% summarise(Generation_Avg = mean(Avg_FT, na.rm = T))
graph_tmp <- full_join(graph_tmp, tmp, by = "Generation")

ggplot(graph_tmp, aes(x = Avg_FT)) +
  geom_histogram(binwidth = 1.1) +
  geom_vline(aes(xintercept = Generation_Avg, color = "red")) +
  facet_grid(~Generation) +
  labs(x = "Average Flowering Time (Days After Sowing)")

# Levene test - unequal variances between groups

leveneTest(FT ~ as.factor(Generation), Single_2021_2022)

### QQ-plot for Single FT & Shapiro Normality | (SOME GROUPS NOT NORMAL)

h <- Single_2021_2022 %>% select(FT, Generation, Genotype) %>%  arrange(Generation)
h <- h %>% mutate(FT = FT)
h <- h %>% select(c("Generation", "Genotype", "FT"))
h$Genotype <- factor(h$Genotype, levels = h$Genotype[order(h$FT)])

par(mfrow = c(2, 3))
for (i in unique(h$Generation)) {
  p <- h %>% filter(Generation == i) %>% select(FT)
  qqnorm(p$FT, main = paste0("Generation ", i))
  qqline(p$FT)
}

for(i in unique(h$Generation)){
  p <- h %>% filter(Generation == i) %>% select(FT)
  s <- shapiro.test(p$FT)
  print(s)
}

### Kruskal Wallis and Dunn test

kruskal.test(FT ~ Generation, Single_2021_2022)
dunn.test(Single_2021_2022$FT, Single_2021_2022$Generation)


### Average Fecundity Distributions

graph_tmp <- Single_2021_2022 %>% group_by(Generation, Genotype, Exp_year) %>% summarise(Avg_Fec = mean(FECUNDITY)) %>% ungroup()
tmp <- graph_tmp %>% group_by(Generation) %>% summarise(Generation_Avg = mean(Avg_Fec, na.rm = T))
graph_tmp <- full_join(graph_tmp, tmp, by = "Generation")

ggplot(graph_tmp, aes(x = Avg_Fec)) +
  geom_histogram(binwidth = .09) +
  geom_vline(aes(xintercept = Generation_Avg, color = "red")) +
  facet_grid(~Generation) +
  scale_y_continuous(breaks = seq(0, 15, 1)) +
  labs(x = "Centered Fecundity")

# Levene Test

leveneTest(FECUNDITY ~ as.factor(Generation), Single_2021_2022)

### Homogeneity of Variance and AVOVA/Tukey Post-hoc

ANOVA_Fec <- aov(Fecundity ~ as.factor(Generation), Single_2021_2022)
summary(ANOVA_Fec)
TukeyHSD(ANOVA_Fec)

### QQ-plot for Single Fec & Shapiro Normality (NOT NORMAL)

h <- Single_2021_2022 %>% select(FECUNDITY, Generation, Genotype) %>%  arrange(Generation)
h <- h %>% mutate(Fec = FECUNDITY)
h <- h %>% select(c("Generation", "Genotype", "Fec"))
h$Genotype <- factor(h$Genotype, levels = h$Genotype[order(h$Fec)])

par(mfrow = c(2, 3))
for (i in unique(h$Generation)) {
  p <- h %>% filter(Generation == i) %>% select(Fec)
  qqnorm(p$Fec, main = paste0("Generation ", i))
  qqline(p$Fec)
}

for(i in unique(h$Generation)){
  p <- h %>% filter(Generation == i) %>% select(Fec)
  s <- shapiro.test(p$Fec)
  print(s)
}

### Kruskal Wallis and Dunn test

kruskal.test(FECUNDITY ~ Generation, Single_2021_2022)
dunn.test(Single_2021_2022$FEC, Single_2021_2022$Generation)


### Average Fitness Distributions

graph_tmp <- Single_2021_2022 %>% group_by(Generation, Genotype, Exp_year) %>% summarise(Avg_Fit = mean(ABS_FITNESS, na.rm = T)) %>% ungroup()
tmp <- graph_tmp %>% group_by(Generation) %>% summarise(Generation_Avg = mean(Avg_Fit, na.rm = T))
graph_tmp <- full_join(graph_tmp, tmp, by = "Generation")

ggplot(graph_tmp, aes(x = Avg_Fit)) +
  geom_histogram(binwidth = .08) +
  geom_vline(aes(xintercept = Generation_Avg, color = "red")) +
  facet_grid(~Generation) +
  labs(x = "Centered Fitness")

# Levene Test

leveneTest(ABS_FITNESS ~ as.factor(Generation), Single_2021_2022)

### QQ-plot for Fitness & Shapiro Normality (NOT ALL NORMAL)

h <- Single_2021_2022 %>% select(ABS_FITNESS, Generation, Genotype) %>%  arrange(Generation)
h <- h %>% mutate(Fit = ABS_FITNESS)
h <- h %>% select(c("Generation", "Genotype", "Fit"))
h$Genotype <- factor(h$Genotype, levels = h$Genotype[order(h$Fit)])

par(mfrow = c(2, 3))
for (i in unique(h$Generation)) {
  p <- h %>% filter(Generation == i) %>% select(Fit)
  qqnorm(p$Fit, main = paste0("Generation ", i))
  qqline(p$Fit)
}

for(i in unique(h$Generation)){
  p <- h %>% filter(Generation == i) %>% select(Fit)
  s <- shapiro.test(p$Fit)
  print(s)
}

# Kruskal-Wallis Test and Dunn test

kruskal.test(ABS_FITNESS ~ Generation, Single_2021_2022)
dunn.test(Single_2021_2022$ABS_FITNESS, Single_2021_2022$Generation)

# 100 SW Distributions

graph_tmp <- Single_2021_2022 %>% group_by(Generation, Genotype, Exp_year) %>% summarise(Avg_100_SW = mean(SEED_WEIGHT_100)) %>% ungroup()
tmp <- graph_tmp %>% group_by(Generation) %>% summarise(Generation_Avg = mean(Avg_100_SW, na.rm = T))
graph_tmp <- full_join(graph_tmp, tmp, by = "Generation")

ggplot(graph_tmp, aes(x = Avg_100_SW)) +
  geom_histogram(binwidth = .07) +
  geom_vline(aes(xintercept = Generation_Avg, color = "red")) +
  facet_grid(~Generation) +
  scale_x_continuous(breaks = seq(0, 60000, 10000)) +
  theme(axis.text.x = element_text(angle = 45)) +
  labs(x = "Average 100 Seed Weight (grams)")

# Levene Test (EQUAL VARIANCE)

leveneTest(SEED_WEIGHT_100 ~ as.factor(Generation), Single_2021_2022)

### QQ-plot & Shapiro Normality for 100 SW (NORMAL)

h <- Single_2021_2022 %>% select(SEED_WEIGHT_100, Generation, Genotype) %>%  arrange(Generation)
h <- h %>% mutate(SW = SEED_WEIGHT_100)
h <- h %>% select(c("Generation", "Genotype", "SW"))
h$Genotype <- factor(h$Genotype, levels = h$Genotype[order(h$SW)])

par(mfrow = c(2, 3))
for (i in unique(h$Generation)) {
  p <- h %>% filter(Generation == i) %>% select(SW)
  qqnorm(p$SW, main = paste0("Generation ", i))
  qqline(p$SW)
}

for(i in unique(h$Generation)){
  p <- h %>% filter(Generation == i) %>% select(SW)
  s <- shapiro.test(p$SW)
  print(s)
}

# ANOVA & TUKEY for 100 SW

ANOVA_100 <- aov(SEED_WEIGHT_100 ~ as.factor(Generation), Single_2021_2022)
summary(ANOVA_100)
TukeyHSD(ANOVA_100)


