library(tidyverse)
library(readr)
library(dplyr)
library(ggpubr)
library(ggplot2)
library(tidyr)
library(car)
library(gridExtra)

### Combining Years 2021-2022 and 2022-2023 into one dataframe for plotting (Single)

tmp <- Single_2021_2022
colnames(tmp)[4:19] <- paste0(colnames(tmp)[4:19], sep = "_", "2021_2022")
cmp <- Rep_Single
Combined_Single_Years <- full_join(cmp, tmp)

### Reformatting the dataframe to be better suited at distributions

tmp <- Single_2021_2022 
cmp <- Rep_Single %>% select(!c("Pedigree", "Numbers")) %>% mutate(Year = "2") 
correct <- c("Genotype", "Condition", "Generation", "Brown Bag Weight", "100 seed weight", "Fecundity", "Fitness", "FT", "Haplotype", "Atlas_Avg_Fec", "Atlas_Avg_Fitness", "Atlas_Avg_Total_Weight", "Centered_Fit", "Centered_FT", "Centered_Fec", "Centered_TW", "Exp_Fec_Per_Plant",
             "Exp_Fit_Per_Plant", "Exp_TW_Per_Plant", "Year")
cmp2 <- cmp[, correct]
tmp <- tmp %>% mutate(Year = "1")
names(tmp) <- names(cmp2)
Combined_Single_Years_Distribution <- rbind(cmp, tmp)
Combined_Single_Years_Distribution <- Combined_Single_Years_Distribution %>% mutate(Average_TW = mean(`Brown Bag Weight`),
                                                                                                       Average_FT = mean(FT),
                                                                                                       Average_Fec = mean(Fecundity, na.rm = TRUE),
                                                                                                       Average_Fit = mean(Fitness, na.rm = TRUE)) %>% group_by(Year)

### Correlation graphs of the years for each phenotype

### TW
ggplot(Combined_Single_Years, aes(TOTAL_MASS_2021_2022, `Brown Bag Weight`), add = "reg.line") +
  geom_point() +
  stat_cor(label.y = 170) +
  geom_smooth(method = 'lm') +
  geom_abline(slope =1, intercept= 0, color = "red") +
  labs(x = "Year 1 Total Seed Weight (g)",
       y = "Year 2 Total Seed Weight (g)")

### Flowering
ggplot(Combined_Single_Years, aes(FT_2021_2022, FT), add = "reg.line") +
  geom_point() +
  stat_cor(label.y = 120) +
  geom_smooth(method = 'lm') +
  geom_abline(slope = 1, intercept = 0, color = 'red') +
  labs(x = "Year 1 Days to Flower",
       y = "Year 2 Days to Flower")

### Fecundity
ggplot(Combined_Single_Years, aes(Fecundity_2021_2022, Fecundity), add = "reg.line") +
  geom_point() +
  stat_cor(label.y = 500, label.x = 3000) +
  geom_smooth(method = 'lm') +
  geom_abline(slope = 1, intercept = 0, color = 'red') +
  labs(x = "Year 1 Fecundity",
       y = "Year 2 Fecundity")
  
### Fitness

ggplot(Combined_Single_Years, aes(Fitness_2021_2022, Fitness), add = "reg.line") +
  geom_point() +
  stat_cor(label.y = 10000, label.x = 30000) +
  geom_smooth(method = 'lm') +
  geom_abline(slope = 1, intercept = 0, color = 'red') +
  labs(x = "Year 1 Fitness",
       y = "Year 2 Fitness")

### Distributions Between 2021-2022 and 2022-2023 seasons for all phenotypes

### TW Distribution

ggplot(Combined_Single_Years_Distribution, aes(x = `Brown Bag Weight`, fill = Year)) +
  geom_density(alpha = .5) +
  labs(x = "Average Total Seed Weight") +
  geom_vline(xintercept = Combined_Single_Years_Distribution$Average_TW, color = 'red')

### FT Distribution

ggplot(Combined_Single_Years_Distribution, aes(x = FT, fill = Year)) +
  geom_density(alpha = .5) +
  labs(x = "Average Days To Flower") +
  geom_vline(xintercept = Combined_Single_Years_Distribution$Average_FT, color = 'red')

### Fec Distribution

ggplot(Combined_Single_Years_Distribution, aes(x = Fecundity, fill = Year)) +
  geom_density(alpha = .5) +
  labs(x = "Average Fecundity") +
  geom_vline(xintercept = Combined_Single_Years_Distribution$Average_Fec, color = 'red')

### Fit Distributions 

ggplot(Combined_Single_Years_Distribution, aes(x = Fitness, fill = Year)) +
  geom_density(alpha = .5) +
  labs(x = "Average Fitness") +
  geom_vline(xintercept = Combined_Single_Years_Distribution$Average_Fit, color = 'red')

