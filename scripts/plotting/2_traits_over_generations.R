#!/usr/bin/env Rscript
#SBATCH --ntasks=1
#SBATCH --mem=30G
#SBATCH --time=02:00:00
#SBATCH --output=/rhome/jmarz001/bigdata/Ag-Competition/scripts/2_traits_over_generations.stdout
#SBATCH -p koeniglab

# This script plots scatterplots with linear regressions for each trait in the single subpopulation, as well as distributions for these traits

library(tidyverse)
library(ggpubr)
library(ggplot2)
library(RColorBrewer)

setwd("/bigdata/koeniglab/jmarz001/Ag-Competition")
source("scripts/CUSTOM_FNS.R")
df <- read_delim("data/trait_BLUPs.tsv")
df <- add_generation(df)
df <- df %>% filter(Generation != 50)
#df$Generation <- as.factor(df$Generation)

# the lowest color value is too light, so adjust the color scale down one
#display.brewer.pal(6, "Blues")
adjusted_blues <- brewer.pal(7, "Blues")[3:7]


### PLOTTING
## arrange data for facet plotting
df_long <- df %>% pivot_longer(cols=-c('Genotype', 'Generation'), values_to="VALUE", names_to="TRAIT")

# substitute trait names w/ tidy text versions
df_long$TRAIT <- gsub("FT_blup", "Flowering Time BLUP", df_long$TRAIT)
df_long$TRAIT <- gsub("TOTAL_MASS_blup", "Total Seed Mass BLUP", df_long$TRAIT)
df_long$TRAIT <- gsub("GERMINATION_blup", "Germination BLUP", df_long$TRAIT)
df_long$TRAIT <- gsub("SEED_WEIGHT_100_blup", "100-Seed Weight BLUP", df_long$TRAIT)
df_long$TRAIT <- gsub("FECUNDITY_blup", "Fecundity BLUP", df_long$TRAIT)
df_long$TRAIT <- gsub("FITNESS_blup", "Fitness BLUP", df_long$TRAIT)


# check normality & plot trait distributions
g <- ggplot(df_long, aes(VALUE)) +
  geom_density(color="#084594", linewidth=1) +
  facet_wrap(~TRAIT, scales="free") +  
  labs(x="", y="density") +
  theme_bw(base_size=20) #+
  #stat_summary(fun = median, geom = "vline", orientation = "y", aes(xintercept = after_stat(x), y = 0), color="#eca50b", linewidth=1) 
ggsave("results/trait_distributions.png", g, width=12)


# one density line per generation
g <- ggplot(df_long, aes(VALUE, group=Generation, color=as.factor(Generation))) +
  geom_density(linewidth=1) +
  scale_color_manual(values=adjusted_blues, name="Generation") + 
  facet_wrap(~TRAIT, scales="free") +
  theme_bw(base_size=20) +
  labs(x="", y="density") +
  stat_summary(fun = median, geom = "vline", orientation = "y", aes(xintercept = after_stat(x), y = 0), linewidth=1) 
ggsave("results/trait_distributions_Wgeneration.png", g, width=16)



## boxplots comparing traits over generations
ggplot(single, aes(y=FT, x=as.factor(Generation), group=as.factor(Generation), fill=as.factor(Generation))) +
    geom_boxplot() +
    theme_minimal() +
    facet_wrap(~trait, scales="free")

ggplot(single, aes(y=FT, x=as.factor(Generation), group=as.factor(Generation), fill=as.factor(Generation))) +
    geom_boxplot() +
    theme_minimal() +
    facet_wrap(~trait, scales="free")







## plotting traits over generations
# w standardized data

s_fit <- df %>% mutate(stand_fit = scale(FECUNDITY),
                          stand_fec = scale(FITNESS),
                          stand_tw = scale(TOTAL_MASS),
                          stand_100 = scale(SEED_WEIGHT_100),
                          stand_ft = scale(FT))

### 02_standard_fit_over_gen

a <- ggplot(s_fit, aes(Generation, stand_fit, add = "reg.line")) +
geom_jitter(alpha = .5) +
  geom_smooth(method = lm) +
  geom_boxplot(aes(Generation, stand_fit, group = Generation), width = 1.5, alpha = .5) +
  stat_regline_equation() +
  scale_y_continuous(breaks = seq(-4, 4, 2), limits = c(-4,4)) +
  labs(x = "Generation",
       y = "Average Fitness") +
  theme_bw()

### 02_standard_fec_over_gen

b <- ggplot(s_fit, aes(Generation, stand_fec, add = "reg.line")) +
  geom_jitter(alpha = .5) +
  geom_smooth(method = lm) +
  geom_boxplot(aes(Generation, stand_fec, group = Generation), width = 1.5, alpha = .5) +
  stat_regline_equation() +
  labs(x = "Generation",
       y = "Average Fecundity") +
  scale_y_continuous(breaks = seq(-4, 4, 2), limits = c(-4,4)) +
  theme_bw()

### 02_standard_fit_over_gen

c <- ggplot(s_fit, aes(Generation, stand_ft, add = "reg.line")) +
  geom_jitter(alpha = .5) +
  geom_smooth(method = lm) +
  stat_regline_equation() +
  geom_boxplot(aes(Generation, stand_ft, group = Generation), width = 1.5, alpha = .5) +
  labs(x = "Generation",
       y = "Average Flowering Time (Days)") +
  scale_y_continuous(breaks = seq(-4, 4,2), limits = c(-4,4)) +
  theme_bw()

### 02_standard_tw_over_gen

d <- ggplot(s_fit, aes(Generation, stand_tw, add = "reg.line")) +
  geom_jitter(alpha = .5) +
  geom_smooth(method = lm) +
  geom_boxplot(aes(Generation, stand_tw, group = Generation), width = 1.5, alpha = .5)+
  stat_regline_equation() +
  labs(x = "Generation",
       y = "Average Total Seed Weight (g)") +
  scale_y_continuous(breaks = seq(-4, 4,2), limits = c(-4,4)) +
  theme_bw()

### 02_standard_100sw_over_gen

e <- ggplot(s_fit, aes(Generation, stand_100, add = "reg.line")) +
  geom_jitter(alpha = .5) +
  geom_smooth(method = lm) +
  geom_boxplot(aes(Generation, stand_100, group = Generation), width = 1.5, alpha = .5)+
  stat_regline_equation() +
  scale_y_continuous(breaks = seq(-4, 4, 2), limits = c(-4,4)) +
  labs(x = "Generation",
       y = "Average 100 Seed Weight (g)") +
  theme_bw()

### combine scatterplots & save
y <- ggarrange(a, b, c, d, e, top = "Evolution of Our Four Measured Phenotypes", nrow = 2, ncol = 3)
ggsave("/bigdata/koeniglab/jmarz001/Ag-Competition/results/traits_over_generations_scatterplots_scaled.png",y, width = 14, height = 10)



## Unscaled data


# make this a loop instead of a function
traitlist <- c('FITNESS', 'FECUNDITY', 'FT', '...')
labellist <- c("Average Fitness", "Average Fecundity")

for(i in traitlist) {


unscaled_graphs <- function(x = df, y = trait){
    # set label per plot based on graphed trait
  l <- paste("Average", y)


  ggplot(x, aes(Generation, y)) +
    geom_jitter() +
    geom_boxplot(aes(group = Generation), width = 1.5, alpha = .5) +
    theme_bw() 

}



  x %>% ggplot(aes(Generation, get(y))) +
    geom_point() +
    labs(x="Generation",
          y = l)
}


plot_list <- list(c(3,2))
loop_over <- names(fitness_df)[c(3,4,5)]
for (i in loop_over){
  a <- ggplot(fitness_df, aes(Generation, fitness_df[i])) + geom_point()
  
}

unscaled_graphs(sin, FECUNDITY)



### 02a_fit_over_gen

a <- ggplot(sin, aes(Generation, FITNESS, add = "reg.line")) +
  geom_jitter(alpha = .5) +
  geom_smooth(method = lm) +
  #geom_hline(aes(yintercept = mean(FITNESS)), color = "red") +
  geom_boxplot(aes(Generation, FITNESS, group = Generation), width = 1.5, alpha = .5) +
  stat_regline_equation() +
  labs(x = "Generation",
       y = "Average Fitness") +
  theme_bw()

### 02a_fec_over_gen

b <- ggplot(sin, aes(Generation, FECUNDITY, add = "reg.line")) +
  geom_jitter(alpha = .5) +
  geom_smooth(method = lm) +
  #geom_hline(aes(yintercept = mean(FECUNDITY)), color = "red") +
  geom_boxplot(aes(Generation, FECUNDITY, group = Generation), width = 1.5, alpha = .5) +
  stat_regline_equation() +
  labs(x = "Generation",
       y = "Average Fecundity") +
  theme_bw()

### 02a_ft_over_gen

c <- ggplot(sin, aes(Generation, FT, add = "reg.line")) +
  geom_jitter(alpha = .5) +
  geom_smooth(method = lm) +
  #geom_hline(aes(yintercept = mean(FT)), color = "red") +
  stat_regline_equation() +
  geom_boxplot(aes(Generation, FT, group = Generation), width = 1.5, alpha = .5) +
  labs(x = "Generation",
       y = "Average Flowering Time (Days)") +
  theme_bw()

### 02a_tw_over_gen

d <- ggplot(sin, aes(Generation, TOTAL_MASS, add = "reg.line")) +
  geom_jitter(alpha = .5) +
  geom_smooth(method = lm) +
  #geom_hline(aes(yintercept = mean(TOTAL_MASS)), color = "red") +
  geom_boxplot(aes(Generation, TOTAL_MASS, group = Generation), width = 1.5, alpha = .5)+
  stat_regline_equation() +
  labs(x = "Generation",
       y = "Average Total Seed Weight (g)") +
  theme_bw()

### 02a_100sw_over_gen

e <- ggplot(sin, aes(Generation, SEED_WEIGHT_100, add = "reg.line")) +
  geom_jitter(alpha = .5) +
  geom_smooth(method = lm) +
  #geom_hline(aes(yintercept = mean(SEED_WEIGHT_100)), color = "red") +
  geom_boxplot(aes(Generation, SEED_WEIGHT_100, group = Generation), width = 1.5, alpha = .5)+
  stat_regline_equation() +
  labs(x = "Generation",
       y = "Average 100 Seed Weight (g)") +
  theme_bw()

### 02_combined_single_evolution_scatterplots

gg <- ggarrange(a, b, c, d, e) #, top = "Evolution of Our Four Measured Phenotypes", nrow = 2, ncol = 3)
ggsave("/bigdata/koeniglab/jmarz001/Ag-Competition/results/traits_over_generations_scatterplot.png", gg, width = 14, height = 10)

list_df <- c(TOTAL_MASS, SEED_WEIGHT_100)



for (i in colnames(fitness_df)[c(3,4,5,7,8)]){
  p <- ggplot(fitness_df, aes(Generation, i))
  print(p)
}

### 02b_Fecundity_Distributions_Over_Generations

ggplot(fitness_df, aes(x = FECUNDITY, color = Condition, fill = Condition)) +
  geom_histogram(position = "identity", bins = 45, alpha = .3) +
  labs(x = "Fecundity",
       y = "Frequency",
       title = "Fecundity Over Generations") +
  facet_wrap(~Generation, scales = "free_x")
ggsave("/bigdata/koeniglab/jmarz001/Ag-Competition/results/02b_Fecundity_Over_Generations_Distributions.png")

Average_Haplo_rep$Generation <- as.factor(Average_Haplo_rep$Generation)
ggplot(Average_Haplo_rep, aes(x = Fecundity, fill = Generation, group = Generation)) +
  geom_histogram(alpha = .5, position = "identity", binwidth = 70) +
  scale_fill_brewer(palette = "Blues")
ggsave("/bigdata/koeniglab/jmarz001/Ag-Competition/results/02b_Overlapping_Distributions_Fecundity_Over_Generations.png")

### 02bi_100_SW_Distributions

ggplot(Average_Haplo_rep, aes(x = `100 seed weight`, group = Generation, fill = Generation)) +
  geom_histogram(alpha =.5, position = "identity") +
  scale_fill_brewer(palette = "Blues")
ggsave("/bigdata/koeniglab/jmarz001/Ag-Competition/results/02bi_Overlapping_Distributions_for_Generational_Change_100SW_Distributions.png")








AT_fitness <- sin %>% filter(Genotype == "48_5")
AT_fitness <- AT_fitness$FITNESS

a + geom_hline(aes(yintercept=AT_fitness))

### 02c_Single_Relative_Fitness_to_Atlas

sf <- ggplot(sin, aes(x=Generation, y=FITNESS


, add = "reg.line")) +
  geom_jitter() + #alpha = .4) +
  geom_boxplot(aes(x=Generation, y=FITNESS)) +
  stat_regline_equation() + #label.x = 40
  #ylim(0, 40000) +
  geom_smooth(method = "lm") +
  geom_hline(yintercept = AT_fitness, color = "red") +
  labs(title = "Relative Fitness of Atlas Compared to Single Condition") 
ggsave("/bigdata/koeniglab/jmarz001/Ag-Competition/results/02c_Single_Relative_Fitness_to_Atlas.png")





### 02ci_Mixed_Relative_Fitness_to_Atlas

mf <- ggplot(Rep_Mixed, aes(Generation, FITNESS, add = "reg.line")) +
  geom_jitter(alpha = .4) +
  ylim(0, 40000) +
  stat_regline_equation() +
  geom_hline(yintercept = 21347.22, color = "red") +
  geom_smooth(method = "lm") +
  labs(title = "Relative Fitness of Atlas Compared to Mixed Condition")
ggsave("/bigdata/koeniglab/jmarz001/Ag-Competition/results/02ci_Mixed_Relative_Fitness_to_Atlas.png")

jf <- ggplot(Average_Haplo_rep, aes(Generation, FITNESS, add = "reg.line", color = Condition)) +
  geom_jitter(alpha = .6) +
  stat_regline_equation() +
  geom_smooth(method = "lm" ) +
  geom_hline(yintercept = 21347.22, color = "red") +
  labs(title = "Relative Fitness of Atlas Compared to Both Conditions")
ggsave("/bigdata/koeniglab/jmarz001/Ag-Competition/results/02cii_Relative_Fitness_to_Atlas_Combined.png", width = 12, height = 7)

### 02cii_Relative_Fitness_to_Atlas_Combined

g <- ggarrange(sf, mf, jf, nrow = 2)
ggsave("/bigdata/koeniglab/jmarz001/Ag-Competition/results/02cii_Relative_Fitness_to_Atlas_Combined.png", g, width = 12, height = 7)
