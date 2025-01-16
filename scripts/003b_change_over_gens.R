#!/usr/bin/env Rscript
#SBATCH --mem=50G
#SBATCH --time=02:00:00
#SBATCH --output=/rhome/jmarz001/bigdata/Ag-Competition/scripts/003b_change_over_gens.stdout
#SBATCH -p koeniglab

setwd("/rhome/jmarz001/bigdata/Ag-Competition/data")
source("../scripts/CUSTOM_FNS.R")

library(tidyverse)
library(data.table)
library(lme4)
library(car)
library(dunn.test)


# read in filtered but un-transformed data
df <- fread("DERIVED_PHENOTYPES.tsv")
df <- df %>% 
        filter(Condition == "single") %>% select(-c(Condition, Replicate, SEED_COUNT)) %>%
        group_by(Genotype, Exp_year) %>% summarise(across(where(is.numeric), mean)) %>%
        ungroup() %>% select(-Exp_year) %>% 
        group_by(Genotype) %>% summarise(across(where(is.numeric), mean)) %>%
        ungroup() 

dfb <- fread("trait_BLUPs.tsv")
df <- full_join(df, dfb, by=c('Genotype'), suffix=c("", "_blup"))
df <- add_generation(df)
df <- df %>% select(c('Generation', 'FT', 'TOTAL_MASS', 'SEED_WEIGHT_100', 'GERMINATION', 'FECUNDITY', 'FITNESS', 'FT_blup', 'TOTAL_MASS_blup', 'GERMINATION_blup', 'SEED_WEIGHT_100_blup', 'FECUNDITY_blup', 'FITNESS_blup'))


## BASE STATISTICS
# summarise mean & variance
x <- df %>% 
    group_by(Generation) %>% 
    summarise(across(where(is.numeric), list(mean=mean, var=var), .names="{.col}_{.fn}")) 
write_delim(x, "generations_trait_avg_var.tsv", "\t")
# write table out with generations' trait summary statistics


# filter traits with no variance
variant_phenos <- names(which(apply(X=df[,2:ncol(df)], 2, var) > 0))
df <- df %>% select(c(Generation, all_of(variant_phenos)))


# set up for normality and variance equity tests
collist <- colnames(df)[2:ncol(df)]
generationlist <- x$Generation


pdf("../results/normality_test.pdf")

normality_table <- tibble("Generation"=generationlist)

# test for normal distributions 
# w shapiro test
# and visualize w qq-plot
for(i in collist) {

  #print(i)
  tmp <- df %>% select(c(Generation, all_of(i))) %>% tibble
  # empty list for p-values
  val <- c()

  for(g in generationlist) {

    #print(g)
    tmp_gen <- tmp %>% filter(Generation == g)
    s <- shapiro.test(tmp_gen[,2][[1]])
    
    # add p-val to list
    val <- c(val, s$p.value)
    qqnorm(tmp_gen[,2][[1]], main = paste0(i, ": ", "Generation ", g))
    qqline(tmp_gen[,2][[1]])
  }

  normality_table <- cbind(normality_table, val)
}

dev.off()

colnames(normality_table)[2:ncol(normality_table)] <- paste0(collist, "_normality_pval")
write_delim(normality_table, "trait_and_blup_normality_pvals.tsv", "\t")
normality_table %>% reframe(across(where(is.numeric), \(x) x < 0.05)) %>% print
normality_table %>% reframe(across(where(is.numeric), \(x) x < 0.05)) %>% colSums %>% print
## germination has the least normal distribution(s), which is expected since it's a 0-1 decimal that is ideally near 1, not really continuously distributed trait


# test for equality of variance between groups before anova
# w Levene test
for(i in collist){
  print(i)
  leveneTest(get(i) ~ as.factor(Generation), df) %>% print
}
## only flowering time has unequal variance between generations





#### signif changes in generation mean
df$Generation <- as.factor(df$Generation)
levels(df$Generation)


# is there a difference between using anova & kruskal test on these phenotypes?
for(m in 2:ncol(df)) {
  print(colnames(df[,m]))
  kruskal.test(unlist(df[,m]) ~ Generation, df) %>% print

  aov(unlist(df[,m]) ~ as.factor(Generation), df) %>% summary %>% print
}






### limit statistics of signif mean differences to blup phenotypes
df <- df %>% select(c(contains("blup"), Generation))
trait_list <- colnames(df)[1:4]

## do phenotypes signif.ly change from parents to F18?
## do phenotypes signif.ly change from F18 compared to F58?

# set up storing results
signif_tab <- tibble("trait"=trait_list, "P_F18"=as.numeric(NA), "F18_F58"=as.numeric(NA)) 

## first, test for significant result w Kruskal-wallis test

for(i in 1:length(signif_tab$trait)){
  # filter to only relevant generations
  tmp1 <- df %>% select(c(all_of(i), Generation)) %>% filter(Generation == 18 | Generation == 0)
  tmp2 <- df %>% select(c(all_of(i), Generation)) %>% filter(Generation == 18 | Generation == 58)

  res1 <- kruskal.test(unlist(tmp1[,1]) ~ Generation, tmp1)
  res2 <- kruskal.test(unlist(tmp2[,1]) ~ Generation, tmp2)
  signif_tab[i, 2] <- res1$p.value
  signif_tab[i, 3] <- res2$p.value
}



## for traits w signif K-W p-val, run posthoc test & store results
signif_tab_KW <- signif_tab %>% pivot_longer(cols=c(P_F18, F18_F58), names_to="generations_compared", values_to="Kruskal_pval")
signif_tab_KW$Kruskal_sig <- signif_tab_KW$Kruskal_pval < 0.05



# set up storing posthoc test results
signif_tab$P_F18_meandiff <- as.numeric(NA) 
signif_tab$F18_F58_meandiff <- as.numeric(NA) 
signif_tab$P_F18_pval <- as.numeric(NA) 
signif_tab$F18_F58_pval <- as.numeric(NA) 



# one posthoc test for each trait that requires it
smth <- signif_tab_KW[which(signif_tab_KW$Kruskal_sig), 1:2]

trt <- unique(unlist(smth[,1]))

for( j in 1:length(trt))  {

    test_trt <- trt[j]
    test_trt_df <- df %>% select(c(all_of(test_trt), Generation))

    rownum <- which(signif_tab$trait == test_trt)

    dunres <- dunn.test(unlist(test_trt_df[,1]), test_trt_df$Generation)

    # extract posthoc values via position in posthoc table
    #dunres$comparisons[1]
    p18_meandiff <- -1*(dunres$Z[1])
    p18_pval <- dunres$P[1]
    #dunres$comparisons[8]
    f1858_meandiff <- -1*(dunres$Z[8])
    f1858_pval <- dunres$P[8]

    # store results in table
    signif_tab[rownum, 4] <- p18_meandiff
    signif_tab[rownum, 5] <- f1858_meandiff
    signif_tab[rownum, 6] <- p18_pval
    signif_tab[rownum, 7] <- f1858_pval

}


## format table by pivoting 
signif_tab2 <- signif_tab %>% pivot_longer(cols=c(P_F18, F18_F58), names_to="generations_compared", values_to="Kruskal_pval")

signif_tab2$Kruskal_sig <- signif_tab2$Kruskal_pval < 0.05


tmp <- signif_tab2 %>% 
  pivot_longer(cols=c(P_F18_meandiff, F18_F58_meandiff), names_to="gensmeandiff", values_to="mean_difference") %>%
  pivot_longer(cols=c(P_F18_pval, F18_F58_pval), names_to="genspostpval", values_to="Dunn_posthoc_pval")

tmp$gensmeandiff <- gsub("_meandiff", "", tmp$gensmeandiff)
tmp$genspostpval <- gsub("_pval", "", tmp$genspostpval)

tmp <- tmp %>% filter(generations_compared == gensmeandiff & gensmeandiff == genspostpval) %>% select(-c(gensmeandiff, genspostpval))



write_delim(tmp, "blups_anova_posthoc_results.tsv", "\t")



