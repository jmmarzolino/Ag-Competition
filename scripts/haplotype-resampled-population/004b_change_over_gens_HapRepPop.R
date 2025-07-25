#!/usr/bin/env Rscript
#SBATCH --mem=32G
#SBATCH --time=01:00:00
#SBATCH --output=/rhome/jmarz001/bigdata/Ag-Competition/scripts/004b_change_over_gens_HapRepPop.stdout
#SBATCH -p koeniglab

setwd("/rhome/jmarz001/bigdata/Ag-Competition/data")
source("../scripts/CUSTOM_FNS.R")

library(tidyverse)
library(data.table)
library(lme4)
library(car)
library(dunn.test)
library(ggpubr)


### statistically test for trait changes between generations
# with phenotypes sampled to more accurately represent their proportion in the orignal population
hap <- fread("trait_BLUPs_HapRepPop.tsv")


###################
## plot frequency of trait values after haplotype-frequency adjustment
tmp <- hap %>% pivot_longer(cols=c(-Generation), names_to="trait", values_to="value") 
tmp$Generation <- as.factor(tmp$Generation)

g1 <- ggplot(tmp, aes(x=value, group=Generation)) + geom_histogram() + facet_wrap(~trait, scales="free_x")
ggsave("../results/haplotype_frequency_trait_vals_histogram.png", g1)

g2 <- ggplot(tmp, aes(x=value, group=Generation, fill=as.factor(Generation))) + geom_histogram() + facet_wrap(~trait, scales="free_x")
ggsave("../results/haplotype_frequency_trait_vals_histogram_colgeneration.png", g2)




#########################
# set up for normality and variance equity tests
collist <- colnames(hap)[1:ncol(hap)-1]
generationlist <- unique(hap$Generation)

pdf("../results/happop_normality_test.pdf")
normality_table <- tibble("Generation"=generationlist)

# test for normal distributions 
# w shapiro test
# and visualize w qq-plot
for(i in collist) {

  #print(i)
  tmp <- hap %>% select(c(Generation, all_of(i))) %>% tibble
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
write_delim(normality_table, "happop_blup_normality_pvals.tsv", "\t")
normality_table %>% reframe(across(where(is.numeric), \(x) x < 0.05)) %>% print
normality_table %>% reframe(across(where(is.numeric), \(x) x < 0.05)) %>% colSums %>% print
## germination has the least normal distribution(s), which is expected since it's a 0-1 decimal that is ideally near 1, not really continuously distributed trait



# summarise mean & variance
x <- hap %>% 
    group_by(Generation) %>% 
    summarise(across(where(is.numeric), list(mean= \(x) mean(x,na.rm=T), var=\(x) var(x, na.rm=T)), .names="{.col}_{.fn}")) 
write_delim(x, "happop_gens_trait_avg_var.tsv", "\t")
# write table out with generations' trait summary statistics



#############  Variance Change   #########
x2 <- hap %>% 
    group_by(Generation) %>% 
    summarise(across(where(is.numeric), \(x) var(x, na.rm=T))) 
x3 <- rbind(x2[2,]-x2[1,], x2[5,]-x2[2,])

x4 <- data.frame(t(x3))[-1,]
x4c <- rownames(x4)
x4 <- tibble(x4c, x4)
colnames(x4) <- c("trait", "F1_F18_vardiff", "F18_F58_vardiff")


##########   are var differences significant?   #########
var_tab <- tibble("trait"=collist, "all_gens"=as.numeric(NA), "F1_F18"=as.numeric(NA), "F18_F58"=as.numeric(NA)) 
# and store change in var & tests p-vals

# w Levene test
for(i in collist){
  #print(i)
  k <- leveneTest(get(i) ~ as.factor(Generation), hap)$`Pr(>F)`[1]
  tmp1 <- hap %>% filter(Generation == 1 | Generation == 18)
  tmp2 <- hap %>% filter(Generation == 18 | Generation == 58)
  m <- leveneTest(get(i) ~ as.factor(Generation), tmp1)$`Pr(>F)`[1]
  n <- leveneTest(get(i) ~ as.factor(Generation), tmp2)$`Pr(>F)`[1]

  # save values to table in correct trait-row
  var_tab[which(var_tab$trait == i), 2] <- k
  var_tab[which(var_tab$trait == i), 3] <- m
  var_tab[which(var_tab$trait == i), 4] <- n
}


# join var change table with p-val table
var_tab_sig <- full_join(var_tab, x4, by="trait")
# save a copy of variance change + significance
write_delim(var_tab_sig, "happop_traits_var_change_sig_test.tsv", "\t")
# remove 'all_gens' comparison test before spivoting, for clarity
var_tab_sig <- var_tab_sig %>% select(-all_gens)

## format table w pivot
var_tab_sig <- var_tab_sig %>% pivot_longer(cols=c(F1_F18, F18_F58), names_to="generations_compared", values_to="Levene_pval")
var_tab_sig$Levene_sig <- var_tab_sig$Levene_pval < 0.05


tmp <- var_tab_sig %>% 
  pivot_longer(cols=c(F1_F18_vardiff, F18_F58_vardiff), names_to="gensvardiff", values_to="var_difference") 

tmp$gensvardiff <- gsub("_vardiff", "", tmp$gensvardiff)
tmp <- tmp %>% filter(generations_compared == gensvardiff) %>% select(-c(gensvardiff))

write_delim(tmp, "happop_traits_var_diff_btwn_gens_sig.tsv", "\t")






#####################
#### signif changes in generation mean
hap <- tibble(hap)
hap$Generation <- as.factor(hap$Generation)
levels(hap$Generation)

######
# is there a difference between using anova & kruskal test on these phenotypes?
for(m in 1:(ncol(hap)-1)) {
  print(colnames(hap)[m])
  kruskal.test(unlist(hap[,m]) ~ Generation, hap) %>% print

  aov(unlist(hap[,m]) ~ as.factor(Generation), hap) %>% summary %>% print
}

#####
## do phenotypes signif.ly change from parents to F18?
## do phenotypes signif.ly change from F18 compared to F58?

# set up storing results
signif_tab <- tibble("trait"=collist, "F1_F18"=as.numeric(NA), "F18_F58"=as.numeric(NA)) 

## first, test for significant result w Kruskal-wallis test

for(i in 1:length(signif_tab$trait)){
  # filter to only relevant generations
  tmp1 <- hap %>% select(c(all_of(i), Generation)) %>% filter(Generation == 18 | Generation == 1)
  tmp2 <- hap %>% select(c(all_of(i), Generation)) %>% filter(Generation == 18 | Generation == 58)

  res1 <- kruskal.test(unlist(tmp1[,1]) ~ Generation, tmp1)
  res2 <- kruskal.test(unlist(tmp2[,1]) ~ Generation, tmp2)
  signif_tab[i, 2] <- res1$p.value
  signif_tab[i, 3] <- res2$p.value
}



## for traits w signif K-W p-val, run posthoc test & store results
signif_tab_KW <- signif_tab %>% pivot_longer(cols=c(F1_F18, F18_F58), names_to="generations_compared", values_to="Kruskal_pval")
signif_tab_KW$Kruskal_sig <- signif_tab_KW$Kruskal_pval < 0.05

# set up storing posthoc test results
signif_tab$F1_F18_pval <- as.numeric(NA) 
signif_tab$F18_F58_pval <- as.numeric(NA) 

# one posthoc test for each trait that requires it
smth <- signif_tab_KW[which(signif_tab_KW$Kruskal_sig), 1:2]
trt <- unique(unlist(smth[,1]))

for( j in 1:length(trt))  {

    test_trt <- trt[j]
    test_trt_joined_happops <- hap %>% select(c(all_of(test_trt), Generation))

    rownum <- which(signif_tab$trait == test_trt)

    dunres <- dunn.test(unlist(test_trt_joined_happops[,1]), test_trt_joined_happops$Generation)

    # extract posthoc values via position in posthoc table
    #dunres$comparisons[1]
    p18_pval <- dunres$P[1]
    #dunres$comparisons[8]
    f1858_pval <- dunres$P[8]

    # store results in table
    signif_tab[rownum, 4] <- p18_pval
    signif_tab[rownum, 5] <- f1858_pval
}


# claculate trait change between generations of interest
x2 <- hap %>% 
    group_by(Generation) %>% 
    summarise(across(where(is.numeric), \(x) mean(x, na.rm=T))) 
x3 <- rbind(x2[2,]-x2[1,], x2[5,]-x2[2,])

x4 <- data.frame(t(x3))[-1,]
x4c <- rownames(x4)
x4 <- tibble(x4c, x4)
colnames(x4) <- c("trait", "F1_F18_meandiff", "F18_F58_meandiff")

# and join with p-val table
signif_tab <- full_join(signif_tab, x4, by="trait")


## format table by pivoting 
signif_tab2 <- signif_tab %>% pivot_longer(cols=c(F1_F18, F18_F58), names_to="generations_compared", values_to="Kruskal_pval")
signif_tab2$Kruskal_sig <- signif_tab2$Kruskal_pval < 0.05


tmp <- signif_tab2 %>% 
  pivot_longer(cols=c(F1_F18_meandiff, F18_F58_meandiff), names_to="gensmeandiff", values_to="mean_difference") %>%
  pivot_longer(cols=c(F1_F18_pval, F18_F58_pval), names_to="genspostpval", values_to="Dunn_posthoc_pval")

tmp$gensmeandiff <- gsub("_meandiff", "", tmp$gensmeandiff)
tmp$genspostpval <- gsub("_pval", "", tmp$genspostpval)

tmp <- tmp %>% filter(generations_compared == gensmeandiff & gensmeandiff == genspostpval) %>% select(-c(gensmeandiff, genspostpval))


# mark significant dunn p-values based on adjusted alpha value
tmp$Dunn_posthoc_sig <- tmp$Dunn_posthoc_pval < (0.05 / 2)

write_delim(tmp, "happop_blups_anova_posthoc_results.tsv", "\t")
