#!/usr/bin/env Rscript
#SBATCH --mem=50G
#SBATCH --time=02:00:00
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

hap <- fread("hap_assign.txt")
df <- fread("trait_BLUPs.tsv")


hap$Haplotype <- as.factor(hap$Haplotype)
# cut 3-digit genotype codes down to 2 digits
# to match family lines
df$Genotype <- gsub("(\\d+_\\d+)_\\d", "\\1", df$Genotype)
grep("\\d+_\\d+_\\d+", df$Genotype)


# how common was each haplotype in the 4 sampled generations?
## calculate haplotype frequency per generation
f18_hap <- hap[grep("^1_\\d+", hap$Genotype),]
f18_hap_table <- data.frame(table(f18_hap$Haplotype))
colnames(f18_hap_table) <- c("Haplotype", "Frequency")
f18_hap_table$fraction <- f18_hap_table$Frequency / sum(f18_hap_table$Frequency)
f18_hap_table <- f18_hap_table %>% filter(Frequency > 0)

f28_hap <- hap[grep("^2_\\d+", hap$Genotype),]
f28_hap_table <- data.frame(table(f28_hap$Haplotype))
colnames(f28_hap_table) <- c("Haplotype", "Frequency")
f28_hap_table$fraction <- f28_hap_table$Frequency / sum(f28_hap_table$Frequency)
f28_hap_table <- f28_hap_table %>% filter(Frequency > 0)

f50_hap <- hap[grep("^3_\\d+", hap$Genotype),]
f50_hap_table <- data.frame(table(f50_hap$Haplotype))
colnames(f50_hap_table) <- c("Haplotype", "Frequency")
f50_hap_table$fraction <- f50_hap_table$Frequency / sum(f50_hap_table$Frequency)
f50_hap_table <- f50_hap_table %>% filter(Frequency > 0)

f58_hap <- hap[grep("^7_\\d+", hap$Genotype),]
f58_hap_table <- data.frame(table(f58_hap$Haplotype))
colnames(f58_hap_table) <- c("Haplotype", "Frequency")
f58_hap_table$fraction <- f58_hap_table$Frequency / sum(f58_hap_table$Frequency)
f58_hap_table <- f58_hap_table %>% filter(Frequency > 0)

# combine genotype, phenotype, and haplotype
hap_join <- inner_join(hap, df, by = "Genotype")


## calculate generations' trait averages weighted by haplotype

# average phenotype for each haplotype
hap_trait_avg <- hap_join %>% group_by(Haplotype) %>% summarise(across(where(is.numeric), mean)) 

# join haplotype, generation frequency, and haplotype phenotypes
f18_hap_join <- right_join(hap_trait_avg, f18_hap_table, by="Haplotype")
f28_hap_join <- right_join(hap_trait_avg, f28_hap_table, by="Haplotype")
f50_hap_join <- right_join(hap_trait_avg, f50_hap_table, by="Haplotype")
f58_hap_join <- right_join(hap_trait_avg, f58_hap_table, by="Haplotype")


# copy rows to reflect their generation frequency number

for(j in c("f18_hap_join", "f28_hap_join", "f50_hap_join", "f58_hap_join")) {
  x <- get(j)
  out_df <- tibble(.rows=sum(x$Frequency))
  for(i in c(2:7)) {
    p <- rep((x[,i][[1]]), x$Frequency)
    out_df <- cbind(out_df, p)
  }
  colnames(out_df) <- colnames(x)[c(2:7)]
  assign(paste0(j, "_poprepd"), out_df)
}


# add col indicating phenotype's generation
f18_hap_join_poprepd$Generation <- 18
f28_hap_join_poprepd$Generation <- 28
f50_hap_join_poprepd$Generation <- 50
f58_hap_join_poprepd$Generation <- 58
## represent parental generation phenotypes at frequency
## with one entry for each haplotype
#table(hap_trait_avg$Haplotype)
f0_hap_join_poprepd <- hap_trait_avg %>% select(-Haplotype)
f0_hap_join_poprepd$Generation <- 0



x1 <- rbind(f18_hap_join_poprepd, f28_hap_join_poprepd)
x2 <- rbind(f50_hap_join_poprepd, f58_hap_join_poprepd)
joined_happops <- rbind(x1,x2)
joined_happops <- rbind(joined_happops,f0_hap_join_poprepd)


###################
## plot frequency of trait values after haplotype-frequency adjustment
tmp <- joined_happops %>% pivot_longer(cols=c(-Generation), names_to="trait", values_to="value") 
tmp$Generation <- as.factor(tmp$Generation)

g1 <- ggplot(tmp, aes(x=value, group=Generation)) + geom_histogram() + facet_wrap(~trait, scales="free_x")
ggsave("../results/haplotype_frequency_trait_vals_histogram.png", g1)

g2 <- ggplot(tmp, aes(x=value, group=Generation, fill=as.factor(Generation))) + geom_histogram() + facet_wrap(~trait, scales="free_x")
ggsave("../results/haplotype_frequency_trait_vals_histogram_colgeneration.png", g2)


## BASE STATISTICS
# summarise mean & variance
x <- joined_happops %>% 
    group_by(Generation) %>% 
    summarise(across(where(is.numeric), list(mean= \(x) mean(x,na.rm=T), var=\(x) var(x, na.rm=T)), .names="{.col}_{.fn}")) 
write_delim(x, "happop_gens_trait_avg_var.tsv", "\t")
# write table out with generations' trait summary statistics



#########################
# set up for normality and variance equity tests
collist <- colnames(joined_happops)[1:6]
generationlist <- x$Generation


pdf("../results/normality_test.pdf")

normality_table <- tibble("Generation"=generationlist)

# test for normal distributions 
# w shapiro test
# and visualize w qq-plot
for(i in collist) {

  #print(i)
  tmp <- joined_happops %>% select(c(Generation, all_of(i))) %>% tibble
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
write_delim(normality_table, "blup_normality_pvals.tsv", "\t")
normality_table %>% reframe(across(where(is.numeric), \(x) x < 0.05)) %>% print
normality_table %>% reframe(across(where(is.numeric), \(x) x < 0.05)) %>% colSums %>% print
## germination has the least normal distribution(s), which is expected since it's a 0-1 decimal that is ideally near 1, not really continuously distributed trait

#######
# test for equality of variance between groups before anova
# w Levene test
for(i in collist){
  print(i)
  leveneTest(get(i) ~ as.factor(Generation), joined_happops) %>% print
}
## all traits have (at least one) sig diff in variance btwn generations




#####################
#### signif changes in generation mean
joined_happops$Generation <- as.factor(joined_happops$Generation)
levels(joined_happops$Generation)

######
# is there a difference between using anova & kruskal test on these phenotypes?
for(m in 1:6) {
  print(colnames(joined_happops)[m])
  kruskal.test(unlist(joined_happops[,m]) ~ Generation, joined_happops) %>% print

  aov(unlist(joined_happops[,m]) ~ as.factor(Generation), joined_happops) %>% summary %>% print
}

#####
## do phenotypes signif.ly change from parents to F18?
## do phenotypes signif.ly change from F18 compared to F58?

# set up storing results
trait_list <- colnames(joined_happops)[1:6]

signif_tab <- tibble("trait"=trait_list, "P_F18"=as.numeric(NA), "F18_F58"=as.numeric(NA)) 

## first, test for significant result w Kruskal-wallis test

for(i in 1:length(signif_tab$trait)){
  # filter to only relevant generations
  tmp1 <- joined_happops %>% select(c(all_of(i), Generation)) %>% filter(Generation == 18 | Generation == 0)
  tmp2 <- joined_happops %>% select(c(all_of(i), Generation)) %>% filter(Generation == 18 | Generation == 58)

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
    test_trt_joined_happops <- joined_happops %>% select(c(all_of(test_trt), Generation))

    rownum <- which(signif_tab$trait == test_trt)

    dunres <- dunn.test(unlist(test_trt_joined_happops[,1]), test_trt_joined_happops$Generation)

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



