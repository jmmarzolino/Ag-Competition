#!/usr/bin/env Rscript
#SBATCH --mem=30G
#SBATCH --time=02:00:00
#SBATCH --job-name='posthoc-HOC'
#SBATCH --output=/rhome/jmarz001/bigdata/Ag-Competition/scripts/004b_anova_posthoc_hapreppop.stdout
#SBATCH -p koeniglab

setwd("/bigdata/koeniglab/jmarz001/Ag-Competition")
library(tidyverse)
library(FSA)

anova_result <- read_delim("results/anova.tsv", "\t")
# subset trait list to those w significant anova result
anova_result <- anova_result %>% filter(unequal_means)
# divide traits based on un/equal variance
equ <- anova_result %>% filter(equal_vars) %>% select(trait)
equ <- equ[[1]]
uneq <- anova_result %>% filter(!equal_vars) %>% select(trait)
uneq <- uneq[[1]]
# divide data based on un/equal variance
df <- read_delim("data/trait_BLUPs_HapRepPop.tsv") 
df_equ <- df %>% select(c(Generation, all_of(equ)))
df_uneq <- df %>% select(c(Generation, all_of(uneq)))

##  anova posthoc test comparing generation 0 to 18, and 18 to 58
# record p-values for each trait tested
# traits from equal variance list
post_hoc_tukey <- function(data) {

  output <- data.frame()

  # loop over each trait column
  for(i in 2:ncol(data)) {

    # subset to trait of interest & save trait name
    tmp <- data[,c(1, i)]
    trait_storage <- colnames(tmp)[2]

    # save results of anova model
    fm1 <- aov(unlist(tmp[,2]) ~ as.factor(Generation), data = tmp)
    #summary(fm1)

    # run posthoc Tukey test
    x <- TukeyHSD(fm1)
    #plot(x)
    # save posthoc test result as data frame
    y <- data.frame(x$`as.factor(Generation)`)
    #z <- y["p.adj"]
    z <- rownames_to_column(y)
    z_sub <- z[c(1,7) , c(1,5)]

    #z_sub$post_hoc_testing_sig <- z_sub$`p.adj` < 0.05
    #z_sub <- z_sub %>% select(c(rowname, post_hoc_testing_sig))

    w <- pivot_wider(z_sub, values_from="p.adj", names_from="rowname", names_prefix="post_hoc_")
    colnames(w) <- gsub("-", "_", colnames(w))
    w$trait <- trait_storage

    output <- bind_rows(output, w)
    }

  return(output)
}

if(ncol(df_equ)>1){
    posthoc <- post_hoc_tukey(df_equ)
}


##  non-parametric anova posthoc test
# record p-values for each trait tested
# test traits with unequal variance
output <- data.frame()
for(i in 2:ncol(df_uneq)){
        trait_storage <- colnames(df_uneq)[i]
        result <- dunnTest(unlist(df_uneq[,i]) ~ as.factor(Generation),
                  data = df_uneq, method = "holm")

        result <- result$res %>% tibble %>% select(c(Comparison, P.adj))
        #result$posthoc_sig <- result$P.adj < 0.05
        result <- result %>% select(c(Comparison, P.adj)) %>% filter(Comparison == c("0 - 18", "18 - 58"))

        result$Comparison <- gsub("(\\d+) - (\\d+)", "\\2_\\1", result$Comparison)

        m <- pivot_wider(result, values_from="P.adj", names_from="Comparison", names_prefix="post_hoc_")
        m$trait <- trait_storage

        output <- bind_rows(output, m)
}

# combine results of tukey and dunn tests

if(exists(quote(posthoc))){
    posthoc2 <- full_join(output, posthoc)
    anova_posthoc <- full_join(anova_result, posthoc2, by='trait')
} else {
    anova_posthoc <- full_join(anova_result, output, by='trait')
}
write_delim(anova_posthoc, "results/anova_posthoc.tsv", delim="\t")
