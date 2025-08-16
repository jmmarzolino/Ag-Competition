#!/usr/bin/env Rscript
#SBATCH --job-name=gwas
#SBATCH -o /rhome/jmarz001/bigdata/Ag-Competition/scripts/gwas/greenhouse-comparison/4_manhattan_plot_cross_exp_plots.stdout
#SBATCH --mem=40G
#SBATCH -t 02:00:00
#SBATCH -p koeniglab

# Manhattan plots of GEMMA results
library(pacman)
p_load(tidyverse, data.table, gridExtra, ggsci, Cairo)
options(stringsAsFactors = F)

setwd("/rhome/jmarz001/bigdata/Ag-Competition/results/gwas/CCII_greenhouse_exp_gwas")
source("../../../scripts/CUSTOM_FNS.R")

# load file with trait names and combine
phenotype_names <- read_delim("CCII_GH_trait_file_nums.tsv", col_names=T)
phenotype_names <- phenotype_names[which(phenotype_names$trait_names %in% c("Mass_100","Seed_Estimate", "Flowering_2018_Median", "Flowering_days_2017", "avg_FT")), ]

field_exp <- read_delim("../trait_name_to_col_numbers.tsv", col_names=T)
field_exp <- field_exp[which(field_exp$trait_names %in% c("SEED_WEIGHT_100", "FECUNDITY", "FT")), ]



addPlot <- function(FileName, TraitName, Experiment){
  # convert computer-style trait names to human-readable
  AssocTraitName <- tidy_text_substitution(TraitName)
  AssocTraitName <- paste0(AssocTraitName, " measured in ", Experiment)
  # put line breaks into long trait names so they're readable
  AssocTraitName <- paste(strwrap(AssocTraitName, width=60), collapse="\n")

  # import data
  df <- fread(FileName)

  # Bonferroni threshold
  threshold <- 0.05/nrow(df)

  # parse locus
  names(df)[1] <- "CHR"

  # following code adapted from:
    # https://www.r-graph-gallery.com/wp-content/uploads/2018/02/Manhattan_plot_in_R.html
    # format for plotting
    df$BP <- as.numeric(df$ps)

    result <- df %>%
      # Compute chromosome size
      group_by(CHR) %>%
      summarise(chr_len = max(BP)) %>%
      # Calculate cumulative position of each chromosome
      mutate(tot = cumsum(as.numeric(chr_len))-as.numeric(chr_len)) %>%
      select(-chr_len) %>%
      # Add this info to the initial dataset
      left_join(df, ., by=c("CHR" = "CHR")) %>%
      # Add a cumulative position of each SNP
      arrange(CHR, BP) %>%
      mutate(BPcum = BP+tot)

    #result <- result %>% filter(-log10(p_lrt)>2)
    axisdf <- result %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

    # Manhattan plot
    ggplot(result, aes(x=BPcum, y=-log10(p_lrt))) +
        # Show all points
        geom_point(aes(color=as.factor(CHR)), alpha=0.5, size=1.3) +
        geom_hline(aes(yintercept=-log10(threshold)), color = "firebrick1", linetype="dashed", alpha=0.7) +
        scale_color_manual(values = rep(c("dodgerblue4", "deepskyblue"), 22 )) +
        # custom X axis:
        scale_x_continuous(label = axisdf$CHR, breaks= axisdf$center) +
        scale_y_continuous(expand = c(0, 0.5)) +     # remove space between plot area and x axis
        # Custom theme:
        theme_classic() +
        theme(legend.position="none",
              panel.border = element_blank(),
              panel.grid.major.x = element_blank(),
              panel.grid.minor.x = element_blank(),
              text=element_text(size=16),
              plot.title = element_text(size=12),
              plot.subtitle = element_text(size=10)) +
        xlab("Chromosome") +
        ylab(expression(-log[10](italic(p)))) +
        ggtitle(AssocTraitName)

    #print(g)
}


# trait pairs
## FT - Flowering_days_2017
## FT - Flowering_2018_Median
## SEED_WEIGHT_100 - Mass_100
## FECUNDITY - Seed_Estimate


# set variables
# manually compile matched lists of traits & files
matched_exp_ft <- tibble(
          "experiment"=c("field", "greenhouse", "greenhouse", "greenhouse") ,
          "trait"=c("FT", "Flowering_days_2017", "Flowering_2018_Median", "avg_FT") , 
          "file"=c("../ASSOC_6_lmm.assoc.txt", "ASSOC_10_lmm.assoc.txt", "ASSOC_11_lmm.assoc.txt", "ASSOC_12_lmm.assoc.txt")
          )

matched_exp_fec <- tibble(
          "experiment"=c("field", "greenhouse") ,
          "trait"=c("FECUNDITY", "Seed_Estimate") , 
          "file"=c("../ASSOC_8_lmm.assoc.txt", "ASSOC_11_lmm.assoc.txt")
          )

matched_exp_sw <- tibble(
          "experiment"=c("field", "greenhouse") ,
          "trait"=c("SEED_WEIGHT_100", "Mass_100") , 
          "file"=c("../ASSOC_7_lmm.assoc.txt", "ASSOC_9_lmm.assoc.txt")
          )

combo_traits <- bind_rows(matched_exp_ft, matched_exp_sw, matched_exp_fec, )


pdf("FLD_GH_GWAS_manhattan.pdf")
for(i in 1:length(combo_traits$trait)){

  gggg <- addPlot(FileName = combo_traits[[i, 3]], TraitName=combo_traits[i, 2], Experiment=combo_traits[i, 1])
  print(gggg)

}
dev.off()

#lst <- phenotype_names$file
#test <- lapply(lst, addPlot)
#marrangeGrob(grobs=test, nrow=2, ncol=1)