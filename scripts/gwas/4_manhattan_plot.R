#!/usr/bin/env Rscript
#SBATCH --job-name=gwas
#SBATCH -o /rhome/jmarz001/bigdata/Ag-Competition/scripts/gwas/4_manhattan_plot.stdout
#SBATCH --mem=40G
#SBATCH -t 02:00:00
#SBATCH -p koeniglab

# Manhattan plots of GEMMA results
library(pacman)
p_load(tidyverse, data.table, gridExtra, ggsci, Cairo)
options(stringsAsFactors = F)

setwd("/rhome/jmarz001/bigdata/Ag-Competition/results/gwas")

# load file with trait names and combine
phenotype_names <- read_delim("trait_name_to_col_numbers.tsv", col_names=T)


addPlot <- function(FileName){
  # look up file name in df that connects trait name and file
  AssocTraitName <- as.character(phenotype_names[grep(FileName, phenotype_names$file), 1])
  # change periods in trait name to spaces so it can be split by strwrap
  AssocTraitName <- gsub("\\.", " ", AssocTraitName)
  #AssocTraitName <- str_to_upper(AssocTraitName)
  # put line breaks into long trait names so they're readable
  AssocTraitName <- paste(strwrap(AssocTraitName, width=60), collapse="\n")

  # import data
  df <- fread(FileName)

  # Bonferroni threshold
  threshold <- 0.05/nrow(df)

  # parse locus
  names(df)[1] <- "CHR"
  df$CHR <- gsub("(chr\\w+)_\\w+_\\d+", "\\1", df$CHR)

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

    result <- result %>% filter(-log10(p_lrt)>2)
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

# set variables
#df <- read.table(args[1])
lst <- phenotype_names$file
test <- lapply(lst, addPlot)

pdf("GWAS_manhattan.pdf")
marrangeGrob(grobs=test, nrow=3, ncol=1)
dev.off()

lst_lmm <- phenotype_names$file_lmm
test_lmm <- lapply(lst_lmm, addPlot)

pdf("GWAS_manhattan_lmm.pdf")
marrangeGrob(grobs=test_lmm, nrow=3, ncol=1)
dev.off()