#!/usr/bin/env Rscript
#SBATCH --mem=50G
#SBATCH --time=02:00:00
#SBATCH --output=/rhome/jmarz001/bigdata/Ag-Competition/scripts/000_row_and_col_boxplots.stdout
#SBATCH -p short

setwd("/rhome/jmarz001/bigdata/Ag-Competition/data")
library(tidyverse)

# load FT data & tidy
ft_22 <- read_delim("FT_2021_2022.tsv")
  ft22 <- ft_22 %>% select(-c(Germination, Generation))
  colnames(ft22) <- c("Genotype", "Condition", "Replicate", "Bed", "Row", "FT")
  ft22$Exp_year <- 2022
ft_23 <- read_delim("FT_2023.tsv")
  ft23 <- ft_23 %>% select(-c(Plot_germination, PLOT_ID, Generation))
  colnames(ft23) <- c("Genotype", "Condition", "Replicate", "Bed", "Row", "FT")
  ft23$Exp_year <- 2023
## join ft year data frames
ft <- full_join(ft22, ft23) %>% filter(!is.na(FT))



## Flowering Time Traits First
trait = c("Flowering Time", "100 Seed Weight")
Replicate = c("Rep 1", "Rep 2")
Year = c(2022, 2023)

## plot flowering time for both years & Replicates
tt = trait[1]
#print(tt)

# add plots to a pdf
#pdf("Competition_row_and_col_boxplots.pdf", width=14, height=7)

stat_box_data <- function(y, upper_limit = max(tmp$FT) * 1.1) {
  return(
    data.frame(
      y = 0.95 * upper_limit,
      label = paste("n =", length(y), "\n",
                    "avg =", round(mean(y)), "\n",
                    "med =", round(median(y)), "\n")
    )
  )
}
stat_box_data_row <- function(y, upper_limit = max(tmp$FT) * 1.1) {
  return(
    data.frame(
      y = 0.95 * upper_limit,
      label = length(y)
    )
  )
}

for(yr in Year){
  #print(yr)
  for(rp in Replicate){
    #print(rp)

    # filter data to the right year & Replicate
    tmp <- ft %>% filter(Exp_year == yr) %>% filter(str_to_title(Replicate) == rp)
    print(head(tmp))

    ## title plots with rep # / year # / trait
    title = paste(tt, rp, yr, sep=" - ")
    print(title)


    b <- ggplot(tmp, aes(as.factor(Bed), y=FT)) +
          geom_boxplot() +
          stat_summary(
              fun.data = stat_box_data,
              geom = "text",
              hjust = 0.5,
              vjust = 0.9
            ) +
          labs(title= title, x = "Bed", y=tt) +
          theme_minimal()
    print(b)

    save_name_bed <- gsub(" ", "_", paste0(paste(tt, rp, yr, "BED", sep="_"), ".png"))
    ggsave(plot = b, filename = save_name_bed, width=13.4, height= 6.9, units = "in")



    b <- ggplot(tmp, aes(x=as.factor(Row), y=FT, group=as.factor(Row))) +
          geom_boxplot() +
          stat_summary(
              fun.data = stat_box_data_row,
              geom = "text",
              hjust = 0.5,
              vjust = 0.9
            ) +
          labs(title= title, x = "Row", y=tt) +
          theme_minimal()
    print(b)

    save_name_row <- gsub(" ", "_", paste0(paste(tt, rp, yr, "ROW", sep="_"), ".png"))
    ggsave(plot = b, filename = save_name_row, width=13.4, height= 6.9, units = "in")

  }
}

dev.off()
########stat_summary(fun = count, geom = "text")
#########stat_summary(fun.data=tmp$Row, fun=count, geom="text")

### Load seed weight phenotypes
wt <- read_delim("Ag-Comp Seed Weights - Sheet3.csv")
## tidy data
  wt <- wt %>% select(c(Genotype, Condition, Replicate, `BED_2021`, `ROW_2021`, FT, SEED_WEIGHT_100)) %>% filter(SEED_WEIGHT_100 > 1) %>% filter(SEED_WEIGHT_100 < 10)
  colnames(wt) <- c("Genotype", "Condition", "Replicate", "Bed", "Row", "FT", "SEED_WEIGHT_100")



  stat_box_data <- function(y, upper_limit = max(tmp$SEED_WEIGHT_100) * 1.1) {
    return(
      data.frame(
        y = 0.95 * upper_limit,
        label = paste("n =", length(y), "\n",
                      "avg =", round(mean(y)), "\n",
                      "med =", round(median(y)), "\n")
      )
    )
  }
  stat_box_data_row <- function(y, upper_limit = max(tmp$SEED_WEIGHT_100) * 1.1) {
    return(
      data.frame(
        y = 0.95 * upper_limit,
        label = length(y)
      )
    )
  }


yr=Year[1]
print(yr)
tt = trait[2]
print(tt)

## plot seed weight for 2022 Replicates
for(rp in Replicate){
    #print(rp)
    title = paste(tt, rp, yr, sep=" - ")

    # filter data to the right year & Replicate
    tmp <- wt %>% filter(str_to_title(Replicate) == rp)
    print(head(tmp))

    ## title plots with rep # / year # / trait
    title = paste(tt, rp, yr, sep=" - ")
    print(title)

    b <- ggplot(tmp, aes(as.factor(Bed), y=SEED_WEIGHT_100)) +
          geom_boxplot() +
          stat_summary(
              fun.data = stat_box_data,
              geom = "text",
              hjust = 0.5,
              vjust = 0.9
            ) +
          labs(title= title, x = "Bed", y=tt) +
          theme_minimal()
    print(b)

    save_name_bed <- gsub(" ", "_", paste0(paste(tt, rp, yr, "BED", sep="_"), ".png"))
    ggsave(plot = b, filename = save_name_bed, width=13.4, height= 6.9, units = "in")

    b <- ggplot(tmp, aes(x=as.factor(Row), y=SEED_WEIGHT_100, group=as.factor(Row))) +
          geom_boxplot() +
          stat_summary(
              fun.data = stat_box_data_row,
              geom = "text",
              hjust = 0.5,
              vjust = 0.9
            ) +
          labs(title= title, x = "Row", y=tt) +
          theme_minimal()
    print(b)

    save_name_row <- gsub(" ", "_", paste0(paste(tt, rp, yr, "ROW", sep="_"), ".png"))
    ggsave(plot = b, filename = save_name_row, width=13.4, height= 6.9, units = "in")
}

dev.off()
