#!/usr/bin/env Rscript
#SBATCH --ntasks=1
#SBATCH --mem=30G
#SBATCH --time=02:00:00
#SBATCH --output=/rhome/jmarz001/bigdata/Ag-Competition/scripts/plotting/2_blup_correlations.stdout
#SBATCH -p koeniglab

library(tidyverse)
library(ggpubr)
library(RColorBrewer)
library(data.table)
library(corrplot)

setwd("/rhome/jmarz001/bigdata/Ag-Competition/results")
source("../scripts/CUSTOM_FNS.R")

df <- fread("../data/BLUPs.tsv")
df_gen <- add_generation(df) %>% select(-Genotype)
df <- df %>% select(-c(Genotype))


# save table of trait correlation values
x_cor <- cor(df, use="na.or.complete", method="spearman")
dimnames(x_cor)[[1]] <- tidy_text_substitution(dimnames(x_cor)[[1]])
dimnames(x_cor)[[2]] <- tidy_text_substitution(dimnames(x_cor)[[2]])
x <- data.table(x_cor)
x$TRAIT <- colnames(x)
write_delim(tibble(x), "correlation_table.tsv", "\t")


# test minimal trait correlations with spearmans rho
# and output p-values
cor.test(df$FT, df$SEED_WEIGHT_100, method="spearman") %>% print
#getOption("na.action")
#na.action = c('no.omit', 'na.fail', 'na.exclude', 'na.pass')
cor.test(df$FT, df$FECUNDITY, method="spearman") %>% print
cor.test(df$FECUNDITY, df$SEED_WEIGHT_100, method="spearman") %>% print


# plot trait relationships 
png("blup_correlations.png", width=5, height=5, units="in", res=300)
corrplot(x_cor, method="color", type="upper", order="original", title="", mar=c(0,0,4,0), addCoef.col = "black", tl.col = "black")
dev.off()


# plot trait relationships minus flowering time which might skew corr values
tmp <- df %>% select(-FT)
tmp2 <- cor(tmp, use="na.or.complete", method="spearman")

png("blup_correlations_minusFT.png", height=10, width=10, units="in", res=200)
corrplot(tmp2, method="color", type="upper", order="original", title="", mar=c(0,0,4,0), addCoef.col = "black")
dev.off()



########## plot trait correlations for each generation in turn
pdf("blup_correlations_per_generation.pdf")
# check correlations between traits for each generation
for(i in c(0, 18, 28, 50, 58)) {

  traits_df <- df_gen %>% filter(Generation == i) %>% select(-c(Generation)) 
  x <- cor(traits_df, use="na.or.complete", method="spearman")
  corrplot(x, method="color", type="upper", order="original", title=paste0("Generation ", i), mar=c(0,0,4,0), addCoef.col = "black", tl.col = "black")

}
dev.off()



########### plot trait relationships w scatterplots
df_scat <- df
colnames(df_scat) <- tidy_text_substitution(colnames(df_scat))
png("blup_correlations_scatterplot.png", height=10, width=10, units="in", res=200)
plot(df_scat, main="Parents & Progeny")
dev.off()

# plot w progeny only
progeny <- df_gen %>% filter(Generation != 0) %>% select(-c(Generation))
colnames(progeny) <- tidy_text_substitution(colnames(progeny))

png("blup_correlations_scatterplot_progeny.png", height=10, width=10, units="in", res=200)
plot(progeny, main="Progeny Only")
dev.off()




########## trait relationships over generations
# to see more clearly if the relationship between traits are changing over generations
# plot trait v. trait, colored & approx line by generation
# and look to see if value range or relationship change much

# the lowest color value is too light, so adjust the color scale down one
#display.brewer.pal(6, "Blues")
adjusted_blues <- brewer.pal(7, "Blues")[3:7]

# flowering time plots
df_ft <- df_gen %>% 
        pivot_longer(-c(FT, Generation), values_to = "VALUE", names_to="TRAIT")
df_ft$TRAIT <- tidy_text_substitution(df_ft$TRAIT)


gee <- df_ft %>% 
        ggplot(aes(x=FT, y=VALUE,)) + 
          geom_point(alpha=0.7) + 
          geom_smooth(method="lm", se=T) +
          scale_color_manual(values=adjusted_blues) + 
          facet_wrap(~TRAIT, scales="free_y") +
          labs(x="Flowering Time BLUP", y="", color="Generation") +
          theme_bw()
ggsave("blups_vs_FT.png", gee, width=(7*3), height=(7*2))


gee <- df_ft %>% 
        ggplot(aes(x=FT, y=VALUE,)) + 
          geom_point(aes(color=as.factor(Generation)), alpha=0.7) + 
          geom_smooth(method="lm", se=T) +
          scale_color_manual(values=adjusted_blues) + 
          facet_wrap(~TRAIT, scales="free_y") +
          labs(x="Flowering Time BLUP", y="", color="Generation") +
          theme_bw()
ggsave("blups_vs_FT_color_gen.png", gee, width=(7*3), height=(7*2))

gee <- df_ft %>% 
        ggplot(aes(x=FT, y=VALUE, group=as.factor(Generation))) + 
          geom_point(aes(color=as.factor(Generation)), alpha=0.7) + 
          geom_smooth(method="lm", se=T, aes(color=as.factor(Generation))) +
          scale_color_manual(values=adjusted_blues) + 
          facet_wrap(~TRAIT, scales="free_y") +
          labs(x="Flowering Time BLUP", y="", color="Generation") +
          theme_bw()
ggsave("blups_vs_FT_by_generation.png", gee, width=(7*3), height=(7*2))


# seed-weight plots
df_sw <- df_gen %>% 
        pivot_longer(-c(SEED_WEIGHT_100, Generation), values_to = "VALUE", names_to="TRAIT")
df_sw$TRAIT <- tidy_text_substitution(df_sw$TRAIT)


golly <- df_sw %>% 
        ggplot(aes(x=SEED_WEIGHT_100, y=VALUE)) + 
          geom_point(alpha=0.7) + 
          geom_smooth(method="lm", se=T) +
          scale_color_manual(values=adjusted_blues) + 
          facet_wrap(~TRAIT, scales="free_y") +
          labs(x="100-Seed Weight BLUP", y="", color="Generation") +
          theme_bw()
ggsave("blups_vs_100sw.png", golly, width=(7*3), height=(7*2))

golly <- df_sw %>% 
        ggplot(aes(x=SEED_WEIGHT_100, y=VALUE)) + 
          geom_point(aes(color=as.factor(Generation)), alpha=0.7) + 
          geom_smooth(method="lm", se=T) +
          scale_color_manual(values=adjusted_blues) + 
          facet_wrap(~TRAIT, scales="free_y") +
          labs(x="100-Seed Weight BLUP", y="", color="Generation") +
          theme_bw()
ggsave("blups_vs_100sw_color_gen.png", golly, width=(7*3), height=(7*2))

golly <- df_sw %>% 
        ggplot(aes(x=SEED_WEIGHT_100, y=VALUE, group=Generation)) + 
          geom_point(aes(color=as.factor(Generation)), alpha=0.7) + 
          geom_smooth(method="lm", se=T, aes(color=as.factor(Generation)), alpha=0.5) +
          scale_color_manual(values=adjusted_blues) + 
          facet_wrap(~TRAIT, scales="free_y") +
          labs(x="100-Seed Weight BLUP", y="", color="Generation") +
          theme_bw()
ggsave("blups_vs_100sw_by_generation.png", golly, width=(7*3), height=(7*2))




# FECUNDITY plots
df_fec <- df_gen %>% 
        pivot_longer(-c(FECUNDITY, Generation), values_to = "VALUE", names_to="TRAIT")
df_fec$TRAIT <- tidy_text_substitution(df_fec$TRAIT)

gee_wiz <- df_fec %>% 
        ggplot(aes(x=FECUNDITY, y=VALUE)) + 
          geom_point(alpha=0.7) + 
          geom_smooth(method="lm", se=T) +
          scale_color_manual(values=adjusted_blues) + 
          facet_wrap(~TRAIT, scales="free_y") +
          labs(x="Mass per Plant", y="", color="Generation") +
          theme_bw()
ggsave("blups_vs_FEC.png", gee_wiz, width=(7*3), height=(7*2))

gee_wiz <- df_fec %>% 
        ggplot(aes(x=FECUNDITY, y=VALUE)) + 
          geom_point(aes(color=as.factor(Generation)), alpha=0.7) + 
          geom_smooth(method="lm", se=T) +
          scale_color_manual(values=adjusted_blues) + 
          facet_wrap(~TRAIT, scales="free_y") +
          labs(x="Mass per Plant", y="", color="Generation") +
          theme_bw()
ggsave("blups_vs_FEC_color_gen.png", gee_wiz, width=(7*3), height=(7*2))


gee_wiz <- df_fec %>% 
        ggplot(aes(x=FECUNDITY, y=VALUE, group=Generation)) + 
          geom_point(aes(color=as.factor(Generation)), alpha=0.7) + 
          geom_smooth(method="lm", se=T, aes(color=as.factor(Generation))) +
          scale_color_manual(values=adjusted_blues) + 
          facet_wrap(~TRAIT, scales="free_y") +
          labs(x="Fecundity", y="", color="Generation") +
          theme_bw()
ggsave("blups_vs_FEC_by_generation.png", gee_wiz, width=(7*3), height=(7*2))









###############################
# fit linear & quadratic regressions to trait relationships
# and report AIC values in table

#FT ~ FEC 
#FT ~ SW
#SW ~ FEC

# check AIC of FEC ~ FT vs MASSPER ~ FT^2
ft_fec_lm_aic <- print(lm(`FT` ~ `FECUNDITY`, data=traits_df) %>% AIC)
ft_fec_lmpoly_aic <- print(lm(`FT` ~ poly(`FECUNDITY`, 2), data=traits_df) %>% AIC)

ft_sw_lm_aic <- print(lm(`FT` ~ `SEED_WEIGHT_100`, data=traits_df) %>% AIC)
ft_sw_lmpoly_aic <- print(lm(`FT` ~ poly(`SEED_WEIGHT_100`, 2), data=traits_df) %>% AIC)

sw_fec_lm_aic <- print(lm(`SEED_WEIGHT_100` ~ `FECUNDITY`, data=traits_df) %>% AIC)
sw_fec_lmpoly_aic <- print(lm(`SEED_WEIGHT_100` ~ poly(`FECUNDITY`, 2), data=traits_df) %>% AIC)

aic_table <- tibble("regression"=c("linear", "quadratic"), "FT~FEC"=c(ft_fec_lm_aic, ft_fec_lmpoly_aic), "FT~SW"=c(ft_sw_lm_aic, ft_sw_lmpoly_aic), "SW~FEC"=c(sw_fec_lm_aic, sw_fec_lmpoly_aic))
fwrite(aic_table, "trait_relationship_regression_model_fit_AICs.tsv")


## try to fit a quadratic relationship on FT ~ FEC 
g1 <- ggplot(traits_df, aes(x=FT, y=FECUNDITY)) +
        geom_point() +
        geom_smooth(method="lm", formula = 'y ~ poly(x, 2)') +
        theme_bw(base_size = 18)
ggsave("ftxFEC_quadratic.png", g1)

g1 <- ggplot(traits_df, aes(x=FT, y=SEED_WEIGHT_100)) +
        geom_point() +
        geom_smooth(method="lm", formula = 'y ~ poly(x, 2)') +
        theme_bw(base_size = 18)
ggsave("ftxSW_quadratic.png", g1)

g1 <- ggplot(traits_df, aes(x=SEED_WEIGHT_100, y=FECUNDITY)) +
        geom_point() +
        geom_smooth(method="lm", formula = 'y ~ poly(x, 2)') +
        theme_bw(base_size = 18)
ggsave("SWxFEC_quadratic.png", g1)

# geom_smooth() method=
#‘"lm"’, ‘"glm"’, ‘"gam"’, ‘"loess"’
# or a function, e.g. ‘MASS::rlm’ or ‘mgcv::gam’, ‘stats::lm’,
