FT_FITNESS <- read_sheet('https://docs.google.com/spreadsheets/d/15-7DX0YVGhldTwaW6nkKnNhryFmxEwo2ZHiZwAtBu58/edit#gid=1001803440')

### Create function to unlist and convert vectors to numeric
unlist_numeric <- function(x){
  unlist(x) %>% 
    as.numeric(x)
}

### Replace negative TW with NA, unlist all vectors and convert them to numeric, fill in missing generation values
FT_FITNESS <- FT_FITNESS %>% relocate(Exp_year, .after = replicate)
FT_FITNESS[5:14] <- lapply(FT_FITNESS[5:14], unlist_numeric)
FT_FITNESS[which(FT_FITNESS$total_seed_mass_g < 0), (7:ncol(FT_FITNESS))] <- NA

FT_FITNESS$Generation <- gsub("^1_.*", 18, FT_FITNESS$Genotypes)
FT_FITNESS$Generation <- gsub("^2_.*", 28, FT_FITNESS$Generation)
FT_FITNESS$Generation <- gsub("^3_.*", 50, FT_FITNESS$Generation)
FT_FITNESS$Generation <- gsub("^7_.*", 58, FT_FITNESS$Generation)
FT_FITNESS$Generation <- gsub("^*.*_.*", 0, FT_FITNESS$Generation)

outlier_cutoff = quantile(test$FECUNDITY,0.75, na.rm = TRUE) + (1.5 * IQR(test$FECUNDITY, na.rm = TRUE))
ggplot(FT_FITNESS, aes(x = FECUNDITY)) +
  geom_histogram(bins = 50) +
  geom_vline(xintercept = outlier_cutoff, color = 'red')
