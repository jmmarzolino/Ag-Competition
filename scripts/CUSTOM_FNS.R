# add a generation column based on genotype UCRKL numbers

add_generation <- function(x){
  x$Generation <- gsub("^1_.*", 18, x$Genotypes)
  x$Generation <- gsub("^2_.*", 28, x$Generation)
  x$Generation <- gsub("^3_.*", 50, x$Generation)
  x$Generation <- gsub("^7_.*", 58, x$Generation)
  x$Generation <- gsub("^*.*_.*", 0, x$Generation)
  return(x)
}

Exp_Single <- function(x){
  result_single <- x/10
  return(result_single)
}

Exp_Fec_Mixed <- function(x){
  TW_mix <- (x/2) + (Averaged_Full_2021_2022$Atlas_Avg_Fec/2)
  Exp_Fec_mix <- TW_mix/10
  return(Exp_Fec_mix)
}

Exp_Fit_Mixed <- function(x){
  fit_mix <- (x/2) + (Averaged_Full_2021_2022$Atlas_Avg_Fit/2)
  Exp_Fit_mix <- fit_mix/10
  return(Exp_Fit_mix)
}

Exp_TW_mix <- function(x){
  TW_mix <- (x/2) + (Average_Haplo_rep$Atlas_Avg_Total_Weight/2)
  Exp_TW <- TW_mix/10
  return(Exp_TW)
}
