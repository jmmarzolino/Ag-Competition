# add a generation column based on genotype UCRKL numbers

add_generation <- function(x){

  # parent IDs of 1 number get generation 0
  x$Generation <- gsub("^\\d+$", 0, x$Genotype)
  x$Generation <- gsub("^1_.*", 18, x$Generation)
  x$Generation <- gsub("^2_.*", 28, x$Generation)
  x$Generation <- gsub("^3_.*", 50, x$Generation)
  x$Generation <- gsub("^7_.*", 58, x$Generation)
  # any lines left with an _ should be parents of generation 0
  x$Generation <- gsub("^*.*_.*", 0, x$Generation)

  x$Generation <- as.numeric(x$Generation)
  return(x)
}

# substitute trait names w/ tidy text versions
tidy_text_substitution <- function(x) {

  x <- gsub("_blup", " BLUP", x)
  x <- gsub("FT", "Flowering Time", x)
  x <- gsub("TOTAL_MASS", "Total Seed Mass", x)
  x <- gsub("GERMINATION", "Germination", x)
  x <- gsub("Plants", "Germination", x)
  x <- gsub("SEED_WEIGHT_100", "100-Seed Weight", x)
  x <- gsub("FECUNDITY", "Fecundity", x)
  x <- gsub("FITNESS", "Fitness", x)
  x <- gsub("MASS_PER_PLANT", "Mass per Plant", x)
  
  x <- gsub("X100_seed_mass", "greenhouse 100-Seed Weight", x)
  x <- gsub("seed_estimate", "greenhouse Fecundity", x)
  x <- gsub("days_to_heading", "greenhouse Flowering Time", x)
  return(x)
  }


Exp_Single <- function(x){
  result_single <- x/10
  return(result_single)
}

Exp_Mixed <- function(x, y){
  # note - x and y must be table columns or vectors
  Exp_mix <- ((x/2) + (y/2))/10
  return(Exp_mix)
}



#Expected_seed <- function(x, y){
#  separate the rows
#  if(grep("single")){
#    aapoigapoiefd
#  }
#  if(grep("mixed")){}
#    else("print error")
#}
