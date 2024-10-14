# add a generation column based on genotype UCRKL numbers

add_generation <- function(x){
  x$Generation <- gsub("^1_.*", 18, x$Genotype)
  x$Generation <- gsub("^2_.*", 28, x$Generation)
  x$Generation <- gsub("^3_.*", 50, x$Generation)
  x$Generation <- gsub("^7_.*", 58, x$Generation)
  x$Generation <- gsub("^*.*_.*", 0, x$Generation)
  x$Generation <- as.numeric(x$Generation)
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
