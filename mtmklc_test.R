# Mehmet Gonen (mehmetgonen@ku.edu.tr)

mtmklc_test <- function(Km, eta, alpha, b) {
  T <- length(Km)
  Keta <<- vector("list", T)
  for (t in 1:T) {
    Keta[[t]] <<- calculate_Keta(Km[[t]], eta[[t]])
  }
  
  f <- vector("list", T)
  for (t in 1:T) {
    f[[t]] <- Keta[[t]] %*% alpha[[t]] + b[[t]] 
  }
  
  prediction <- list(f = f)
}