# Mehmet Gonen (mehmetgonen@ku.edu.tr)
pdist <- function(X1, X2) {
  if (identical(X1, X2) == TRUE) {
    D <- as.matrix(dist(X1))
  }
  else {
    D <- as.matrix(dist(rbind(X1, X2)))
    D <- D[1:nrow(X1), (nrow(X1) + 1):(nrow(X1) + nrow(X2))]
  }
  return(D)
}

read_pathways <- function(name) {
  symbols_lines <- read.table(sprintf("%s.gmt", name), header = FALSE, sep = ",", stringsAsFactor = FALSE)
  pathways <- vector("list", nrow(symbols_lines))
  for (line in 1:nrow(symbols_lines)) {
    symbols_entries <- strsplit(symbols_lines[line, 1], "\t")
    pathways[[line]]$name <- symbols_entries[[1]][1]
    pathways[[line]]$link <- symbols_entries[[1]][2]
    pathways[[line]]$symbols <- sort(symbols_entries[[1]][-2:-1])
  }
  return(pathways)
}

calculate_Keta <- function(Km, eta) {
  P <- dim(Km)[3]
  
  Keta <- eta[1] * Km[,,1]
  for (m in 2:P) {
    Keta <- Keta + eta[m] * Km[,,m]
  }
  
  return(Keta)
}

calucalte_kappa <- function(alpha, Km, P)
{
  kappa <- numeric(P)
  
  for(m in 1:P)
  {
    kappa[m] <- 0.5*t(alpha) %*% Km[,,m] %*% alpha
  }
  return(kappa)
}

normalize_eta <- function(eta, epsilon)
{
  et <- eta / sum(eta)
  et[et < epsilon] <- 0
  et <- et / sum(et)
  return(et)
}

get_cross_validation_tuples <- function()
{
  C_set <- 10^(0:6)
  
  tuples <- matrix(C_set, ncol = 1, nrow = length(C_set))
  colnames(tuples) <- c("C")
  
  return(tuples)
}

get_core_point(TN, K, phi)
{
  core_point <- matrix(0, TN, TN)
  for(s in 1:(TN-1))
    for(q in (s+1):TN)
      core_point[s,q] <- phi^(q-s)
  core_point <- core_point/sum(score_point)*(TN-K)
  
  return(core_point)
}
