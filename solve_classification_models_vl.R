# Mehmet Gonen (mehmetgonen@ku.edu.tr)

library(Rcplex)

solve_classification_svm <- function(K, y, C, epsilon, scale_free = FALSE) {
  N <- length(y)
  yyK <- (y %*% t(y)) * K
  if(scale_free)
  {
    yyK <- C*yyK #Now yyK is multiplied by C
  }
  
  opts <- list()
  opts$trace <- 0
  opts$maxcalls <- 41 * 200
  
  start_time <- Sys.time()
  result <- Rcplex(cvec = rep(-1, N), 
                   Amat = matrix(y, nrow = 1, byrow = TRUE), 
                   bvec = 0, 
                   Qmat = yyK, 
                   lb = rep(0, N),
                   ub = rep(C, N),
                   control = opts,
                   objsense = "min",
                   sense = "E")
  end_time <- Sys.time()
  elapsed_time <- as.numeric(end_time-start_time, units="secs")
  
  alpha <- result$xopt[1:N]
  alpha_original <- alpha
  if(scale_free)
  {
    alpha[alpha < epsilon] <- 0
    alpha[alpha > (1 - epsilon)] <- +1
  }else
  {
    alpha[alpha < +C * epsilon] <- 0
    alpha[alpha > +C * (1 - epsilon)] <- +C 
  }
  
  objective <- sum(alpha) - 0.5 * (t(alpha) %*% yyK) %*% (alpha)
  objective <- objective * (objective >= 0)
  
  objective_original <- sum(alpha_original) - 0.5 * (t(alpha_original) %*% yyK) %*% (alpha_original)
  objective_original <- objective_original * (objective_original >= 0)
  
  if(scale_free)
  {
    objective <- C*objective #the overall objective is a multiple of C
    objective_original <- C*objective_original
  }
  
  support_indices <- which(alpha != 0)
  active_indices <- which(alpha != 0 & alpha < C)
  if(scale_free)
  {
    active_indices <- which(alpha != 0 & alpha < 1) #Now alpha <=1
  }
  
  if (length(active_indices) > 0) {
    b <- mean(y[active_indices] * (1 - yyK[active_indices, support_indices] %*% alpha[support_indices]))
  } else {
    b <- 0
  }
  
  model <- list(alpha = alpha * y, b = b, objective = objective, alpha_original = alpha_original*y, objective_original = objective_original, CPU = elapsed_time)
  
  return(model)
}

solve_vl_master_problem <- function(obj_coef, constraints_matrix, rhs, senses, lb, ub, vtype, TN, P, K, variable_indexer) {
  
  opts <- list()
  opts$trace <- 0
  opts$maxcalls <- 41 * 200
  
  start_time <- Sys.time()
  result <- Rcplex(cvec = obj_coef, 
                   Amat = constraints_matrix, 
                   bvec = rhs, 
                   lb = lb,
                   ub = ub,
                   vtype = vtype,
                   control = opts,
                   objsense = "min",
                   sense = senses)
  end_time <- Sys.time()
  elapsed_time <- as.numeric(end_time-start_time, units="secs")
  
  sln <- result$xopt
  Gamma <- numeric(TN)
  for(s in 1:TN)
    Gamma[s] <- sln[variable_indexer$gamma(s)]
  
  eta_l <- vector("list", K)
  for(l in 1:K)
  {
    et <- numeric(P)
    for(m in 1:P)
      et[m] <- sln[variable_indexer$eta_lm(l,m)]
    eta_l[[l]] <- et
  }
  
  z_sl <- vector("list", TN)
  for(s in 1:TN)
  {
    z_s <- numeric(K)
    for(l in 1:K)
      z_s[l] <- sln[variable_indexer$z_sl(s,l)]
    z_sl[[s]] <- z_s
  }
  
  objective <- result$obj
  
  output <- list(Gamma = Gamma, eta_l = eta_l, z_sl = z_sl, objective = objective, CPU = elapsed_time)
  return(output)
}