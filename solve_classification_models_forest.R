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

solve_mtmklc_forest_master_problem <- function(obj_coef, constraints_matrix, rhs, senses, lb, ub, 
                                            vtype, TN, P, variable_indexer, LP_relaxation = FALSE) {
  opts <- list()
  opts$trace <- 0
  opts$maxcalls <- 41 * 200
  
  if(LP_relaxation)
    vtype <- rep("C", length(vtype))
  
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
  eta <- list()
  for(s in 1:TN)
  {
    Gamma[s] <- sln[variable_indexer$gamma_s(s)]
    
    et <- numeric(P)
    for(m in 1:P)
      et[m] <- sln[variable_indexer$eta_sm(s,m)]
    eta[[s]] <- et
  }
  
  u_sq <- matrix(0, TN, TN)
  for(s in 1:(TN-1))
  {
    for(q in (s+1):TN)
    {
      u_sq[s,q] <- sln[variable_indexer$u_sq(s,q)]
    }
  }
  
  objective <- sum(sln*obj_coef)
  
  output <- list(Gamma = Gamma, eta = eta, u_sq = u_sq, objective = objective, CPU = elapsed_time)
  
  if(LP_relaxation)
  {
    dual_values =  result$extra$lambda
    rc_sq <- matrix(0, TN, TN)
    for(s in 1:(TN-1))
    {
      for(q in (s+1):TN)
      {
        i <- variable_indexer$u_sq(s,q)
        rc_sq[s,q] <- obj_coef[i] - sum(constraints_matrix[,i]*dual_values)
      }
    }
    output$dual_info = list(dual_values = dual_values, rc_sq = rc_sq)
  }
  
  output$result <- result
  output$input = list(obj_coef = obj_coef,
                      constraints_matrix = constraints_matrix,
                      rhs = rhs, 
                      senses = senses, lb = lb, ub = ub, vtype = vtype, TN = TN, P = P, 
                      variable_indexer = variable_indexer, LP_relaxation = LP_relaxation)
  return(output)
}

solve_mtmklc_bdf_master_problem <- function(obj_coef, constraints_matrix, rhs, senses, lb, ub, 
                                               vtype, TN, variable_indexer, LP_relaxation = FALSE) {
  
  opts <- list()
  opts$trace <- 0
  opts$maxcalls <- 41 * 200
  
  if(LP_relaxation)
    vtype <- rep("C", length(vtype))
  
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
  
  Mu <- sln[variable_indexer$mu()]
  
  u_sq <- matrix(0, TN, TN)
  for(s in 1:(TN-1))
  {
    for(q in (s+1):TN)
    {
      u_sq[s,q] <- sln[variable_indexer$u_sq(s,q)]
    }
  }
  
  # print(sln)
  
  objective <- sum(sln*obj_coef)
  
  output <- list(Mu = Mu, u_sq = u_sq, objective = objective, CPU = elapsed_time)
  
  if(LP_relaxation)
  {
    dual_values =  result$extra$lambda
    rc_sq <- matrix(0, TN, TN)
    for(s in 1:(TN-1))
    {
      for(q in (s+1):TN)
      {
        i <- variable_indexer$u_sq(s,q)
        rc_sq[s,q] <- obj_coef[i] - sum(constraints_matrix[,i]*dual_values)
      }
    }
    output$dual_info = list(dual_values = dual_values, rc_sq = rc_sq)
  }
  
  output$result <- result
  output$input = list(obj_coef = obj_coef,
                      constraints_matrix = constraints_matrix,
                      rhs = rhs, 
                      senses = senses, lb = lb, ub = ub, vtype = vtype, TN = TN, 
                      variable_indexer = variable_indexer, LP_relaxation = LP_relaxation)
  
  return(output)
}

solve_mtmklc_bdf_DSP_problem <- function(obj_coef, constraints_matrix, rhs, senses, lb, ub, 
                                            vtype, TN, P, alphas, variable_indexer, row_indexer, 
                                            initial_sol, strengthend_with_min_grouping = FALSE) {
  input <- as.list(environment())
  
  opts <- list()
  opts$trace <- 0
  opts$maxcalls <- 41 * 200
  
  # if(length(initial_sol) > 0) #warm start is not implemented in rcplex package
  #   problem$sol <- list(itr = list(xx = initial_sol))
  #   opts$usesol <- TRUE
  # }
  
  start_time <- Sys.time()
  result <- Rcplex(cvec = obj_coef, 
                   Amat = constraints_matrix, 
                   bvec = rhs, 
                   lb = lb,
                   ub = ub,
                   vtype = vtype,
                   control = opts,
                   objsense = "max",
                   sense = senses)
  end_time <- Sys.time()
  elapsed_time <- as.numeric(end_time-start_time, units="secs")
  
  sln <- result$xopt
  dual_values =  result$extra$lambda
  
  objective <- sum(sln*obj_coef)
  
  Theta <- numeric(TN)
  for(s in 1:TN)
    Theta[s] <- sln[variable_indexer$theta_s(s)]
  
  Beta <- vector("list", TN)
  for(s in 1:TN)
  {
    bt <- numeric(alphas[s])
    if(alphas[s] > 0)
    {
      for(a in 1:alphas[s])
      {
        vi <- variable_indexer$beta_sa(s,a)
        if(vi > 0)
          bt[a] <- sln[vi]
      }
    }
    Beta[[s]] <- bt
  }
  
  dual_solution <- list(Sol = sln, Theta = Theta, Beta = Beta)
  
  Gamma = numeric(TN)
  Eta <- vector("list", TN)
  for(s in 1:TN)
  {
    Gamma[s] <- dual_values[row_indexer$gamma_s(s)]
    
    eta <- numeric(P)
    for(m in 1:P)
    {
      eta[m] <- dual_values[row_indexer$eta_sm(s,m)]
    }
    Eta[[s]] <- eta
  }
  primal_solution <- list(Gamma = Gamma, Eta = Eta)
  
  if(strengthend_with_min_grouping)
  {
    Lambda <- matrix(0, TN, TN)
    for(s in 1:(TN-1))
    {
      for(q in (s+1):TN)
      {
        vi <- variable_indexer$lambda_sq(s,q)
        Lambda[s,q] <- sln[vi]
      }
    }
    dual_solution$Lambda <- Lambda
    
    Nu <- vector("list", P)
    for(m in 1:P)
    {
      nu <- matrix(0, TN, TN)
      for(s in 1:TN)
      {
        for(q in 1:TN)
        {
          if(s != q)
          {
            vi <- variable_indexer$nu_sqm(s,q,m)
            nu[s,q] <- sln[vi] 
          }
        }
      }
      Nu[[m]] <- nu
    }
    dual_solution$Nu <- Nu
    
    Phi <- vector("list", P)
    for(m in 1:P)
    {
      phi <- matrix(0, TN, TN)
      for(s in 1:(TN-1))
      {
        for(q in (s+1):TN)
        {
          phi[s,q] <- dual_values[row_indexer$phi_sqm(s,q,m)] 
        }
      }
      Phi[[m]] <- phi
    }
    primal_solution$Phi <- Phi
    
  }else
  {
    Lambda <- vector("list", P)
    for(m in 1:P)
    {
      lmbda <- matrix(0, TN, TN)
      for(s in 1:TN)
      {
        for(q in 1:TN)
        {
          vi <- variable_indexer$lambda_sqm(s,q,m)
          if(vi > 0) #otherwise it's zero
            lmbda[s,q] <- sln[vi]
        }
      }
      Lambda[[m]] <- lmbda
    }
    dual_solution$Lambda <- Lambda
  }
  
  output <- list(dual_solution = dual_solution, primal_solution = primal_solution, 
                 objective = objective, CPU = elapsed_time)
  
  output$input <- input
  
  return(output)
}