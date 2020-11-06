# Mehmet Gonen (mehmetgonen@ku.edu.tr)
mtmklc_vl_train <- function(Km, y, parameters, settings) {
  # parameters should have: K, optimality_gap
  
  TN <- length(Km)
  P <- dim(Km[[1]])[3]
  K <- parameters$K
  
  terminate <- FALSE
  UB <- Inf
  LB <- -Inf
  
  best_eta <- NULL
  best_alpha <- NULL
  best_b <- NULL
  best_grouping <- NULL
  
  iteration <- 1
  
  #define constraint matrix here, add the equality constraint. Also define the objective function
  #order of variables: Gamma (TN), eta_l (K,P), eta_sl (TN,K,P), z_sl (TN,K)
  NVars <- TN + K*P+TN*K*P + TN*K
  
  variable_indexer <- list()
  variable_indexer$gamma <- function(s)
  {
    return(s)
  }
  variable_indexer$eta_lm <- function(l,m)
  {
    return(TN+(l-1)*P+m)
  }
  variable_indexer$eta_slm <- function(s,l,m)
  {
    return(TN+K*P+(s-1)*K*P+(l-1)*P+m)
  }
  variable_indexer$z_sl <- function(s,l)
  {
    return(TN+K*P+TN*K*P+(s-1)*K+l)
  }
  
  master_obj_coef <- numeric(NVars)
  master_obj_coef[1:TN] <- 1
  
  master_constraint_matrix <- matrix(0, nrow = ifelse(settings$force_nonempty_groups,  K+TN+TN*K*P+TN*K+K,  K+TN+TN*K*P+TN*K), ncol = NVars)
  
  #Constraint 3:
  row <- 0
  for(l in 1:K)
  {
    for(m in 1:P)
    {
      master_constraint_matrix[row+l,variable_indexer$eta_lm(l,m)] <- 1
    }
  }
  master_rhs <- rep(1,K)
  master_senses <- rep("E",K)
  
  #Constraint 4:
  row <- K
  for(s in 1:TN)
  {
    for(l in 1:K)
    {
      master_constraint_matrix[row+s,variable_indexer$z_sl(s,l)] <- 1
    }
  }
  master_rhs <- c(master_rhs, rep(1,TN))
  master_senses <- c(master_senses,rep("E",TN))
  
  #Constraint 13:
  row <- K+TN
  for(s in 1:TN)
  {
    for(l in 1:K)
    {
      for(m in 1:P)
      {
        row <- row + 1
        master_constraint_matrix[row,variable_indexer$eta_slm(s,l,m)] <- 1
        master_constraint_matrix[row,variable_indexer$eta_lm(l,m)] <- -1
      }
    }
  }
  master_rhs <- c(master_rhs, rep(0,TN*K*P))
  master_senses <- c(master_senses,rep("L",TN*K*P))
  
  #Constraint 15:
  row <- K + TN + TN*K*P
  for(s in 1:TN)
  {
    for(l in 1:K)
    {
      row <- row + 1
      for(m in 1:P)
      {
        master_constraint_matrix[row,variable_indexer$eta_slm(s,l,m)] <- 1
      }
      master_constraint_matrix[row,variable_indexer$z_sl(s,l)] <- -1
    }
  }
  master_rhs <- c(master_rhs, rep(0,TN*K))
  master_senses <- c(master_senses,rep("E",TN*K))
  
  if(settings$force_nonempty_groups)
  {
    row <- K + TN + TN*K*P + TN*K 
    for(l in 1:K)
    {
      for(s in 1:TN)
      {
        master_constraint_matrix[row+l,variable_indexer$z_sl(s,l)] <- 1
      }
    }
    master_rhs <- c(master_rhs,rep(1,K))
    master_senses <- c(master_senses,rep("G",K))
  }
  
  master_lb <- numeric(NVars)
  master_ub <- rep(Inf, NVars)
  
  master_vtype <- rep("C", NVars)
  for(s in 1:TN)
  {
    for(l in 1:K)
    {
      master_vtype[variable_indexer$z_sl(s,l)] <- "B"
      master_ub[variable_indexer$z_sl(s,l)] <- 1
    }
  }
  
  Gamma <- numeric(TN)
  J <- numeric(TN)
  alpha <- vector("list", TN)
  b <- vector("list", TN)
  
  eta_l <- vector("list", K)
  for (l in 1:K) {
    eta_l[[l]] <- rep(1 / P, P)
  }
  
  cluster_size <- ceiling(TN/K)
  z_sl <- list()
  for (s in 1:TN) {
    z_initial <- rep(0,K)
    z_initial[floor((s-1)/cluster_size)+1] <- 1
    z_sl[[s]] <- z_initial
  }
  
  objectives <- c()
  total_time <- 0
  
  info_cols <- c("obj", "UB", "LB", "Gap","Gamma","J","Time(M)", "Time(S)")
  iteration_info <- matrix(0,0,length(info_cols))
  colnames(iteration_info) <- info_cols
  
  min_iterations <- 6
  if(parameters$iteration_count > 0 & min_iterations > parameters$iteration_count)
    min_iterations <- parameters$iteration_count/2
  
  max_nonimproving_iterations <- 100
  non_improving_iterations <- 0
  latest_gap <- 1e+10
  
  while(terminate == FALSE)
  {
    grouping <- rep(0, TN)
    CPU_Master <- 0
    CPU_SVMs <- rep(0,TN)
    if(iteration > 1)
    {
      #SOLVE MP and obtain Gamma[s], eta_l[[l]], z_sl[[s]]
      master_result <- solve_model3P1_master_problem(master_obj_coef, master_constraint_matrix, master_rhs, master_senses, master_lb, master_ub, master_vtype, TN, P, K, variable_indexer)
      Gamma <- master_result$Gamma
      
      CPU_Master <- master_result$CPU
      
      eta_l <- master_result$eta_l
      
      for(s in 1:TN)
      {
        z_sl[[s]] <- master_result$z_sl[[s]]
      }
    }
    
    for(s in 1:TN)
    {
      grouping[s] <- which(z_sl[[s]] >= 0.99)
    }
    
    kappa <- vector("list", TN)
    
    for (s in 1:TN) {
      Keta <- calculate_Keta(Km[[s]], eta_l[[grouping[s]]])
      C <- parameters$C
      if(parameters$normalize_C)
        C <- C / length(y[[s]])
      svm_result <- solve_classification_svm(Keta, y[[s]], C, parameters$epsilon)
      
      CPU_SVMs[s] <- svm_result$CPU
      
      alpha[[s]] <- svm_result$alpha
      b[[s]] <- svm_result$b
      
      alpha_cut <- svm_result$alpha_original
      J[s] <- svm_result$objective_original
      
      #model$alpha is already multiplied by y
      alpha_sum <- t(alpha_cut) %*% y[[s]]
      kappa[[s]] <- calucalte_kappa(alpha_cut, Km[[s]], P)
      
      #add constraint Gamma_s + Sum_{lm}Eta_slm*kappa_sm >= alpha_sum
      constraint_alpha <- numeric(NVars)
      constraint_alpha[variable_indexer$gamma(s)] <- 1
      for(l in 1:K)
      {
        for(m in 1:P)
        {
          constraint_alpha[variable_indexer$eta_slm(s,l,m)] <- kappa[[s]][m]
        }
      }
      master_constraint_matrix <- rbind(master_constraint_matrix, constraint_alpha)
      master_rhs <- c(master_rhs, alpha_sum)
      master_senses <- c(master_senses, "G")
    }
    
    obj <- sum(J)
    if(obj < UB)
    {
      UB <- obj
      best_eta <- eta_l
      best_alpha <- alpha
      best_b <- b
      best_grouping <- grouping
    }
    
    LB <- sum(Gamma)
    
    Gap <- (UB-LB)/max(abs(LB),parameters$epsilon)
    
    if(abs(Gap-latest_gap) < 1e-5)
    {
      non_improving_iterations <- non_improving_iterations + 1
    }else
    {
      non_improving_iterations <- 0
    }
    
    objectives <- c(objectives, obj)
    
    info <- c(obj, UB, LB, Gap, sum(Gamma), sum(J), CPU_Master, sum(CPU_SVMs))
    iteration_info <- rbind(iteration_info, info)
    
    # if(iteration %% 5 == 1)
    print(sprintf("%d %f %f %f %f%% %f %f %f %f",iteration, obj, UB, LB, Gap*100, sum(Gamma), sum(J), CPU_Master, sum(CPU_SVMs)))
    print(grouping)
    
    if(Gap <= parameters$optimality_gap & iteration >= min_iterations)
    {
      terminate <- TRUE
    }
    
    if(parameters$iteration_count > 0 & iteration >= parameters$iteration_count)
    {
      terminate <- TRUE
    }
    
    total_time <- total_time + CPU_Master + sum(CPU_SVMs)
    
    if(parameters$time_limit > 0)
    {
      expected_end_time <- total_time + mean(iteration_info[, "Time(M)"] + iteration_info[, "Time(S)"])
      if(expected_end_time >= parameters$time_limit)
      {
        terminate <- TRUE
        print(sprintf("Time elapsed: %f, terminated at iteration %d due to time limit of %d", total_time, iteration, parameters$time_limit))
      }
    }
    
    if(non_improving_iterations > max_nonimproving_iterations)
    {
      terminate <- TRUE
      print(sprintf("Non-improving iterations: %d, terminated at iteration %d due to non-improving iterations limit of %d", non_improving_iterations, iteration, max_nonimproving_iterations))
    }
    
    # if(terminate & iteration %% 5 != 1)
    #   print(sprintf("%d %f %f %f %f%% %f %f %f %f", iteration, obj, UB, LB, Gap*100, sum(Gamma), sum(J), CPU_Master, sum(CPU_SVMs)))
    
    iteration <- iteration + 1
    
    latest_gap <- Gap
  }
  
  for(l in 1:K)
  {
    best_eta[[l]] <- normalize_eta(best_eta[[l]], parameters$epsilon)
  }
  
  state <- list(alpha = best_alpha, b = best_b, eta = best_eta, grouping = best_grouping, total_time = total_time,
                objectives = objectives, iteration_info = iteration_info, parameters = parameters)
  
  output <- list(state = state)
  
  return(output)
}
