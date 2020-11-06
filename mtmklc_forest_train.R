# Mehmet Gonen (mehmetgonen@ku.edu.tr)

mtmklc_forest_train <- function(Km, y, parameters, settings) {
  # parameters should have: K, optimality_gap
  
  TN <- length(Km)
  P <- dim(Km[[1]])[3]
  K <- parameters$K
  
  symmetry_breaking <- settings$symmetry_breaking
  
  terminate <- FALSE
  UB <- Inf
  LB <- -Inf
  
  best_eta <- NULL
  best_alpha <- NULL
  best_b <- NULL
  best_grouping <- NULL
  
  iteration <- 1
  
  #define constraint matrix here, add the equality constraint. Also define the objective function
  #order of variables: Gamma (TN), eta (TN,P), u_sq (TN,TN)
  NVars <- TN + TN*P + TN*(TN-1)/2
  
  variable_indexer <- list()
  variable_indexer$gamma_s <- function(s)
  {
    return(s)
  }
  variable_indexer$eta_sm <- function(s,m)
  {
    return(TN+(s-1)*P+m)
  }
  variable_indexer$u_sq <- function(s,q) #s<q
  {
    return( ifelse(s<q, TN+TN*P+(s-1)*(TN-s/2)+q-s, -1) )
  }
  
  master_obj_coef <- numeric(NVars)
  master_obj_coef[1:TN] <- 1
  
  master_constraint_matrix <- matrix(0, nrow = TN+TN*(TN-1)/2*P+1, ncol = NVars)
  
  #Eta's sum up to 1
  row <- 0
  for(s in 1:TN)
  {
    for(m in 1:P)
    {
      master_constraint_matrix[row+s,variable_indexer$eta_sm(s,m)] <- 1
    }
  }
  master_rhs <- rep(1,TN)
  master_senses <- rep("E",TN)
  
  #Eta and u relationship: eta_sm - eta_qm <= 1-u_sq
  row <- TN
  for(s in 1:(TN-1))
  {
    for(q in (s+1):TN)
    {
      for(m in 1:P)
      {
        row <- row+1
        master_constraint_matrix[row,variable_indexer$u_sq(s,q)] <- 1
        master_constraint_matrix[row,variable_indexer$eta_sm(s,m)] <- 1
        master_constraint_matrix[row,variable_indexer$eta_sm(q,m)] <- -1
      }
    }
  }
  master_rhs <- c(master_rhs, rep(1,TN*(TN-1)/2*P))
  master_senses <- c(master_senses,rep("L",TN*(TN-1)/2*P))
  
  #forest edges: sum_{s,q} u_sq = TN-K
  row = TN + TN*(TN-1)/2*P+1
  for(s in 1:(TN-1))
  {
    for(q in (s+1):TN)
    {
      master_constraint_matrix[row,variable_indexer$u_sq(s,q)] <- 1
    }
  }
  master_rhs <- c(master_rhs,TN-K)
  master_senses <- c(master_senses,"E")
  
  if(symmetry_breaking)
  {
    for(s in 1:TN)
    {
      if(s > 2) #sum_{q<s}u_{q,s} <=1
      {
        constraint_sb <- numeric(NVars)
        for(q in 1:(s-1))
        {
          constraint_sb[variable_indexer$u_sq(q,s)] <- 1
        }
        master_constraint_matrix <- rbind(master_constraint_matrix, constraint_sb)
        master_rhs <- c(master_rhs, 1)
        master_senses <- c(master_senses, "L")
      }
      
      if(s < TN-1) #sum_{q>s}u_{s,q} <=1
      {
        constraint_sb <- numeric(NVars)
        for(q in (s+1):TN)
        {
          constraint_sb[variable_indexer$u_sq(s,q)] <- 1
        }
        master_constraint_matrix <- rbind(master_constraint_matrix, constraint_sb)
        master_rhs <- c(master_rhs, 1)
        master_senses <- c(master_senses, "L")
      }
    }
  }
  
  master_lb <- numeric(NVars)
  master_ub <- rep(Inf, NVars)
  
  master_vtype <- rep("C", NVars)
  for(s in 1:(TN-1))
  {
    for(q in (s+1):TN)
    {
      master_vtype[variable_indexer$u_sq(s,q)] <- "B"
      master_ub[variable_indexer$u_sq(s,q)] <- 1
    }
  }
  
  Gamma <- numeric(TN)
  J <- numeric(TN)
  alpha <- vector("list", TN)
  b <- vector("list", TN)
  
  eta <- vector("list", TN)
  for (s in 1:TN) {
    eta[[s]] <- rep(1 / P, P)
  }
  
  cluster_size <- ceiling(TN/K)
  u_sq <- matrix(0, nrow = TN, ncol = TN)
  for (s in 1:(TN-1)) {
    if(floor((s-1)/cluster_size) == floor((s)/cluster_size))
      u_sq[s,s+1] <- 1
  }
  
  objectives <- c()
  total_time <- 0
  
  info_cols <- c("obj", "UB", "LB", "Gap","Gamma","J","Time(M)", "Time(S)", "Cycles")
  iteration_info <- matrix(0, 0, length(info_cols))
  colnames(iteration_info) <- info_cols
  
  min_iterations <- 6
  if(parameters$iteration_count > 0 & min_iterations > parameters$iteration_count)
    min_iterations <- parameters$iteration_count/2
  
  while(terminate == FALSE)
  {
    CPU_Master <- 0
    CPU_SVMs <- rep(0,TN)
    if(iteration > 1)
    {
      #SOLVE MP and obtain Gamma[s], eta[[s]], u_sq[s,q]
      master_result <- solve_mtmklc_forest_master_problem(master_obj_coef, master_constraint_matrix, 
                                                          master_rhs, master_senses, master_lb, master_ub, 
                                                          master_vtype, TN, P, variable_indexer)
      Gamma <- master_result$Gamma
      
      CPU_Master <- master_result$CPU
      
      eta <- master_result$eta
      
      u_sq <- master_result$u_sq
    }
    
    components <- find_components(u_sq)
    grouping <- rep(0, TN)
    for(l in 1:length(components))
    {
      grouping[components[[l]]] <- l
    }
    
    #subtour elimination
    cycles <- 0
    if(symmetry_breaking == FALSE){
      for(component in components)
      {
        if(sum(u_sq[component, component]) >= length(component)) #a cycle detected
        {
          cycles <- cycles + 1
          
          #add constraint Sum_{sq in cyle}u_sq <= |cycle| - 1
          constraint_subtour <- numeric(NVars)
          for(s in component)
          {
            for(q in component)
            {
              if(s < q)
              {
                constraint_subtour[variable_indexer$u_sq(s,q)] <- 1
              }
            }
          }
          master_constraint_matrix <- rbind(master_constraint_matrix, constraint_subtour)
          master_rhs <- c(master_rhs, length(component)-1)
          master_senses <- c(master_senses, "L")
        }
      }
    }
    
    kappa <- vector("list", TN)
    
    for (s in 1:TN) {
      Keta <- calculate_Keta(Km[[s]], eta[[s]])
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
      
      #add constraint Gamma_s + Sum_{m}Eta_sm*kappa_sm >= alpha_sum
      constraint_alpha <- numeric(NVars)
      constraint_alpha[variable_indexer$gamma_s(s)] <- 1
      for(m in 1:P)
      {
        constraint_alpha[variable_indexer$eta_sm(s,m)] <- kappa[[s]][m]
      }
      master_constraint_matrix <- rbind(master_constraint_matrix, constraint_alpha)
      master_rhs <- c(master_rhs, alpha_sum)
      master_senses <- c(master_senses, "G")
    }
    
    obj <- sum(J)
    if(obj < UB & cycles == 0)
    {
      UB <- obj
      best_eta <- eta
      best_alpha <- alpha
      best_b <- b
      best_grouping <- grouping
    }
    
    LB <- sum(Gamma)
    
    Gap <- (UB-LB)/max(abs(LB), parameters$epsilon)
    
    objectives <- c(objectives, obj)
    
    info <- c(obj, UB, LB, Gap, sum(Gamma), sum(J), CPU_Master, sum(CPU_SVMs), cycles)
    iteration_info <- rbind(iteration_info, info)
    
    # if(iteration %% 5 == 1)
    print(sprintf("%d %f %f %f %f%% %f %f %f %f %d", iteration, obj, UB, LB, Gap*100, sum(Gamma), sum(J), CPU_Master, sum(CPU_SVMs), cycles))
    
    if(Gap <= parameters$optimality_gap & iteration >= min_iterations & cycles == 0)
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
      expected_end_time <- total_time + sum(iteration_info[,c("Time(M)", "Time(S)")])/nrow(iteration_info)
      if(expected_end_time >= parameters$time_limit)
      {
        terminate <- TRUE
        print(sprintf("Time elapsed: %f, terminated at iteration %d due to time limit of %d", total_time, iteration, parameters$time_limit))
      }
    }
    
    if(terminate & iteration %% 5 != 1)
      print(sprintf("%d %f %f %f %f%% %f %f %f %f %d", iteration, obj, UB, LB, Gap*100, sum(Gamma), sum(J), CPU_Master, sum(CPU_SVMs), cycles))
    
    iteration <- iteration + 1
  }
  
  for(s in 1:TN)
  {
    best_eta[[s]] <- normalize_eta(best_eta[[s]], parameters$epsilon)
  }
  
  state <- list(alpha = best_alpha, b = best_b, eta = best_eta, grouping = best_grouping, total_time = total_time,
                objectives = objectives, iteration_info = iteration_info, parameters = parameters)
  
  output <- list(state = state)
  
  return(output)
}

mtmklc_forest_bdf_cumulative_train <- function(Km, y, parameters, settings) {
  
  two_sided_eta_diff <- settings$two_sided_eta_diff
  cumulative_eta_diff <- settings$cumulative_eta_diff
  
  optimality_cut_frequency <- settings$optimality_cut_frequency
  
  warmup_iterations <- settings$warmup_iterations
  
  beta_positivity_threshold <- settings$beta_positivity_threshold
  beta_steps_before_elimination <- settings$beta_steps_before_elimination
  
  compute_po_cut <- settings$compute_po_cut
  add_both_normal_and_po_cut <- settings$add_both_normal_and_po_cut
  core_point <- settings$core_point
  
  symmetry_breaking <- settings$symmetry_breaking
  
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
  #order of variables: \mu, u_sq (TN,TN)
  MP_NVar <- 1 + TN*(TN-1)/2
  
  MP_variable_indexer <- list()
  MP_variable_indexer$mu <- function()
  {
    return(1)
  }
  MP_variable_indexer$u_sq <- function(s,q) #s<q
  {
    return( ifelse(s<q, 1 + (s-1)*(TN-s/2)+q-s, -1) )
  }
  
  MP_obj <- numeric(MP_NVar)
  MP_obj[MP_variable_indexer$mu()] <- 1
  
  #forest edges: sum_{s,q} u_sq = TN-K
  MP_constraints <- matrix(1, nrow = 1, ncol = MP_NVar)
  MP_constraints[1, MP_variable_indexer$mu()] <- 0
  MP_rhs <- c(TN-K)
  MP_senses <- c("E")
  
  if(symmetry_breaking)
  {
    for(s in 1:TN)
    {
      if(s > 2) #sum_{q<s}u_{q,s} <=1
      {
        constraint_sb <- numeric(MP_NVar)
        for(q in 1:(s-1))
        {
          constraint_sb[MP_variable_indexer$u_sq(q,s)] <- 1
        }
        MP_constraints <- rbind(MP_constraints, constraint_sb)
        MP_rhs <- c(MP_rhs, 1)
        MP_senses <- c(MP_senses, "L")
      }
      
      if(s < TN-1) #sum_{q>s}u_{s,q} <=1
      {
        constraint_sb <- numeric(MP_NVar)
        for(q in (s+1):TN)
        {
          constraint_sb[MP_variable_indexer$u_sq(s,q)] <- 1
        }
        MP_constraints <- rbind(MP_constraints, constraint_sb)
        MP_rhs <- c(MP_rhs, 1)
        MP_senses <- c(MP_senses, "L")
      }
    }
  }
  
  MP_lb <- rep(0, MP_NVar)
  MP_ub <- rep(1, MP_NVar)
  MP_vtype <- rep("B", MP_NVar)
  
  MP_vtype[MP_variable_indexer$mu()] <- "C"
  MP_ub[MP_variable_indexer$mu()] <- Inf
  
  svm_solutions_info = vector("list", TN)
  for (s in 1:TN) {
    svm_solutions_info[[s]] <- matrix(0, nrow = 0, ncol = P+1) #1 column per each k_m, 1 column for sum_alpha: k_1, ..., k_P, sum_alpha
  }
  
  DSP_Nvar <- TN + P*TN*(TN-1) #+ betas which will be added iteratively
  if(!two_sided_eta_diff)
    DSP_Nvar <- TN + P*TN*(TN-1)/2
  
  DSP_Alpha_Index <- vector("list", TN)
  for (s in 1:TN) {
    alphas <- nrow(svm_solutions_info[[s]])
    indices <- numeric(alphas)
    
    if(alphas > 0)
    {
      indices <- DSP_Nvar + (1:alphas)
      DSP_Nvar <- DSP_Nvar + alphas
    }
    
    DSP_Alpha_Index[[s]] <- indices
  }
  
  DSP_variable_indexer <- list()
  DSP_variable_indexer$theta_s <- function(s)
  {
    return(s)
  }
  DSP_variable_indexer$lambda_sqm <- function(s,q,m)
  {
    if(two_sided_eta_diff)
    {
      if(s == q) 
        return(-1)
      
      if(s < q)
        return(TN + (m-1)*TN*(TN-1) + (s-1)*(TN-1) + q-1)
      
      return(TN + (m-1)*TN*(TN-1) + (s-1)*(TN-1) + q)
    }else  #we must have s<q
    {
      if(q <= s)
        return(-1)
      
      return(TN + (m-1)*TN*(TN-1)/2 + (s-1)*(TN-s/2)+q-s)
    }
  }
  
  DSP_variable_indexer$beta_sa <- function(s,a)
  {
    if(a > length(DSP_Alpha_Index[[s]]))
      return(-1)
    
    return(DSP_Alpha_Index[[s]][a])
  }
  
  DSP_row_indexer <- list()
  DSP_row_indexer$gamma_s <- function(s)
  {
    return(s)
  }
  DSP_row_indexer$eta_sm <- function(s,m)
  {
    return(TN + (s-1)*P+m)
  }
  DSP_row_indexer$po_optimality <- function()
  {
    return(TN + TN*P + 1)
  }
  
  DSP_obj <- numeric(DSP_Nvar)
  for(s in 1:TN)
  {
    DSP_obj[DSP_variable_indexer$theta_s(s)] <- 1
    
    if(length(DSP_Alpha_Index[[s]]) > 0)
    {
      for(a in 1:length(DSP_Alpha_Index[[s]]))
        DSP_obj[DSP_variable_indexer$beta_sa(s,a)] <- svm_solutions_info[[s]][a,P+1]
    }
  }
  #coefficients of lambda's are not known yet
  
  DSP_Nrow <- TN + TN*P
  DSP_constraints <- matrix(0, nrow = DSP_Nrow, ncol = DSP_Nvar)
  DSP_rhs <- rep(0, DSP_Nrow)
  DSP_senses <- rep("L", DSP_Nrow)
  
  for(s in 1:TN)
  { 
    row <- DSP_row_indexer$gamma_s(s)
    if(length(DSP_Alpha_Index[[s]]) > 0)
    {
      for(a in 1:length(DSP_Alpha_Index[[s]]))
        DSP_constraints[row, DSP_variable_indexer$beta_sa(s,a)] <- 1
    }
    DSP_rhs[row] <- 1 
  }
  
  for(s in 1:TN)
  {
    for(m in 1:P)
    {
      row <- DSP_row_indexer$eta_sm(s,m)
      DSP_constraints[row, DSP_variable_indexer$theta_s(s)] <- 1
      
      m_range <- c(m)
      if(cumulative_eta_diff)
      {
        m_range <- m:P
      }
      
      for(n in m_range)
      {
        for(q in 1:TN)
        {
          if(two_sided_eta_diff)
          {
            if(s!=q)
            {
              DSP_constraints[row, DSP_variable_indexer$lambda_sqm(s,q,n)] <- -1
              DSP_constraints[row, DSP_variable_indexer$lambda_sqm(q,s,n)] <- +1
            }
          }else
          {
            if(q > s)
            {
              DSP_constraints[row, DSP_variable_indexer$lambda_sqm(s,q,m)] <- -1
            }else if(q < s)
            {
              DSP_constraints[row, DSP_variable_indexer$lambda_sqm(q,s,m)] <- +1
            }
          }
        }
      }
      
      if(length(DSP_Alpha_Index[[s]]) > 0)
      {
        for(a in 1:length(DSP_Alpha_Index[[s]]))
          DSP_constraints[row, DSP_variable_indexer$beta_sa(s,a)] <- svm_solutions_info[[s]][a,m]
      }
    }
  }
  
  DSP_lb <- rep(0, DSP_Nvar)
  for(s in 1:TN)
    DSP_lb[DSP_variable_indexer$theta_s(s)] <- -Inf
  DSP_ub <- rep(Inf, DSP_Nvar)
  DSP_vtype <- rep("C", DSP_Nvar)
  
  DSP_variables_per_iteration <- matrix(0, nrow = 0, ncol = 3)
  colnames(DSP_variables_per_iteration) <-c("All", "U0 Elimination", "Beta Elimination")
  
  Mu <- 0
  Gamma <- numeric(TN)
  J <- numeric(TN)
  alpha <- vector("list", TN)
  b <- vector("list", TN)
  
  eta <- vector("list", TN)
  for (s in 1:TN) {
    eta[[s]] <- rep(1 / P, P)
  }
  
  cluster_size <- ceiling(TN/K)
  u_sq <- matrix(0, nrow = TN, ncol = TN)
  for (s in 1:(TN-1)) {
    if(floor((s-1)/cluster_size) == floor((s)/cluster_size))
      u_sq[s,s+1] <- 1
  }
  
  betas_track <-  vector("list", TN)
  betas_zero_count <- vector("list", TN)
  for(s in 1:TN)
  {
    betas_track[[s]] <- matrix(0, nrow = 0, ncol = 1)
    betas_zero_count[[s]] <- c(0)
  }
  
  MP_objectives <- c()
  total_time <- 0
  
  info_cols <- c("Obj", "UB", "LB", "Gap","Gamma","J","Time(M)", "Time(DSP)", "Time(PO)", "Time(SVM)", "Cycles")
  iteration_info <- matrix(0, 0, length(info_cols))
  colnames(iteration_info) <- info_cols
  
  min_iterations <- 6
  if(parameters$iteration_count > 0 & min_iterations > parameters$iteration_count)
    min_iterations <- parameters$iteration_count/2
  
  initial_DSP_sol <- c()
  
  while(terminate == FALSE)
  {
    CPU_MP <- 0
    CPU_SVMs <- rep(0,TN)
    CPU_DSP <- 0
    CPU_PO <- 0
    
    solve_mp <- iteration > warmup_iterations
    
    if(solve_mp)
    {
      #SOLVE MP and obtain mu and u_sq[s,q]
      master_result <- solve_mtmklc_bdf_master_problem(MP_obj, MP_constraints, MP_rhs, MP_senses, 
                                                       MP_lb, MP_ub, MP_vtype, 
                                                       TN, MP_variable_indexer, LP_relaxation = FALSE)
      
      Mu <- master_result$Mu
      CPU_MP <- master_result$CPU
      u_sq <- round(master_result$u_sq, 0)
    }
    
    components <- find_components(u_sq)
    grouping <- rep(0, TN)
    for(l in 1:length(components))
    {
      grouping[components[[l]]] <- l
    }
    
    #subtour elimination
    cycles <- 0
    if(symmetry_breaking == FALSE){
      for(component in components)
      {
        if(sum(u_sq[component, component]) >= length(component)) #a cycle detected
        {
          cycles <- cycles + 1
          
          #add constraint Sum_{sq in cyle}u_sq <= |cycle| - 1
          constraint_subtour <- numeric(MP_NVar)
          for(s in component)
          {
            for(q in component)
            {
              if(s < q)
              {
                constraint_subtour[MP_variable_indexer$u_sq(s,q)] <- 1
              }
            }
          }
          MP_constraints <- rbind(MP_constraints, constraint_subtour)
          MP_rhs <- c(MP_rhs, length(component)-1)
          MP_senses <- c(MP_senses, "L")
        }
      }
    }
    
    if(cycles == 0) #u is a feasible solution
    {
      print(sprintf("** Groupings: {%s} groups = %d", paste(grouping, collapse = " "), max(grouping)))
      
      inner_iterations <- optimality_cut_frequency
      if(iteration <= warmup_iterations)
        inner_iterations <- 1
      
      for(dsp_i in 1:inner_iterations)
      {
        if(iteration > 1 | dsp_i > 1)
        {
          alphas <- sapply(DSP_Alpha_Index, length)
          
          # edge_indexing <- FALSE
          betas_eliminated <- 0
          if(beta_steps_before_elimination > 0)
          {
            for(s in 1:TN)
            {
              for(a in 1:length(betas_zero_count[[s]]))
              {
                if(betas_zero_count[[s]][a] > beta_steps_before_elimination)
                {
                  DSP_ub[DSP_variable_indexer$beta_sa(s,a)] <- 0
                  betas_eliminated <- betas_eliminated + 1
                }
              }
            }
          }
          
          DSP_u1_columns <- rep(x = TRUE, DSP_Nvar)
          
          for(s in 1:TN)
          {
            for(q in 1:TN)
            {
              if(q!=s)
              {
                usq <- u_sq[s, q] + u_sq[q, s]
                for(m in 1:P)
                {
                  v <- DSP_variable_indexer$lambda_sqm(s,q,m)
                  if(v <= 0)
                  {
                    print(c(s,q,m))
                  }
                  if(v > 0) #it might be -1 under one-sided diff if q > s
                  {
                    DSP_obj[v] <- usq - 1
                    DSP_u1_columns[v] <- usq == 1
                  }
                }
              }
            }
          }
          
          current_DSP_obj <- DSP_obj#[DSP_u1_columns]
          current_DSP_constraints <- DSP_constraints#[,DSP_u1_columns]
          current_DSP_lb <- DSP_lb#[DSP_u1_columns]
          current_DSP_ub <- DSP_ub#[DSP_u1_columns]
          current_DSP_ub[!DSP_u1_columns] <- 0
          current_DSP_vtype <- DSP_vtype#[DSP_u1_columns]
          
          DSP_variables_per_iteration <- rbind(DSP_variables_per_iteration, c(length(DSP_ub), 
                                                                              length(DSP_ub)-sum(!DSP_u1_columns),
                                                                              betas_eliminated))
          
          if(length(initial_DSP_sol) > 0)
          {
            tmp <- numeric(DSP_Nvar)
            tmp[1:length(initial_DSP_sol)] <- initial_DSP_sol
            initial_DSP_sol <- tmp
          }
          
          #SOLVE DSP and obtain (Gamma, eta) and (theta, lambda, beta)
          DSP_result <<- solve_model3Tree_BD_DSP_problem(current_DSP_obj, current_DSP_constraints, DSP_rhs, DSP_senses, 
                                                         current_DSP_lb, current_DSP_ub, current_DSP_vtype, 
                                                         TN, P, alphas, DSP_variable_indexer, DSP_row_indexer, initial_DSP_sol)
          
          initial_DSP_sol <- DSP_result$dual_solution$Sol
          
          eta <- DSP_result$primal_solution$Eta
          Gamma <- DSP_result$primal_solution$Gamma
          CPU_DSP <- CPU_DSP + DSP_result$CPU
          sub_objective <- DSP_result$objective
          
          for(s in 1:TN)
          {
            while(length(DSP_result$dual_solution$Beta[[s]]) > ncol(betas_track[[s]]))
            {
              betas_track[[s]] <- cbind(betas_track[[s]], numeric(nrow(betas_track[[s]])))
            }
            
            betas_track[[s]] <- rbind(betas_track[[s]], DSP_result$dual_solution$Beta[[s]])
            
            B0 <- 1*(DSP_result$dual_solution$Beta[[s]] < beta_positivity_threshold)
            reset <- which(B0 == 0)
            B0[1:length(betas_zero_count[[s]])] <- B0[1:length(betas_zero_count[[s]])] + betas_zero_count[[s]]
            betas_zero_count[[s]] <- B0
            betas_zero_count[[s]][reset] <- 0
          }
          
          if(dsp_i == inner_iterations)
          {
            optimality_cuts <- list()
            
            optimality_cuts[[1]] <- list(Theta = DSP_result$dual_solution$Theta,
                                         Lambda = DSP_result$dual_solution$Lambda,
                                         Beta = DSP_result$dual_solution$Beta)
            
            if(iteration > warmup_iterations)
            {
              if(compute_po_cut)
              {
                #SOLVE PO and obtain (theta, lambda, beta) 
                PO_obj <- DSP_obj
                
                for(s in 1:TN)
                {
                  for(q in 1:TN)
                  {
                    if(q!=s)
                    {
                      iusq <- core_point[min(s,q), max(s,q)]
                      for(m in 1:P)
                      {
                        v <- DSP_variable_indexer$lambda_sqm(s,q,m)
                        if(v > 0) #it might be -1 under one-sided diff
                          PO_obj[v] <- iusq- 1
                      }
                    }
                  }
                }
                
                PO_constraints <- rbind(DSP_constraints, DSP_obj)
                PO_rhs <- c(DSP_rhs, sub_objective-0.00001)
                PO_senses <- c(DSP_senses, "G")
                
                PO_result <- solve_model3Tree_BD_DSP_problem(PO_obj, PO_constraints, PO_rhs, PO_senses, 
                                                             DSP_lb, DSP_ub, DSP_vtype, 
                                                             TN, P, alphas, DSP_variable_indexer, DSP_row_indexer,
                                                             initial_DSP_sol)
                
                PO_eta <- PO_result$primal_solution$Eta
                PO_Gamma <- PO_result$primal_solution$Gamma
                CPU_PO <- PO_result$CPU
                PO_objective <- PO_result$objective
                
                po_cut <- list(Theta = PO_result$dual_solution$Theta,
                               Lambda = PO_result$dual_solution$Lambda,
                               Beta = PO_result$dual_solution$Beta)
                
                if(add_both_normal_and_po_cut)
                {
                  optimality_cuts[[2]] <- po_cut
                }else
                {
                  optimality_cuts[[1]] <- po_cut
                }
              }
              
              for(oc in 1:length(optimality_cuts))
              {
                Theta <- optimality_cuts[[oc]]$Theta
                Lambda <- optimality_cuts[[oc]]$Lambda
                Beta <- optimality_cuts[[oc]]$Beta
                
                ds_constant <- sum(Theta)
                for(s in 1:TN)
                {
                  if(length(DSP_Alpha_Index[[s]]) > 0)
                  {
                    for(a in 1:length(DSP_Alpha_Index[[s]]))
                      ds_constant <- ds_constant + svm_solutions_info[[s]][a,P+1]*Beta[[s]][a]
                  }
                }
                
                bd_cut <- numeric(MP_NVar)
                bd_cut[MP_variable_indexer$mu()] <- 1
                for(s in 1:(TN-1)) #yes still correct
                {
                  for(q in (s+1):TN)
                  {
                    l_sum <- 0
                    for(m in 1:P)
                    {
                      l_sum <- l_sum + Lambda[[m]][s,q] + Lambda[[m]][q,s] #if not two-sided, Lambda[[m]][q,s] will be 0 any way
                    }
                    
                    bd_cut[MP_variable_indexer$u_sq(s,q)] <- -l_sum
                    ds_constant <- ds_constant - l_sum
                  }
                }
                
                MP_constraints <- rbind(MP_constraints, bd_cut)
                # print(c(bd_cut, ds_constant))
                MP_rhs <- c(MP_rhs, ds_constant)
                MP_senses <- c(MP_senses, "G") 
              }
            }
          }
        }
        
        # SOLVE SVMs and obtain alpha and J
        kappa <- vector("list", TN)
        
        for (s in 1:TN) {
          
          Keta <- calculate_Keta(Km[[s]], eta[[s]])
          C <- parameters$C
          if(parameters$normalize_C)
            C <- C / length(y[[s]])
          svm_result <<- solve_classification_svm(Keta, y[[s]], C, parameters$epsilon)
          
          CPU_SVMs[s] <- CPU_SVMs[s] + svm_result$CPU
          
          alpha[[s]] <- svm_result$alpha
          b[[s]] <- svm_result$b
          
          alpha_cut <- svm_result$alpha_original
          J[s] <- svm_result$objective_original
          
          #model$alpha is already multiplied by y
          alpha_sum <- t(alpha_cut) %*% y[[s]]
          kappa[[s]] <- calucalte_kappa(alpha_cut, Km[[s]], P)
          
          svm_solutions_info[[s]] <- rbind(svm_solutions_info[[s]], c(kappa[[s]], alpha_sum))
          
          DSP_Nvar <- DSP_Nvar + 1
          DSP_constraints <- cbind(DSP_constraints, numeric(DSP_Nrow))
          DSP_Alpha_Index[[s]] <- c(DSP_Alpha_Index[[s]], DSP_Nvar)
          DSP_obj <- c(DSP_obj, alpha_sum)
          DSP_constraints[DSP_row_indexer$gamma_s(s), DSP_Nvar] <- 1
          for(m in 1:P)
            DSP_constraints[DSP_row_indexer$eta_sm(s,m), DSP_Nvar] <- kappa[[s]][m]
          DSP_lb <- c(DSP_lb, 0)
          DSP_ub <- c(DSP_ub, Inf)
          DSP_vtype <- c(DSP_vtype, "C")
        }
      }
    }
    
    obj <- sum(J)
    if(obj < UB & cycles == 0)
    {
      UB <- obj
      best_eta <- eta
      best_alpha <- alpha
      best_b <- b
      best_grouping <- grouping
    }
    
    LB <- max(Mu, LB) #for numerical stability
    
    Gap <- (UB-LB)/max(abs(LB), parameters$epsilon)
    
    MP_objectives <- c(MP_objectives, obj)
    
    info <- c(obj, UB, LB, Gap, sum(Gamma), sum(J), CPU_MP, CPU_DSP, CPU_PO, sum(CPU_SVMs), cycles)
    iteration_info <- rbind(iteration_info, info)
    
    # if(iteration %% 5 == 1)
    print(sprintf("%d %f %f %f %f%% %f %f %f %f %f %f %d", iteration, obj, UB, LB, Gap*100, sum(Gamma), sum(J), CPU_MP, CPU_DSP, CPU_PO, sum(CPU_SVMs), cycles))
    
    if(Gap <= parameters$optimality_gap & iteration >= min_iterations & cycles == 0)
    {
      terminate <- TRUE
    }
    
    if(parameters$iteration_count > 0 & iteration >= parameters$iteration_count)
    {
      terminate <- TRUE
    }
    
    total_time <- total_time + CPU_MP + sum(CPU_SVMs) + CPU_DSP + CPU_PO
    
    if(parameters$time_limit > 0)
    {
      expected_end_time <- total_time + sum(iteration_info[,c("Time(M)", "Time(DSP)", "Time(PO)", "Time(SVM)")])/nrow(iteration_info)
      if(expected_end_time >= parameters$time_limit)
      {
        terminate <- TRUE
        print(sprintf("Time elapsed: %f, terminated at iteration %d due to time limit of %d", total_time, iteration, parameters$time_limit))
      }
    }
    
    iteration <- iteration + 1
  }
  
  for(s in 1:TN)
  {
    best_eta[[s]] <- normalize_eta(best_eta[[s]], parameters$epsilon)
  }
  
  state <- list(alpha = best_alpha, b = best_b, eta = best_eta, grouping = best_grouping, total_time = total_time,
                MP_objectives = MP_objectives, iteration_info = iteration_info, parameters = parameters, 
                DSP_variables_per_iteration = DSP_variables_per_iteration)
  
  output <- list(state = state)
  
  return(output)
}

mtmklc_forest_bdf_strengthened_train <- function(Km, y, parameters, initial_u = NULL, settings){
  
  eliminate_beta_permanantly <- TRUE
  
  optimality_cut_frequency <- settings$optimality_cut_frequency
  
  warmup_iterations <- settings$warmup_iterations
  
  beta_positivity_threshold <- settings$beta_positivity_threshold
  beta_steps_before_elimination <- settings$beta_steps_before_elimination
  
  symmetry_breaking <- settings$symmetry_breaking
  
  compute_po_cut <- settings$compute_po_cut
  add_both_normal_and_po_cut <- settings$add_both_normal_and_po_cut
  core_point <- settings$core_point
  
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
  best_u <- NULL
  
  iteration <- 1
  
  #define constraint matrix here, add the equality constraint. Also define the objective function
  #order of variables: \mu, u_sq (TN,TN)
  MP_NVar <- 1 + TN*(TN-1)/2
  
  MP_variable_indexer <- list()
  MP_variable_indexer$mu <- function()
  {
    return(1)
  }
  MP_variable_indexer$u_sq <- function(s,q) #s<q
  {
    return( ifelse(s<q, 1 + (s-1)*(TN-s/2)+q-s, -1) )
  }
  
  MP_obj <- numeric(MP_NVar)
  MP_obj[MP_variable_indexer$mu()] <- 1
  
  #forest edges: sum_{s,q} u_sq = TN-K
  MP_constraints <- matrix(1, nrow = 1, ncol = MP_NVar)
  MP_constraints[1, MP_variable_indexer$mu()] <- 0
  MP_rhs <- c(TN-K)
  MP_senses <- c("E")
  
  if(symmetry_breaking)
  {
    for(s in 1:TN)
    {
      if(s > 2) #sum_{q<s}u_{q,s} <=1
      {
        constraint_sb <- numeric(MP_NVar)
        for(q in 1:(s-1))
        {
          constraint_sb[MP_variable_indexer$u_sq(q,s)] <- 1
        }
        MP_constraints <- rbind(MP_constraints, constraint_sb)
        MP_rhs <- c(MP_rhs, 1)
        MP_senses <- c(MP_senses, "L")
      }
      
      if(s < TN-1) #sum_{q>s}u_{s,q} <=1
      {
        constraint_sb <- numeric(MP_NVar)
        for(q in (s+1):TN)
        {
          constraint_sb[MP_variable_indexer$u_sq(s,q)] <- 1
        }
        MP_constraints <- rbind(MP_constraints, constraint_sb)
        MP_rhs <- c(MP_rhs, 1)
        MP_senses <- c(MP_senses, "L")
      }
    }
  }
  
  MP_lb <- rep(0, MP_NVar)
  MP_ub <- rep(1, MP_NVar)
  MP_vtype <- rep("B", MP_NVar)
  
  MP_vtype[MP_variable_indexer$mu()] <- "C"
  MP_ub[MP_variable_indexer$mu()] <- Inf
  
  svm_info <- matrix(0, nrow = P+3, ncol = 0) #P kappas + 1 |alpha| + 1 type + 1 for nonzero_frequency
  get_beta_val_range <- function(s)
  {
    return(which(svm_info[P+2,] == s))
  }
  
  #        theta + lambda      + omega and nu (each P*TN*(TN-1)/2)
  DSP_Nvar <- TN + TN*(TN-1)/2 + 2*P*TN*(TN-1)/2 #+ betas which will be added iteratively
  
  DSP_variable_indexer <- list()
  DSP_variable_indexer$theta_s <- function(s)
  {
    return(s)
  }
  DSP_variable_indexer$lambda_sq <- function(s,q)
  {
    if(q <= s)
      return(-1)
    
    return(TN + (s-1)*(TN-s/2)+q-s)
  }
  DSP_variable_indexer$nu_sqm <- function(s,q,m)
  {
    if(q == s)
      return(-1)
    
    if(s < q)
      return(TN + TN*(TN-1)/2 + (m-1)*TN*(TN-1) + (s-1)*(TN-1) + q-1) #nu
    
    return(TN + TN*(TN-1)/2 + (m-1)*TN*(TN-1) + (s-1)*(TN-1) + q) #omega
  }
  DSP_variable_indexer$beta_range_s <- function(s)
  {
    all_s <- get_beta_val_range(s)
    if(length(all_s) == 0)
      return(c())
    return(TN + TN*(TN-1)/2 + P*TN*(TN-1) + all_s)
  }
  DSP_variable_indexer$beta_sa <- function(s,a)
  {
    range_s <- DSP_variable_indexer$beta_range_s(s)
    if(a > length(range_s))
      return(-1)
    
    return(range_s[a])
  }
  
  DSP_row_indexer <- list()
  DSP_row_indexer$gamma_s <- function(s)
  {
    return(s)
  }
  DSP_row_indexer$eta_sm <- function(s,m)
  {
    return(TN + (s-1)*P+m)
  }
  DSP_row_indexer$phi_sqm <- function(s,q,m)
  {
    if(q <= s)
      return(-1)
    
    return(TN + TN*P + (m-1)*TN*(TN-1)/2 + (s-1)*(TN-s/2)+q-s)
  }
  DSP_row_indexer$po_optimality <- function()
  {
    return(TN + TN*P + P*TN*(TN-1)/2 + 1)
  }
  
  DSP_obj <- numeric(DSP_Nvar)
  for(s in 1:TN)
  {
    DSP_obj[DSP_variable_indexer$theta_s(s)] <- 1
    
    range_s_var <- DSP_variable_indexer$beta_range_s(s)
    range_s_val <- get_beta_val_range(s)
    if(length(range_s_val) > 0)
    {
      DSP_constraints[row, range_s_var] <- svm_info[P+1, range_s_val]
    }
  }
  #coefficients of lambda's are not known yet
  
  DSP_Nrow <- DSP_row_indexer$po_optimality() - 1
  DSP_constraints <- matrix(0, nrow = DSP_Nrow, ncol = DSP_Nvar)
  DSP_rhs <- rep(0, DSP_Nrow)
  DSP_senses <- rep("L", DSP_Nrow)
  
  for(s in 1:TN)
  { 
    row <- DSP_row_indexer$gamma_s(s)
    range_s <- DSP_variable_indexer$beta_range_s(s)
    if(length(range_s) > 0)
    {
      DSP_constraints[row, range_s] <- 1
    }
    DSP_rhs[row] <- 1 
  }
  
  for(s in 1:TN)
  {
    for(m in 1:P)
    {
      row <- DSP_row_indexer$eta_sm(s,m)
      DSP_constraints[row, DSP_variable_indexer$theta_s(s)] <- 1
      
      for(q in 1:TN)
      {
        if(s != q)
        {
          DSP_constraints[row, DSP_variable_indexer$nu_sqm(s,q,m)] <- +1
        }
      }
      
      range_s_var <- DSP_variable_indexer$beta_range_s(s)
      range_s_val <- get_beta_val_range(s)
      if(length(range_s_val) > 0)
      {
        DSP_constraints[row, range_s_var] <- svm_info[m, range_s_val]
      }
    }
  }
  
  for(s in 1:(TN-1))
  {
    for(q in (s+1):TN)
    {
      for(m in 1:P)
      {
        row <- DSP_row_indexer$phi_sqm(s,q,m)
        DSP_constraints[row, DSP_variable_indexer$nu_sqm(s,q,m)] <- -1
        DSP_constraints[row, DSP_variable_indexer$nu_sqm(q,s,m)] <- -1
        DSP_constraints[row, DSP_variable_indexer$lambda_sq(s,q)] <- +1
        # DSP_rhs[row] <- 0 #already 0
      }
    }
  }
  
  DSP_lb <- rep(0, DSP_Nvar)
  for(s in 1:TN)
    DSP_lb[DSP_variable_indexer$theta_s(s)] <- -Inf
  
  DSP_ub <- rep(Inf, DSP_Nvar)
  DSP_vtype <- rep("C", DSP_Nvar)
  
  DSP_variables_per_iteration <- matrix(0, nrow = 0, ncol = 3)
  colnames(DSP_variables_per_iteration) <-c("All", "U0 Elimination", "Beta Elimination")
  
  Mu <- 0
  Gamma <- numeric(TN)
  J <- numeric(TN)
  alpha <- vector("list", TN)
  b <- vector("list", TN)
  
  eta <- vector("list", TN)
  for (s in 1:TN) {
    eta[[s]] <- rep(1 / P, P)
  }
  
  u_sq <- matrix(0, nrow = TN, ncol = TN)
  if(!is.null(initial_u))
  {
    u_sq <- initial_u
  }else
  {
    cluster_size <- ceiling(TN/K)
    for (s in 1:(TN-1)) {
      if(floor((s-1)/cluster_size) == floor((s)/cluster_size))
        u_sq[s,s+1] <- 1
    }
  }
  
  components <- find_components(u_sq)
  cycles <- 0
  grouping <- rep(0, TN)
  for(l in 1:length(components))
  {
    grouping[components[[l]]] <- l
  }
  
  betas_track <-  vector("list", TN)
  for(s in 1:TN)
  {
    betas_track[[s]] <- matrix(0, nrow = 0, ncol = 1)
  }
  
  MP_objectives <- c()
  total_time <- 0
  
  info_cols <- c("Obj", "UB", "LB", "Gap","Gamma","J","Time(M)", "Time(DSP)", "Time(PO)", "Time(SVM)", "Cycles")
  iteration_info <- matrix(0, 0, length(info_cols))
  colnames(iteration_info) <- info_cols
  
  min_iterations <- 6
  if(min_iterations > parameters$iteration_count)
    min_iterations <- parameters$iteration_count/2
  
  initial_DSP_sol <- c()
  
  new_subtours <- list()
  
  while(terminate == FALSE)
  {
    CPU_MP <- 0
    CPU_SVMs <- rep(0,TN)
    CPU_DSP <- 0
    CPU_PO <- 0
    
    solve_mp <- iteration > warmup_iterations
    
    if(solve_mp)
    {
      #SOLVE MP and obtain mu and u_sq[s,q]
      master_result <- solve_model3Tree_BD_master_problem(MP_obj, MP_constraints, MP_rhs, MP_senses, 
                                                          MP_lb, MP_ub, MP_vtype, 
                                                          TN, MP_variable_indexer, LP_relaxation = FALSE)
      
      Mu <- master_result$Mu
      
      CPU_MP <- master_result$CPU
      u_sq <- round(master_result$u_sq, 0)
      
      components <- find_components(u_sq)
      grouping <- rep(0, TN)
      for(l in 1:length(components))
      {
        grouping[components[[l]]] <- l
      }
      
      cycles <- 0
      if(symmetry_breaking == FALSE) {
        #subtour elimination
        for(component in components)
        {
          if(sum(u_sq[component, component]) >= length(component)) #a cycle detected
          {
            cycles <- cycles + 1
            
            new_subtours <- c(new_subtours, list(component))
            
            #add constraint Sum_{sq in cyle}u_sq <= |cycle| - 1
            constraint_subtour <- numeric(MP_NVar)
            for(s in component)
            {
              for(q in component)
              {
                if(s < q)
                {
                  constraint_subtour[MP_variable_indexer$u_sq(s,q)] <- 1
                }
              }
            }
            MP_constraints <- rbind(MP_constraints, constraint_subtour)
            MP_rhs <- c(MP_rhs, length(component)-1)
            MP_senses <- c(MP_senses, "L")
          }
        }
      }
    }
    
    if(cycles == 0) #u is a feasible solution
    {
      print(sprintf("** Groupings: {%s} groups = %d", paste(grouping, collapse = " "), max(grouping)))
    }else
    {
      print(sprintf("Groupings: {%s} groups = %d", paste(grouping, collapse = " "), max(grouping)))
    }
    
    if(cycles == 0) #u is a feasible solution
    {
      inner_iterations <- optimality_cut_frequency
      if(iteration <= warmup_iterations)
        inner_iterations <- 1
      
      for(dsp_i in 1:inner_iterations)
      {
        if(iteration > 1 | dsp_i > 1)
        {
          betas_eliminated <- 0
          if(beta_steps_before_elimination > 0)
          {
            to_eliminate_val_range <- which(svm_info[P+3,] > beta_steps_before_elimination)
            betas_eliminated <- length(to_eliminate_val_range)
            if(betas_eliminated > 0)
            {
              to_eliminate_var_range <- to_eliminate_val_range + (DSP_Nvar-ncol(svm_info))
              
              #TO DO: Permanantly eliminate these betas
              
              if(eliminate_beta_permanantly)
              {
                remanining_vars <- setdiff(1:DSP_Nvar, to_eliminate_var_range)
                DSP_ub <- DSP_ub[remanining_vars]
                DSP_lb <- DSP_lb[remanining_vars]
                DSP_vtype <- DSP_vtype[remanining_vars]
                DSP_obj <- DSP_obj[remanining_vars]
                DSP_constraints <- DSP_constraints[,remanining_vars]
                # print(c(DSP_Nvar - length(to_eliminate_var_range), length(DSP_obj)))
                DSP_Nvar <- length(DSP_obj)
                
                remanining_vals <- setdiff(1:ncol(svm_info), to_eliminate_val_range)
                svm_info <- svm_info[,remanining_vals]
                
                if(length(initial_DSP_sol) > 0)
                {
                  initial_DSP_sol <- initial_DSP_sol[remanining_vars]
                }
              }else
              {
                DSP_ub[to_eliminate_var_range] <- 0
              }
            }
          }
          
          alphas <- sapply(1:TN, function(s){length(DSP_variable_indexer$beta_range_s(s))})
          
          DSP_u1_columns <- rep(x = TRUE, DSP_Nvar)
          
          for(s in 1:(TN-1))
          {
            for(q in (s+1):TN)
            {
              v <- DSP_variable_indexer$lambda_sq(s,q)
              DSP_obj[v] <- u_sq[s, q]
              if(u_sq[s, q] < 0.0001)
                DSP_u1_columns[v] <- FALSE
            }
          }
          
          for(s in 1:TN)
          {
            for(q in 1:TN)
            {
              if(q!=s)
              {
                if(u_sq[s, q] + u_sq[q, s] < 0.001)
                {
                  for(m in 1:P)
                  {
                    #no objective function coefficient...
                    DSP_u1_columns[DSP_variable_indexer$nu_sqm(s,q,m)] <- FALSE
                  }
                }
              }
            }
          }
          
          current_DSP_obj <- DSP_obj#[DSP_u1_columns]
          current_DSP_constraints <- DSP_constraints#[,DSP_u1_columns]
          current_DSP_lb <- DSP_lb#[DSP_u1_columns]
          current_DSP_ub <- DSP_ub#[DSP_u1_columns]
          current_DSP_ub[!DSP_u1_columns] <- 0
          current_DSP_vtype <- DSP_vtype#[DSP_u1_columns]
          
          DSP_variables_per_iteration <- rbind(DSP_variables_per_iteration, c(length(DSP_ub), 
                                                                              sum(!DSP_u1_columns),
                                                                              betas_eliminated))
          
          if(length(initial_DSP_sol) > 0)
          {
            tmp <- numeric(DSP_Nvar)
            tmp[1:length(initial_DSP_sol)] <- initial_DSP_sol
            initial_DSP_sol <- tmp
          }
          
          #SOLVE DSP and obtain (Gamma, eta) and (theta, beta, lambda, nu)
          DSP_result <<- solve_model3Tree_BD_DSP_problem(current_DSP_obj, current_DSP_constraints, DSP_rhs, DSP_senses, 
                                                         current_DSP_lb, current_DSP_ub, current_DSP_vtype, 
                                                         TN, P, alphas, DSP_variable_indexer, DSP_row_indexer, 
                                                         initial_DSP_sol, strengthend_with_min_grouping = TRUE)
          
          initial_DSP_sol <- DSP_result$dual_solution$Sol
          
          eta <- DSP_result$primal_solution$Eta
          Gamma <- DSP_result$primal_solution$Gamma
          CPU_DSP <- CPU_DSP + DSP_result$CPU
          sub_objective <- DSP_result$objective
          
          for(s in 1:TN)
          {
            # while(length(DSP_result$dual_solution$Beta[[s]]) > ncol(betas_track[[s]]))
            # {
            #   betas_track[[s]] <- cbind(betas_track[[s]], numeric(nrow(betas_track[[s]])))
            # }
            # betas_track[[s]] <- rbind(betas_track[[s]], DSP_result$dual_solution$Beta[[s]])
            
            B0 <- 1*(DSP_result$dual_solution$Beta[[s]] < beta_positivity_threshold)
            range_s_val <- get_beta_val_range(s)
            # print(c(length(range_s_val), length(B0)))
            svm_info[P+3, range_s_val] <- svm_info[P+3, range_s_val] + B0
            reset <- which(B0 == 0)
            svm_info[P+3, range_s_val[reset]] <- 0
          }
          
          if(dsp_i == inner_iterations)
          {
            optimality_cuts <- list()
            
            optimality_cuts[[1]] <- list(Theta = DSP_result$dual_solution$Theta,
                                         Lambda = DSP_result$dual_solution$Lambda,
                                         Beta = DSP_result$dual_solution$Beta)
            
            if(iteration > warmup_iterations)
            {
              if(compute_po_cut)
              {
                #SOLVE PO and obtain (theta, lambda, beta) 
                PO_obj <- DSP_obj
                
                for(s in 1:(TN-1))
                {
                  for(q in (s+1):TN)
                  {
                    v <- DSP_variable_indexer$lambda_sq(s,q)
                    PO_obj[v] <- core_point[s,q]
                  }
                }
                
                PO_constraints <- rbind(DSP_constraints, DSP_obj)
                PO_rhs <- c(DSP_rhs, sub_objective-0.00001)
                PO_senses <- c(DSP_senses, "G")
                
                PO_result <- solve_model3Tree_BD_DSP_problem(PO_obj, PO_constraints, PO_rhs, PO_senses, 
                                                             DSP_lb, DSP_ub, DSP_vtype, 
                                                             TN, P, alphas, DSP_variable_indexer, DSP_row_indexer,
                                                             initial_DSP_sol, strengthend_with_min_grouping = TRUE)
                
                PO_eta <- PO_result$primal_solution$Eta
                PO_Gamma <- PO_result$primal_solution$Gamma
                CPU_PO <- PO_result$CPU
                PO_objective <- PO_result$objective
                
                po_cut <- list(Theta = PO_result$dual_solution$Theta,
                               Lambda = PO_result$dual_solution$Lambda,
                               Beta = PO_result$dual_solution$Beta)
                
                if(add_both_normal_and_po_cut)
                {
                  optimality_cuts[[2]] <- po_cut
                }else
                {
                  optimality_cuts[[1]] <- po_cut
                }
              }
              
              for(oc in 1:length(optimality_cuts))
              {
                Theta <- optimality_cuts[[oc]]$Theta
                Lambda <- optimality_cuts[[oc]]$Lambda
                Beta <- optimality_cuts[[oc]]$Beta
                
                ds_constant <- sum(Theta)
                for(s in 1:TN)
                {
                  range_s <- get_beta_val_range(s)
                  if(length(range_s) > 0)
                  {
                    ds_constant <- ds_constant + sum(svm_info[P+1, range_s]*Beta[[s]])
                  }
                }
                
                bd_cut <- numeric(MP_NVar)
                bd_cut[MP_variable_indexer$mu()] <- 1
                for(s in 1:(TN-1)) #yes still correct
                {
                  for(q in (s+1):TN)
                  {
                    bd_cut[MP_variable_indexer$u_sq(s,q)] <- -Lambda[s,q]
                  }
                }
                
                MP_constraints <- rbind(MP_constraints, bd_cut)
                
                MP_rhs <- c(MP_rhs, ds_constant)
                MP_senses <- c(MP_senses, "G")
              }
            }
          }
        }
        
        # SOLVE SVMs and obtain alpha and J
        kappa <- vector("list", TN)
        
        for (s in 1:TN) {
          
          Keta <- calculate_Keta(Km[[s]], eta[[s]])
          
          C <- parameters$C
          if(parameters$normalize_C)
            C <- C / length(y[[s]])
          
          svm_result <<- solve_classification_svm(Keta, y[[s]], C, parameters$epsilon)
          
          CPU_SVMs[s] <- CPU_SVMs[s] + svm_result$CPU
          
          alpha[[s]] <- svm_result$alpha
          b[[s]] <- svm_result$b
          
          alpha_cut <- svm_result$alpha_original
          J[s] <- svm_result$objective_original
          
          #model$alpha is already multiplied by y
          alpha_sum <- t(alpha_cut) %*% y[[s]]
          kappa[[s]] <- calucalte_kappa(alpha_cut, Km[[s]], P)
          
          svm_info <- cbind(svm_info, c(kappa[[s]], alpha_sum, s, 0))
          
          DSP_Nvar <- DSP_Nvar + 1
          DSP_constraints <- cbind(DSP_constraints, numeric(DSP_Nrow))
          DSP_obj <- c(DSP_obj, alpha_sum)
          DSP_constraints[DSP_row_indexer$gamma_s(s), DSP_Nvar] <- 1
          for(m in 1:P)
            DSP_constraints[DSP_row_indexer$eta_sm(s,m), DSP_Nvar] <- kappa[[s]][m]
          DSP_lb <- c(DSP_lb, 0)
          DSP_ub <- c(DSP_ub, Inf)
          DSP_vtype <- c(DSP_vtype, "C")
        }
      }
    }
    
    obj <- sum(J)
    if(obj < UB & cycles == 0)
    {
      UB <- obj
      best_eta <- eta
      best_alpha <- alpha
      best_b <- b
      best_grouping <- grouping
      best_u <- u_sq
    }
    
    LB <- max(Mu, LB) #for numerical stability
    
    Gap <- (UB-LB)/max(abs(LB), parameters$epsilon)
    
    MP_objectives <- c(MP_objectives, obj)
    
    info <- c(obj, UB, LB, Gap, sum(Gamma), sum(J), CPU_MP, CPU_DSP, CPU_PO, sum(CPU_SVMs), cycles)
    iteration_info <- rbind(iteration_info, info)
    
    # if(iteration %% 5 == 1)
    print(sprintf("%d %f %f %f %f%% %f %f %f %f %f %f %d", iteration, obj, UB, LB, Gap*100, sum(Gamma), sum(J), CPU_MP, CPU_DSP, CPU_PO, sum(CPU_SVMs), cycles))
    
    if(Gap <= parameters$optimality_gap & iteration >= min_iterations & cycles == 0)
    {
      terminate <- TRUE
    }
    
    if(parameters$iteration_count > 0 & iteration >= parameters$iteration_count)
    {
      terminate <- TRUE
    }
    
    total_time <- total_time + CPU_MP + sum(CPU_SVMs) + CPU_DSP + CPU_PO
    
    if(parameters$time_limit > 0)
    {
      
      expected_end_time <- total_time + sum(iteration_info[,c("Time(M)", "Time(DSP)", "Time(PO)", "Time(SVM)")])/nrow(iteration_info)
      if(expected_end_time >= parameters$time_limit)
      {
        terminate <- TRUE
        print(sprintf("Time elapsed: %f, terminated at iteration %d due to time limit of %d", total_time, iteration, parameters$time_limit))
      }
    }
    
    iteration <- iteration + 1
  }
  
  for(s in 1:TN)
  {
    best_eta[[s]] <- normalize_eta(best_eta[[s]], parameters$epsilon)
  }
  
  state <- list(alpha = best_alpha, b = best_b, eta = best_eta, grouping = best_grouping, u = best_u, new_subtours = new_subtours, 
                total_time = total_time, MP_objectives = MP_objectives, iteration_info = iteration_info, parameters = parameters, 
                DSP_variables_per_iteration = DSP_variables_per_iteration)
  
  output <- list(state = state)
  
  return(output)
}

find_components <- function(adjacency_matrix)
{
  N <- nrow(adjacency_matrix)
  components <- list()
  visited <- c()
  unvisited <- 1:N
  while(length(visited) < N)
  {
    v <- unvisited[1]
    unvisited <- unvisited[-1]
    visited <- c(visited, v)
    component <- c(v)
    
    while(TRUE)
    {
      to_remove <- c()
      
      for(w in unvisited)
      {
        if(sum(adjacency_matrix[w, component])+sum(adjacency_matrix[component, w]) >= 1)
        {
          component <- c(component, w)
          to_remove <- c(to_remove, w)
          visited <- c(visited, w)
        }
      }
      if(length(to_remove) > 0)
        unvisited <- setdiff(unvisited, to_remove)
      else
        break
    }
    
    components <- c(components, list(component))
  }
  
  return(components)
}
