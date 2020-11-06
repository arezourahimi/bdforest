# Mehmet Gonen (mehmetgonen@ku.edu.tr)

library(AUC)

args <- commandArgs(trailingOnly = TRUE)
TN <- as.numeric(args[[1]])
K <- as.numeric(args[[2]])

data_path <- "../data"

pathway <- "hallmark"
normalize_C <- TRUE #If true, C is devided by the number of data points of each task

for(m in c("forest", "vl", "ol"))
{
  source(sprintf("mtmklc_%s_train.R", m))
  source(sprintf("solve_classification_models_%s.R", m))
}
source("mtmklc_test.R")
source("classification_helper.R")

models <- c("BDFC", "BDFS", "BDFSPO","BDFCPO", "VL", "OLJ1", "OLJ2")
# BDFS or BDFSPO: BDForest with strengthened formulation (PO: Pareto-optimal cuts)
# BDFC or BDFCPO: BDForest with cumulative formulation (PO: Pareto-optimal cuts)
# OL: Objective linearization
# VL: Variable linearization

C_range <- get_cross_validation_tuples()
C_range <- C_range[seq(2,length(C_range),2)]

iteration_count <- -1 #no limit on the number of iterations
optimality_gap <- 1e-4
time_limit <- 3600*12 #seconds (12 hr)
pathway <- "hallmark"

#sorted in decreasing order of number of data points
cohorts <- c("TCGA-BRCA","TCGA-KIRC","TCGA-LUAD","TCGA-THCA","TCGA-LUSC",
             "TCGA-COAD","TCGA-HNSC","TCGA-STAD","TCGA-LIHC","TCGA-KIRP",
             "TCGA-PAAD","TCGA-READ","TCGA-ESCA","TCGA-TGCT","TCGA-KICH")

cohorts <- cohorts[1:TN]

data_path <- "../data"

settings_forest <- list()
settings_forest$optimality_cut_frequency <- 1
settings_forest$add_both_normal_and_po_cut <- FALSE
settings_forest$warmup_iterations <- 5

settings_forest$beta_positivity_threshold <- 1e-5
settings_forest$beta_steps_before_elimination <- 150

settings_forest$two_sided_eta_diff <- TRUE
settings_forest$cumulative_eta_diff <- TRUE
settings_forest$stregthened_grouping <- TRUE
settings_forest$symmetry_breaking <- TRUE

settings_OL <- list()
settings_OL$multi_cut <- TRUE
settings_OL$force_nonempty_groups <- TRUE
settings_OL$J_up <- NULL

settings_VL <- list()
settings_VL$force_nonempty_groups <- TRUE

result_path <- sprintf("./comparative_N%dK%d", TN, K)

if (dir.exists(sprintf("%s", result_path)) == FALSE) {
  dir.create(sprintf("%s", result_path)) 
}

epsilon <- 1e-5
fold_count <- 4
train_ratio <- 0.8

parameters <- list()
parameters$epsilon <- epsilon
parameters$K <- K
parameters$iteration_count <- iteration_count
parameters$time_limit <- time_limit
parameters$optimality_gap <- optimality_gap
parameters$normalize_C <- normalize_C


pathways <- read_pathways(pathway)
gene_names <- sort(unique(unlist(sapply(1:length(pathways), FUN = function(x) {pathways[[x]]$symbols}))))
P <- length(pathways)

X <- vector("list", TN)
y <- vector("list", TN)
negative_indices <- vector("list", TN)
positive_indices <- vector("list", TN)

for (t in 1:TN) {
  load(sprintf("%s/%s.RData", data_path, cohorts[t]))
  
  common_patients <- intersect(rownames(TCGA$clinical)[which(is.na(TCGA$clinical$pathologic_stage) == FALSE)], rownames(TCGA$mrna))
  
  X[[t]] <- log2(TCGA$mrna[common_patients,] + 1)
  y[[t]] <- rep(NA, length(common_patients))
  
  y[[t]][TCGA$clinical[common_patients, "pathologic_stage"] %in% c("Stage I",  "Stage IA",  "Stage IB",  "Stage IC")] <- +1
  y[[t]][TCGA$clinical[common_patients, "pathologic_stage"] %in% c("Stage II", "Stage IIA", "Stage IIB", "Stage IIC",
                                                                   "Stage III", "Stage IIIA", "Stage IIIB", "Stage IIIC",
                                                                   "Stage IV",  "Stage IVA",  "Stage IVB",  "Stage IVC")] <- -1
  
  valid_patients <- which(is.na(y[[t]]) == FALSE)
  valid_features <- as.numeric(which(apply(X[[t]][valid_patients,], 2, sd) != 0))
  X[[t]] <- X[[t]][valid_patients, valid_features]
  y[[t]] <- y[[t]][valid_patients]
  
  negative_indices[[t]] <- which(y[[t]] == -1)
  positive_indices[[t]] <- which(y[[t]] == +1)

  X[[t]] <- X[[t]][, which(colnames(X[[t]]) %in% gene_names)]
}

for(replication in 1:50) {
  
  random_seed <- 1505 * replication
  
  train_negative_indices <- vector("list", TN)
  train_positive_indices <- vector("list", TN)
  negative_allocation <- vector("list", TN)
  positive_allocation <- vector("list", TN)
  K_train <- vector("list", TN)
  K_test <- vector("list", TN)
  y_train <- vector("list", TN)
  y_test <- vector("list", TN)
  for (t in 1:TN) {
    set.seed(random_seed)
    train_negative_indices[[t]] <- sample(negative_indices[[t]], ceiling(train_ratio * length(negative_indices[[t]])))
    train_positive_indices[[t]] <- sample(positive_indices[[t]], ceiling(train_ratio * length(positive_indices[[t]])))
    
    negative_allocation[[t]] <- sample(rep(1:fold_count, ceiling(length(train_negative_indices[[t]]) / fold_count)), length(train_negative_indices[[t]]))
    positive_allocation[[t]] <- sample(rep(1:fold_count, ceiling(length(train_positive_indices[[t]]) / fold_count)), length(train_positive_indices[[t]]))
 
    train_indices <- c(train_negative_indices[[t]], train_positive_indices[[t]])
    test_indices <- setdiff(1:length(y[[t]]), train_indices)
    
    X_train <- X[[t]][train_indices,]
    X_test <- X[[t]][test_indices,]
    X_train <- scale(X_train)
    X_test <- (X_test - matrix(attr(X_train, "scaled:center"), nrow = nrow(X_test), ncol = ncol(X_test), byrow = TRUE)) / matrix(attr(X_train, "scaled:scale"), nrow = nrow(X_test), ncol = ncol(X_test), byrow = TRUE)
    
    N_train <- nrow(X_train)
    N_test <- nrow(X_test)
    K_train[[t]] <- array(0, dim = c(N_train, N_train, P))
    K_test[[t]] <- array(0, dim = c(N_test, N_train, P))
    for (m in 1:P) {
      feature_indices <- which(colnames(X_train) %in% pathways[[m]]$symbols)
      D_train <- pdist(X_train[, feature_indices], X_train[, feature_indices])
      D_test <- pdist(X_test[, feature_indices], X_train[, feature_indices])
      sigma <- mean(D_train)
      K_train[[t]][,,m] <- exp(-D_train^2 / (2 * sigma^2))
      K_test[[t]][,,m] <- exp(-D_test^2 / (2 * sigma^2))
    }
    
    y_train[[t]] <- y[[t]][train_indices]
    y_test[[t]] <- y[[t]][test_indices] 
  }
  
  for(C in C_range)
  {
    parameters$C <- C
    
    for(m in models)
    {
      state_file <- sprintf("%s/%s_rep_%d_C%g_state.RData", result_path, m, replication, C)
      
      if (file.exists(state_file) == FALSE) {
        print(sprintf("Replication %d ==> running model %s with C = %g", replication, m, C))
        
        if(startsWith(m, "BDF"))
        {
          settings_forest$compute_po_cut <- endsWith(m, "PO")
          settings_forest$core_point <- NA
          if(settings_forest$compute_po_cut)
          {
            settings_forest$core_point <- get_core_point(TN, K, 0.5)
          }
          
          if(startsWith(m, "BDFS"))
          {
            output <- mtmklc_forest_bdf_strengthened_train(K_train, y_train, parameters, settings = settings_forest)
          }else
          {
            output <- mtmklc_forest_bdf_cumulative_train(K_train, y_train, parameters, settings = settings_forest)
          }
          
          state <- output$state
          eta <- state$eta
        }else if(m == "VL")
        {
          output <- mtmklc_vl_train(K_train, y_train, parameters, settings = settings_VL)
          state <- output$state
          eta <- lapply(state$grouping, function(l){state$eta[[l]]})
        }else if(startsWith(m, "OL"))
        {
          settings_OL$J_method <- substr(m,3,4)
          
          output <- mtmklc_ol_train(K_train, y_train, parameters, settings = settings_OL)
          state <- output$state
          eta <- lapply(state$grouping, function(l){state$eta[[l]]})
        }
        
        prediction <- mtmklc_test(K_test, eta, state$alpha, state$b)
        
        auroc <- numeric(TN)
        for (t in 1:TN) {
          auroc[t] <- auc(roc(prediction$f[[t]], as.factor(y_test[[t]])))
        }
        state$AUROC <- auroc
        
        save("state", file = state_file)
      }
    }
  }
}
