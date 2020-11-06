# Mehmet Gonen (mehmetgonen@ku.edu.tr)

library(AUC)

args <- commandArgs(trailingOnly = TRUE)
K <- as.numeric(args[[1]])

source(sprintf("solve_classification_models_%s.R", "forest"))
source(sprintf("mtmklc_%s_train.R", "forest"))
source("mtmklc_test.R")
source("classification_helper.R")

cohorts <- c("TCGA-BRCA", "TCGA-COAD", "TCGA-ESCA", "TCGA-HNSC", "TCGA-KICH",
             "TCGA-KIRC", "TCGA-KIRP", "TCGA-LIHC", "TCGA-LUAD", "TCGA-LUSC",
             "TCGA-PAAD", "TCGA-READ", "TCGA-STAD", "TCGA-TGCT", "TCGA-THCA")

TN <- length(cohorts)

data_path <- "../data"

pathway <- "hallmark"


normalize_C <- TRUE #If true, C is devided by the number of data points of each task

result_path <- sprintf("results_bdforest_K%d", K)
if (dir.exists(sprintf("%s", result_path)) == FALSE) {
  dir.create(sprintf("%s", result_path)) 
}

settings_forest <- list()
settings_forest$compute_po_cut <- TRUE
settings_forest$core_point <- NA
if(settings_forest$compute_po_cut)
{
  settings_forest$core_point <- get_core_point(TN, K, 0.5)
}

settings_forest$optimality_cut_frequency <- 1
settings_forest$add_both_normal_and_po_cut <- FALSE
settings_forest$warmup_iterations <- 5

settings_forest$beta_positivity_threshold <- 1e-5
settings_forest$beta_steps_before_elimination <- 60

settings_forest$two_sided_eta_diff <- TRUE
settings_forest$cumulative_eta_diff <- TRUE
settings_forest$stregthened_grouping <- TRUE
settings_forest$symmetry_breaking <- TRUE

epsilon <- 1e-5
fold_count <- 4
train_ratio <- 0.8

iteration_count <- 300
optimality_gap <- 0.01 #for cross-validation; will be 1e-4 in the final training
time_limit <- 0 #0 means no limit. greater than zero is in seconds (e.g., 12*3600)

replications <- 1:100
state_files <- sprintf("%s/forest_%s_replication_%d_state.RData", result_path, pathway, replications)

if(all(sapply(state_files, function(sf) all(sapply(sf, file.exists)))) == FALSE)
{
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
  }
  
  pathways <- read_pathways(pathway)
  gene_names <- sort(unique(unlist(sapply(1:length(pathways), FUN = function(x) {pathways[[x]]$symbols}))))
  for (t in 1:TN) {
    X[[t]] <- X[[t]][, which(colnames(X[[t]]) %in% gene_names)]
  }
  
  P <- length(pathways)
  
  tuples <- get_cross_validation_tuples()
  
  for(r in 1:length(replications)) {
    
    parameters <- list()
    parameters$epsilon <- epsilon
    parameters$K <- K
    parameters$iteration_count <- iteration_count
    parameters$time_limit <- time_limit
    parameters$optimality_gap <- optimality_gap
    parameters$normalize_C <- normalize_C
    
    state_file <- state_files[r]
    replication <- replications[r]
    
    if(file.exists(state_file) == FALSE){
      random_seed <- 1505 * replication
      
      train_negative_indices <- vector("list", TN)
      train_positive_indices <- vector("list", TN)
      negative_allocation <- vector("list", TN)
      positive_allocation <- vector("list", TN)
      for (t in 1:TN) {
        set.seed(random_seed)
        train_negative_indices[[t]] <- sample(negative_indices[[t]], ceiling(train_ratio * length(negative_indices[[t]])))
        train_positive_indices[[t]] <- sample(positive_indices[[t]], ceiling(train_ratio * length(positive_indices[[t]])))
        
        negative_allocation[[t]] <- sample(rep(1:fold_count, ceiling(length(train_negative_indices[[t]]) / fold_count)), length(train_negative_indices[[t]]))
        positive_allocation[[t]] <- sample(rep(1:fold_count, ceiling(length(train_positive_indices[[t]]) / fold_count)), length(train_positive_indices[[t]]))
      }
      
      auroc_tuples <- matrix(0, nrow = nrow(tuples), ncol = ncol(tuples)+fold_count)
      colnames(auroc_tuples) <- c("C", paste("Auroc",1:fold_count))
      auroc_tuples[,"C"] <- tuples[,"C"]
      
      time_tuples <- matrix(0, nrow = nrow(tuples), ncol = ncol(tuples)+fold_count)
      colnames(time_tuples) <- c("C", paste("Time",1:fold_count))
      
      K_train <- vector("list", TN)
      K_test <- vector("list", TN)
      y_train <- vector("list", TN)
      y_test <- vector("list", TN)
      
      cv_process_time <- 0
      
      for (fold in 1:fold_count) {
        
        cv_process_start <- Sys.time()
        
        for (t in 1:TN) {
          train_indices <- c(train_negative_indices[[t]][which(negative_allocation[[t]] != fold)], train_positive_indices[[t]][which(positive_allocation[[t]] != fold)])
          test_indices <- c(train_negative_indices[[t]][which(negative_allocation[[t]] == fold)], train_positive_indices[[t]][which(positive_allocation[[t]] == fold)]) 
          
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
        
        cv_process_time <- cv_process_time + as.numeric(Sys.time()-cv_process_start, units="secs")
      
        for(tpl in 1:nrow(tuples))
        {
          C <- tuples[tpl, "C"]
          parameters$C <- C
          
          if(settings_forest$stregthened_grouping)
          {
            output <- mtmklc_forest_bdf_strengthened_train(K_train, y_train, parameters, settings = settings_forest)
          }else
          {
            output <- mtmklc_forest_bdf_cumulative_train(K_train, y_train, parameters, settings = settings_forest)
          }
          
          state <- output$state
          eta <- state$eta
          prediction <- mtmklc_test(K_test, eta, state$alpha, state$b)
          auroc <- numeric(TN)
          for (t in 1:TN) {
            auroc[t] <- auc(roc(prediction$f[[t]], as.factor(y_test[[t]])))
          } 
          result <- list(AUROC = auroc, eta = eta, prediction = prediction, y_test = y_test, time = state$total_time)
          
          auroc_tuples[tpl, paste("Auroc",fold)] <- mean(result$AUROC)
          time_tuples[tpl, paste("Time",fold)] <- mean(result$time)
        }
      }
      
      cv_process_time <- cv_process_time/fold_count
      
      average_aurocs <- rowMeans(auroc_tuples[,ncol(tuples)+(1:fold_count)])
      tuple_star <- which(average_aurocs == max(average_aurocs))[1]
      C <- auroc_tuples[tuple_star, "C"]
      
      train_process_start <- Sys.time()
      
      for (t in 1:TN) {
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
      
      train_process_time <- as.numeric(Sys.time()-train_process_start, units="secs")
      
      print(sprintf("Replication %d => running final training with C = %g", replication, C))
      
      parameters$iteration_count <- iteration_count*2
      parameters$optimality_gap <- optimality_gap*0.01
      parameters$C <- C
      
      if(settings_forest$stregthened_grouping)
      {
        output <- mtmklc_forest_bdf_strengthened_train(K_train, y_train, parameters, settings = settings_forest)
      }else
      {
        output <- mtmklc_forest_bdf_cumulative_train(K_train, y_train, parameters, settings = settings_forest)
      }
      
      state <- output$state
      state$cv_process_time <- cv_process_time
      state$train_process_time <- train_process_time
      state$train_time <- state$total_time
      state$time_tuples <- time_tuples
      state$auroc_tuples <- auroc_tuples
      eta <- state$eta
      prediction <- group_lasso_multiple_kernel_classification_model3Tree_test(K_test, eta, state$alpha, state$b)
      auroc <- numeric(TN)
      for (t in 1:TN) {
        auroc[t] <- auc(roc(prediction$f[[t]], as.factor(y_test[[t]])))
      }
      state$AUROC <- auroc
      state$prediction <- prediction
      state$y_test <- y_test
      
      save("state", file = state_file)
    }
  }
}
