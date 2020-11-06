This repository contains R implementation of the algorithms proposed in "Efficient Multitask Multiple Kernel Learning with Application to Cancer Research" (IEEE Trans. Cybern., under review).

* run_forest.R : shows how to solve the MTMKLC problem based on the TCGA cohorts and Hallmark pathways using BDForest
* run_comparative.R : runs various variants of the proposed algorithms (BDF,OL,VL) to compare their algorithmic performances

MTMKLC methods
------------
* classification_helper.R => helper functions
* solve_classification_models_{MODEL}.R (one for each MODEL = forest, ol, or vl) => support vector machine classification solver, and LP/MILP problems asssociated with each model using CPLEX optimization software
* mtmklc_{MODEL}_train.R (one for each MODEL = forest, ol, or vl) => training procedure for MTMKLC problem based on the specified model
* mtmklc_test.R => test procedure for MTMKLC


