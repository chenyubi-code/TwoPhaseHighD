### Scenario: SRS, P(R=1)=0.05, rho=0.3 ####################################################################################

library(MASS)
library(gam)
library(nleqslv)
library(MESS)
library(glmnet)
library(numDeriv)
source("Format_Data.R")
source("Model_Estimation.R")
source("Model_Evaluation.R")

###################################################  Set Parameters ######################################################### 

# Total number of individuals in Phase I
N <- 10000
# Prevalence of y
prob_y <- 0.2
# Number of Phase I predictors, X
dim_X <- 300
# Number of Phase I non-zero predictors, X
nonzero_X <- 10
# Dimension of Phase II predictors, Z
dim_Z <- 3
# proportion of individuals in phase II sample
prob_n_II <- 0.05
# rho in covariance matrix
rho <- 0.3
# set log odds ratios for the nonzero X predictors
nonzeroX_logOR <- c(log(2),-log(1.8),log(1.6),-log(1.4),log(1.2),log(2),-log(1.8),log(1.6),-log(1.4),log(1.2))
# set log odds ratios for phase II predictors, Z
Z_logOR <- c(log(1.5),log(2),-log(1.8))
# number of beta coefficients in prediction model
numbeta <- length(nonzeroX_logOR)+length(Z_logOR)+1
# set log odds ratios for all X predictors
X_logOR <- NULL
for (i in 1:nonzero_X) {
  X_logOR <- c(X_logOR,nonzeroX_logOR[i],rep(log(1),dim_X/nonzero_X-1))
}
# log odds ratios for all predictors
gamma_logOR <- c(X_logOR,Z_logOR)
# number of simulation iterations
m <- 1000
# number of folds for cross-fitting
nfold <- 2
# set the intercepts for Phase I and II model
gamma_0 <- -2.102

###################################################  Fit Models ######################################################### 

results_list <- vector("list", m)

for (iter in 1:m) {
  cat("Running iteration:", iter, "\n")
  
  set.seed(iter)
  
  # generate full data
  data.full1 <- data_gen(N, prob_y, nonzero_X, Z_logOR, rho, gamma_0)
  R <- Dat_gensrs(prob_n_II,N)
  data.full <- cbind(data.full1, R)
  
  # fit the phase I model using X that are available for all 
  P1OnlyFit_List <- Dat_adaLasso_crossfit(data.full)
  
  # fit the phase II model and calculate the auc using spline term only to model missingness
  plFit_list <- Dat_PLFit(data.full, P1OnlyFit_List$X_gamma, P1OnlyFit_List$yhat.full, P1OnlyFit_List$fit)
  plauc_list <- Dat_auc_PL_var(data.full, plFit_list)
  
  # Benchmark: fit the two stage model with no missing data
  two_stage_fit_list <- Dat_twostage_XZ(data.full, P1OnlyFit_List$X_gamma)
  two_stage_auc_list <- Dat_fullauc(two_stage_fit_list$yhat)
  
  results_list[[iter]] <- list(
    # Our Method Model
    coef = plFit_list$beta.pl,
    var_coef = plFit_list$var,
    auc = plauc_list$auc,
    var_auc = plauc_list$var,
    # Benchmark: Two-stage Model
    two_stage_coef = two_stage_fit_list$beta,
    two_stage_auc = two_stage_auc_list$auc
  )
}

Coef <- t(sapply(results_list, function(x) x$coef))
Var_coef <- t(sapply(results_list, function(x) x$var_coef))
Auc <- t(sapply(results_list, function(x) x$auc))
Var_auc <- t(sapply(results_list, function(x) x$var_auc))
Two_stage_coef <- t(sapply(results_list, function(x) x$two_stage_coef))
Two_stage_auc <- t(sapply(results_list, function(x) x$two_stage_auc))

simulation_results <- list(
  Coef = Coef,
  Var_coef = Var_coef,
  Auc = Auc,
  Var_auc = Var_auc,
  Two_stage_coef = Two_stage_coef,
  Two_stage_auc = Two_stage_auc
)
