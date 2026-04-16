### Scenario: Single-index, P(R=1)=0.05, rho=0 ####################################################################################

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
# set sparsity for model for R
q=ceiling(sqrt(dim_X))
# proportion of individuals in phase II sample
prob_n_II <- 0.05
# rho in covariance matrix
rho <- 0
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
gamma_0 <- -2.124
tao_00 <- -7.337
tao_01 <- -7.71
# superlearner libraries
SL.library <- c(
  "SL.glm",
  "SL.glmnet",
  "SL.randomForest",
  "SL.rpart",      
  "SL.nnet"         
)

###################################################  Fit Models ######################################################### 

results_list <- vector("list", m)

for (iter in 1:m) {
  cat("Running iteration:", iter, "\n")
  
  set.seed(iter)
  
  # generate full data
  data.full <- data_gen(N, prob_y, nonzero_X, Z_logOR, rho, gamma_0)
  
  # generate missingness indicator R
  data.full.y0 <- data.full[data.full$Y == 0,]
  data.full.y1 <- data.full[data.full$Y == 1,]
  R_Y0 <- Dat_gensingle_index(data.full.y0[, c(1:dim_X)], tao_00)
  R_Y1 <- Dat_gensingle_index(data.full.y1[, c(1:dim_X)], tao_01)
  R <- R_Y1
  data.full.y1 <- cbind(data.full.y1, R)
  R <- R_Y0
  data.full.y0 <- cbind(data.full.y0, R)
  data.full <- rbind(data.full.y1, data.full.y0)
  
  # fit the phase I model using X that are available for all 
  P1OnlyFit_List <- Dat_adaLasso_crossfit(data.full)
  
  # fit the phase II model and calculate the auc using spline term only to model missingness
  plFit_list <- Dat_PLFit(data.full, P1OnlyFit_List$X_gamma, P1OnlyFit_List$yhat.full, P1OnlyFit_List$fit)
  plauc_list <- Dat_auc_PL_var(data.full, plFit_list)
  
  # Benchmark: fit the two stage model with no missing data
  two_stage_fit_list <- Dat_twostage_XZ(data.full, P1OnlyFit_List$X_gamma)
  two_stage_auc_list <- Dat_fullauc(two_stage_fit_list$yhat)
  
  # CC: fit the two stage model for individuals with Z available
  data.P2 <- data.full[which(data.full$R == 1),]
  P2OnlyFit_List <- Dat_adaLasso_crossfit(data.P2)
  P2Fit_list <- Dat_twostage_CC(data.P2, P2OnlyFit_List)
  P2auc_list <- Dat_fullauc(P2Fit_list$yhat)
  
  # Super Learner
  P1_SL_List <- SL_StageI(data.full, SL.library)
  plFit_SL <- SL_StageII(data.full, P1_SL_List)
  plauc_SL <- Dat_auc_SL(data.full, plFit_SL)
  
  results_list[[iter]] <- list(
    # Our Method Model
    coef = plFit_list$beta.pl,
    var_coef = plFit_list$var,
    auc = plauc_list$auc,
    var_auc = plauc_list$var,
    # Benchmark: Two-stage Model
    two_stage_coef = two_stage_fit_list$beta,
    two_stage_auc = two_stage_auc_list$auc,
    # CC
    cc_coef = P2Fit_list$beta,
    cc_auc = P2auc_list$auc,
    # SuperLearner
    sl_coef = plFit_SL$beta.pl,
    sl_auc = plauc_SL$auc
  )
}

Coef <- t(sapply(results_list, function(x) x$coef))
Var_coef <- t(sapply(results_list, function(x) x$var_coef))
Auc <- t(sapply(results_list, function(x) x$auc))
Var_auc <- t(sapply(results_list, function(x) x$var_auc))
Two_stage_coef <- t(sapply(results_list, function(x) x$two_stage_coef))
Two_stage_auc <- t(sapply(results_list, function(x) x$two_stage_auc))
CC_coef <- t(sapply(results_list, function(x) x$cc_coef))
CC_auc <- t(sapply(results_list, function(x) x$cc_auc))
SL_coef <- t(sapply(results_list, function(x) x$sl_coef))
SL_auc <- t(sapply(results_list, function(x) x$sl_auc))

simulation_results <- list(
  Coef = Coef,
  Var_coef = Var_coef,
  Auc = Auc,
  Var_auc = Var_auc,
  Two_stage_coef = Two_stage_coef,
  Two_stage_auc = Two_stage_auc,
  CC_coef = CC_coef,
  CC_auc = CC_auc,
  SL_coef = SL_coef,
  SL_auc = SL_auc
)
