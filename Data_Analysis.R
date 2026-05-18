library(MASS)
library(gam)
library(nleqslv)
library(MESS)
library(glmnet)
library(numDeriv)
set.seed(1)

# data
# load("data_PL_fit_R35_PCAM.RData")
data_PL_fit$`RASH YES/NO` <- ifelse(data_PL_fit$`RASH YES/NO`==2,1,0)

# define functions
Dat_adaLasso_crossfit <- function(data) {
  fold.ID<-sample(c(1:nfold),size=nrow(data), replace=TRUE)
  pred <- rep(NA,nrow(data))
  X_gamma <- NA
  
  # storage fit in each fold
  coef_mat <- matrix(NA, nrow = nfold, ncol = dim_X + 1)
  colnames(coef_mat) <- c("(Intercept)", colnames(data)[1:dim_X])
  
  for (i in 1:nfold) {
    train_index <- which(fold.ID!=i)
    test_index <- which(fold.ID==i)
    X_train <- as.matrix(data[train_index,c(1:dim_X)])
    X_test <- as.matrix(data[test_index,c(1:dim_X)])
    Y_train <- data$Y[train_index]
    Y_test <- data$Y[test_index]
    ridge <- cv.glmnet(X_train,Y_train,family="binomial",alpha=0)
    best_ridge_coef <- coef(ridge, s = ridge$lambda.min)
    best_ridge_coef <- as.numeric(best_ridge_coef)[-1]
    count_vars <- c("X..BAND.NEUTROPHILS..count","X..MONOCYTES..count","X..NEUTROPHILS.MANUAL..count","X..LYMPHOCYTES.MANUAL..count",
                    "X..LYMPHOCYTES..count","ALT..count","BILIRUBIN.DIRECT..count","CREATININE..count","HEMATOCRIT..count",
                    "POTASSIUM..count","PROTEIN.TOTAL..count","URIC.ACID..count",'has.CA.125',
                    'has.CANCER.ANTIGEN.15.3','has.CANCER.ANTIGEN.19.9','has.CARBOXYHEMOGLOBIN','has.CARCINOEMBRYONIC.ANT',
                    'has.FERRITIN','has.FIBRINOGEN','has.LACTATE.DEHYDROGENASE','has.LACTATE','has.P.Axis','has.PR.Interval',
                    'has.PT','has.PTT','has.Q.T.INTERVAL','has.QRS.DURATION','has.R.AXIS','has.T.AXIS',
                    'has.THYROID.STIMULATING.HORMONE','has.VENTRICULAR.RATE','blood.lab.counts')
    zero_penalty <- which(colnames(data_PL_fit) %in% count_vars)
    penalty <- 1/abs(best_ridge_coef)
    penalty[zero_penalty] <- 0
    adaptive_lasso <- cv.glmnet(x=X_train,y=Y_train,family="binomial",alpha=1,penalty.factor=penalty)
    yhat.cv <- predict(adaptive_lasso,s="lambda.min",X_test,type="response")
    pred[test_index] <- yhat.cv
    X_gamma[test_index] <- predict(adaptive_lasso,s="lambda.min",X_test)
    coef_i <- as.numeric(coef(adaptive_lasso, s = "lambda.min"))
    coef_mat[i, ] <- coef_i
  }
  
  # mean coefficient
  fit <- colMeans(coef_mat, na.rm = TRUE)
  
  return(list(yhat.full=pred,X_gamma=X_gamma,fit=fit,coef_mat=coef_mat))
}

## Two-Stage I
Dat_PLFit_crossfit <- function(data,X_gamma,yhat.full,fit) {
  coef <- fit
  X_gamma <- X_gamma
  prob_initial_Y1 <- yhat.full
  prob_initial_Y0 <- 1-prob_initial_Y1
  data <- cbind(data,prob_initial_Y1,prob_initial_Y0)
  names(data)[211] <- "prob_initial_Y1"
  names(data)[212] <- "prob_initial_Y0"
  coef_gamma <- coef[which(coef!=0)]
  
  # score function for gamma
  X1 <- as.matrix(cbind(rep(1,nrow(data)),data[1:(dim_X)]))
  X1_A <- X1[,which(coef!=0)]
  dim_X_A <- ncol(X1_A)
  S_gamma <- apply(cbind(data$Y,X1_A, prob_initial_Y1), 1, function(x) x[2:(dim_X_A+1)]*(x[1]-x[dim_X_A+2]))
  
  # Fisher information matrix for gamma
  I_gamma <- solve(t(X1_A) %*% diag(as.vector(prob_initial_Y1*(prob_initial_Y0))) %*% X1_A/N)
  
  # influence function for gamma
  h_gamma <-  I_gamma %*% S_gamma
  
  # create dataset of phase I cases and phase I controls
  data_0 <- data[which(data$Y==0),]
  data_1 <- data[which(data$Y==1),]
  
  # fit missingness model for cases and controls
  f1_fit <- glm(R~bs(prob_initial_Y0,knots = (1-prob_y),Boundary.knots = c(0,1)),data=data_1,family=binomial)
  f0_fit <- glm(R~bs(prob_initial_Y1,knots = prob_y,Boundary.knots = c(0,1)),data=data_0,family=binomial)
  
  # calculate predicted probabilities for offset
  f1_hat <- predict(f1_fit,newdata=data,type="response")
  f0_hat <- predict(f0_fit,newdata=data,type="response")
  
  # calculate X matrix for f1 and f0
  Xmat_f1 <- cbind(rep(1,nrow(data)),bs(data$prob_initial_Y0,knots = (1-prob_y),Boundary.knots = c(0,1)))
  Xmat_f0 <- cbind(rep(1,nrow(data)),bs(data$prob_initial_Y1,knots = prob_y,Boundary.knots = c(0,1)))
  
  # calculate score function for tau1 and tau0
  S_tau1 <- apply(cbind(data$Y,data$R,Xmat_f1, f1_hat), 1, function(x) x[1]*x[3:(3+ncol(Xmat_f1)-1)]*(x[2]-x[(3+ncol(Xmat_f1))]))
  S_tau0 <- apply(cbind((1-data$Y),data$R,Xmat_f0, f0_hat), 1, function(x) x[1]*x[3:(3+ncol(Xmat_f0)-1)]*(x[2]-x[(3+ncol(Xmat_f0))]))
  
  # A matrix for tau (derivative of score with respect to tau)
  A_tau1 <- -solve(t(Xmat_f1) %*% diag(data$Y*as.vector(f1_hat*(1-f1_hat))) %*% (Xmat_f1)/N)
  A_tau0 <- -solve(t(Xmat_f0) %*% diag((1-data$Y)*as.vector(f0_hat*(1-f0_hat))) %*% (Xmat_f0)/N)
  
  # B matrix for tau (derivative of score with respect to gamma) - calculate numerically
  b1gammaforjaco <- function(gammader) {
    p0_jaco <- 1-exp(X1_A%*%gammader)/(1+exp(X1_A%*%gammader))
    
    # create dataset of phase I cases and phase I controls
    data_jaco <- cbind(data,p0_jaco)
    
    #calculate X matrix
    Xmat_f1_jaco <- cbind(rep(1,nrow(data_jaco)),bs(data_jaco$p0_jaco,knots=(1-prob_y),Boundary.knots = c(0,1)))
    f1_hat_jaco <- exp(Xmat_f1_jaco%*%f1_fit$coefficients)/(1+exp(Xmat_f1_jaco%*%f1_fit$coefficients))
    
    big_score <- apply(cbind(data_jaco$Y,data_jaco$R,Xmat_f1_jaco, f1_hat_jaco), 1, function(x) x[1]*x[3:(3+ncol(Xmat_f1_jaco)-1)]*(x[2]-x[(3+ncol(Xmat_f1_jaco))]))
    mean_score <- rowMeans(big_score)
    return(mean_score)
  }
  
  b0gammaforjaco <- function(gammader) {
    p1_jaco <- exp(X1_A%*%gammader)/(1+exp(X1_A%*%gammader))
    
    # create dataset of phase I cases and phase I controls
    data_jaco <- cbind(data,p1_jaco)
    
    #calculate X matrix
    Xmat_f0_jaco <- cbind(rep(1,nrow(data_jaco)),bs(data_jaco$p1_jaco,knots=prob_y,Boundary.knots = c(0,1)))
    f0_hat_jaco <- exp(Xmat_f0_jaco%*%f0_fit$coefficients)/(1+exp(Xmat_f0_jaco%*%f0_fit$coefficients))
    
    big_score <- apply(cbind((1-data_jaco$Y),data_jaco$R,Xmat_f0_jaco, f0_hat_jaco), 1, function(x) x[1]*x[3:(3+ncol(Xmat_f0_jaco)-1)]*(x[2]-x[(3+ncol(Xmat_f0_jaco))]))
    mean_score <- rowMeans(big_score)
    return(mean_score)
  }
  
  gammader <- (coef)[which(coef!=0)]
  B_tau1_gamma <- jacobian(b1gammaforjaco,gammader)
  B_tau0_gamma <- jacobian(b0gammaforjaco,gammader)
  
  # Calculate influence function for tau
  h_tau1 <- A_tau1 %*% (-B_tau1_gamma %*% h_gamma - S_tau1)
  h_tau0 <- A_tau0 %*% (-B_tau0_gamma %*% h_gamma - S_tau0)
  
  # calculate variance for tau
  coef_tau1 <- f1_fit$coefficients
  coef_tau0 <- f0_fit$coefficients
  
  # estimating equation for our two phase method
  XZ <- as.matrix(cbind(rep(1,nrow(data)),as.numeric(X_gamma),data[,(dim_X+1):(dim_X+dim_Z)]))
  XZ[which(is.na(XZ))] <- 0
  estequal <- function(x){
    U <-  colSums(XZ*data$R*as.vector(data$Y - (exp(XZ%*%x)*as.vector(f1_hat/f0_hat))/(1 + exp(XZ%*%x)*as.vector(f1_hat/f0_hat))))
    return(U)
  }
  
  # solve the score function
  inibeta=rep(0,dim_Z+2)
  coef <- nleqslv(inibeta,estequal)$x
  
  # yscore.pl is equivalent to u_i
  yscore.pl <- exp(XZ%*%coef)*as.vector(f1_hat/f0_hat)/(1+exp(XZ%*%coef)*as.vector(f1_hat/f0_hat))
  
  # Calculate U for estimating equation
  U <- apply(cbind(data$R,data$Y,XZ, yscore.pl), 1, function(x) x[1]*x[3:(length(coef)+2)]*(x[2]-x[length(coef)+3]))
  
  # Calculate D_beta_inv
  neg_D_beta_inv <- solve(t(data$R*XZ) %*% diag(as.vector(yscore.pl*(1-yscore.pl))) %*% (XZ)/N)
  
  dgammaforjaco <- function(gammader) {
    p1_jaco <- exp(X1_A%*%gammader)/(1+exp(X1_A%*%gammader))
    p0_jaco <- 1-p1_jaco
    X_gamma_jaco <- X1_A%*%gammader
    
    # create dataset of phase I cases and phase I controls
    data_jaco <- cbind(data,p1_jaco,p0_jaco)
    
    #calculate X matrix
    Xmat_f0_jaco <- cbind(rep(1,nrow(data_jaco)),bs(data_jaco$p1_jaco,knots=prob_y,Boundary.knots = c(0,1)))
    Xmat_f1_jaco <- cbind(rep(1,nrow(data_jaco)),bs(data_jaco$p0_jaco,knots=(1-prob_y),Boundary.knots = c(0,1)))
    f0_hat_jaco <- exp(Xmat_f0_jaco%*%f0_fit$coefficients)/(1+exp(Xmat_f0_jaco%*%f0_fit$coefficients))
    f1_hat_jaco <- exp(Xmat_f1_jaco%*%f1_fit$coefficients)/(1+exp(Xmat_f1_jaco%*%f1_fit$coefficients))
    XZ_jaco <- as.matrix(cbind(rep(1,nrow(data)),as.numeric(X_gamma_jaco),data[,(dim_X+1):(dim_X+dim_Z)]))
    XZ_jaco[which(is.na(XZ_jaco))] <- 0
    yhat_jaco <- exp(XZ_jaco%*%coef)*(as.vector(f1_hat_jaco/f0_hat_jaco))/(1+exp(XZ_jaco%*%coef)*as.vector(f1_hat_jaco/f0_hat_jaco))
    
    #calculate score
    big_U1 <- apply(cbind(data_jaco$R,data_jaco$Y,XZ_jaco,yhat_jaco), 1, function(x) x[1]*x[3:(length(coef)+2)]*(x[2]-x[length(coef)+3]))
    mean_score <- rowMeans(big_U1)
    return(mean_score)
  }
  D_gamma <- jacobian(dgammaforjaco,gammader)
  
  # ystar.pl is equivalent to xi_i
  ystar.pl <- exp(XZ%*%coef)/(1+exp(XZ%*%coef))
  
  # Calculate D_tau1 and D_tau0
  D_tau1 <- t(data$R*XZ) %*% diag(as.vector(-yscore.pl*(1-yscore.pl)*(1-f1_hat))) %*% (Xmat_f1)/N
  D_tau0 <- t(data$R*XZ) %*% diag(as.vector(yscore.pl*(1-yscore.pl)*(1-f0_hat))) %*% (Xmat_f0)/N
  
  # Claculate influence function for beta
  h_beta <- neg_D_beta_inv %*% (U + D_tau1 %*% h_tau1 + D_tau0 %*% h_tau0 + D_gamma %*% h_gamma)
  var_coef <- diag(h_beta %*% t(h_beta)/N^2)
  
  # estimate of f(x,z)
  fhat.pl <- (data$R)/(as.vector(f1_hat)*ystar.pl+as.vector(f0_hat)*(1-ystar.pl))/N
  
  return(list(beta.pl=coef,yhat.pl=ystar.pl,fhat.pl=fhat.pl,h_beta=h_beta,h_tau0=h_tau0,h_tau1=h_tau1,h_gamma=h_gamma,
              f1_hat=f1_hat,f0_hat=f0_hat,Xmat_f1=Xmat_f1,Xmat_f0=Xmat_f0,var=var_coef,
              coef_tau1=coef_tau1,coef_tau0=coef_tau0,XZ=XZ,coef_gamma=coef_gamma,X1_A=X1_A))
}

# Two-Stage II
Dat_PLFit_crossfit_add <- function(data,X_gamma,yhat.full,fit) {
  coef <- fit
  X_gamma <- X_gamma
  prob_initial_Y1 <- yhat.full
  prob_initial_Y0 <- 1-prob_initial_Y1
  data <- cbind(data,prob_initial_Y1,prob_initial_Y0)
  names(data)[211] <- "prob_initial_Y1"
  names(data)[212] <- "prob_initial_Y0"
  coef_gamma <- coef[which(coef!=0)]
  
  # score function for gamma
  X1 <- as.matrix(cbind(rep(1,nrow(data)),data[1:(dim_X)]))
  X1_A <- X1[,which(coef!=0)]
  dim_X_A <- ncol(X1_A)
  S_gamma <- apply(cbind(data$Y,X1_A, prob_initial_Y1), 1, function(x) x[2:(dim_X_A+1)]*(x[1]-x[dim_X_A+2]))
  
  # Fisher information matrix for gamma
  I_gamma <- solve(t(X1_A) %*% diag(as.vector(prob_initial_Y1*(prob_initial_Y0))) %*% X1_A/N)
  
  # influence function for gamma
  h_gamma <-  I_gamma %*% S_gamma
  
  # create dataset of phase I cases and phase I controls
  data_0 <- data[which(data$Y==0),]
  data_1 <- data[which(data$Y==1),]
  
  # fit missingness model for cases and controls
  f1_fit <- gam(R~blood.lab.counts+MEAN.CELLULAR.HEMOGLOBIN.CONCENTRATION..mean+ALT..count+ALBUMIN..min+n_is_Tumor+
                  HEMOGLOBIN..last+X..LYMPHOCYTES..min+n_is_METS+PAT_AGE+ANION.GAP..last+bs(prob_initial_Y0,knots = (1-prob_y),Boundary.knots = c(0,1)),data=data_1,family=binomial)
  f0_fit <- gam(R~blood.lab.counts+MEAN.CELLULAR.HEMOGLOBIN.CONCENTRATION..mean+ALT..count+ALBUMIN..min+n_is_Tumor+
                  HEMOGLOBIN..last+X..LYMPHOCYTES..min+n_is_METS+PAT_AGE+ANION.GAP..last+bs(prob_initial_Y1,knots = prob_y,Boundary.knots = c(0,1)),data=data_0,family=binomial)
  
  # calculate predicted probabilities for offset
  f1_hat <- predict(f1_fit,newdata=data,type="response")
  f0_hat <- predict(f0_fit,newdata=data,type="response")
  
  # calculate X matrix for f1 and f0
  Xmat_f1 <- cbind(rep(1,nrow(data)),data$blood.lab.counts,data$MEAN.CELLULAR.HEMOGLOBIN.CONCENTRATION..mean,data$ALT..count, data$ALBUMIN..min,
                   data$n_is_Tumor,data$HEMOGLOBIN..last,data$X..LYMPHOCYTES..min,data$n_is_METS, data$PAT_AGE,data$ANION.GAP..last,
                   bs(data$prob_initial_Y0,knots = (1-prob_y),Boundary.knots = c(0,1)))
  Xmat_f0 <- cbind(rep(1,nrow(data)),data$blood.lab.counts,data$MEAN.CELLULAR.HEMOGLOBIN.CONCENTRATION..mean,data$ALT..count, data$ALBUMIN..min,
                   data$n_is_Tumor,data$HEMOGLOBIN..last,data$X..LYMPHOCYTES..min,data$n_is_METS, data$PAT_AGE,data$ANION.GAP..last,
                   bs(data$prob_initial_Y1,knots = prob_y,Boundary.knots = c(0,1)))
  
  # calculate score function for tau1 and tau0
  S_tau1 <- apply(cbind(data$Y,data$R,Xmat_f1, f1_hat), 1, function(x) x[1]*x[3:(3+ncol(Xmat_f1)-1)]*(x[2]-x[(3+ncol(Xmat_f1))]))
  S_tau0 <- apply(cbind((1-data$Y),data$R,Xmat_f0, f0_hat), 1, function(x) x[1]*x[3:(3+ncol(Xmat_f0)-1)]*(x[2]-x[(3+ncol(Xmat_f0))]))
  
  # A matrix for tau (derivative of score with respect to tau)
  A_tau1 <- -solve(t(Xmat_f1) %*% diag(data$Y*as.vector(f1_hat*(1-f1_hat))) %*% (Xmat_f1)/N)
  A_tau0 <- -solve(t(Xmat_f0) %*% diag((1-data$Y)*as.vector(f0_hat*(1-f0_hat))) %*% (Xmat_f0)/N)
  
  # B matrix for tau (derivative of score with respect to gamma) - calculate numerically
  b1gammaforjaco <- function(gammader) {
    p0_jaco <- 1-exp(X1_A%*%gammader)/(1+exp(X1_A%*%gammader))
    
    # create dataset of phase I cases and phase I controls
    data_jaco <- cbind(data,p0_jaco)
    
    #calculate X matrix
    Xmat_f1_jaco <- cbind(rep(1,nrow(data_jaco)),data$blood.lab.counts,data$MEAN.CELLULAR.HEMOGLOBIN.CONCENTRATION..mean,data$ALT..count, data$ALBUMIN..min,
                          data$n_is_Tumor,data$HEMOGLOBIN..last,data$X..LYMPHOCYTES..min,data$n_is_METS, data$PAT_AGE,data$ANION.GAP..last,
                          bs(data_jaco$p0_jaco,knots=(1-prob_y),Boundary.knots = c(0,1)))
    f1_hat_jaco <- exp(Xmat_f1_jaco%*%f1_fit$coefficients)/(1+exp(Xmat_f1_jaco%*%f1_fit$coefficients))
    
    big_score <- apply(cbind(data_jaco$Y,data_jaco$R,Xmat_f1_jaco, f1_hat_jaco), 1, function(x) x[1]*x[3:(3+ncol(Xmat_f1_jaco)-1)]*(x[2]-x[(3+ncol(Xmat_f1_jaco))]))
    mean_score <- rowMeans(big_score)
    return(mean_score)
  }
  
  b0gammaforjaco <- function(gammader) {
    p1_jaco <- exp(X1_A%*%gammader)/(1+exp(X1_A%*%gammader))
    
    # create dataset of phase I cases and phase I controls
    data_jaco <- cbind(data,p1_jaco)
    
    #calculate X matrix
    Xmat_f0_jaco <- cbind(rep(1,nrow(data_jaco)),data$blood.lab.counts,data$MEAN.CELLULAR.HEMOGLOBIN.CONCENTRATION..mean,data$ALT..count, data$ALBUMIN..min,
                          data$n_is_Tumor,data$HEMOGLOBIN..last,data$X..LYMPHOCYTES..min,data$n_is_METS, data$PAT_AGE,data$ANION.GAP..last,
                          bs(data_jaco$p1_jaco,knots=prob_y,Boundary.knots = c(0,1)))
    f0_hat_jaco <- exp(Xmat_f0_jaco%*%f0_fit$coefficients)/(1+exp(Xmat_f0_jaco%*%f0_fit$coefficients))
    
    big_score <- apply(cbind((1-data_jaco$Y),data_jaco$R,Xmat_f0_jaco, f0_hat_jaco), 1, function(x) x[1]*x[3:(3+ncol(Xmat_f0_jaco)-1)]*(x[2]-x[(3+ncol(Xmat_f0_jaco))]))
    mean_score <- rowMeans(big_score)
    return(mean_score)
  }
  
  gammader <- (coef)[which(coef!=0)]
  B_tau1_gamma <- jacobian(b1gammaforjaco,gammader)
  B_tau0_gamma <- jacobian(b0gammaforjaco,gammader)
  
  # Calculate influence function for tau
  h_tau1 <- A_tau1 %*% (-B_tau1_gamma %*% h_gamma - S_tau1)
  h_tau0 <- A_tau0 %*% (-B_tau0_gamma %*% h_gamma - S_tau0)
  
  # calculate variance for tau
  coef_tau1 <- f1_fit$coefficients
  coef_tau0 <- f0_fit$coefficients
  
  # estimating equation for our two phase method
  XZ <- as.matrix(cbind(rep(1,nrow(data)),as.numeric(X_gamma),data[,(dim_X+1):(dim_X+dim_Z)]))
  XZ[which(is.na(XZ))] <- 0
  estequal <- function(x){
    U <-  colSums(XZ*data$R*as.vector(data$Y - (exp(XZ%*%x)*as.vector(f1_hat/f0_hat))/(1 + exp(XZ%*%x)*as.vector(f1_hat/f0_hat))))
    return(U)
  }
  
  # solve the score function
  inibeta=rep(0,dim_Z+2)
  coef <- nleqslv(inibeta,estequal)$x
  
  # yscore.pl is equivalent to u_i
  yscore.pl <- exp(XZ%*%coef)*as.vector(f1_hat/f0_hat)/(1+exp(XZ%*%coef)*as.vector(f1_hat/f0_hat))
  
  # Calculate U for estimating equation
  U <- apply(cbind(data$R,data$Y,XZ, yscore.pl), 1, function(x) x[1]*x[3:(length(coef)+2)]*(x[2]-x[length(coef)+3]))
  
  # Calculate D_beta_inv
  neg_D_beta_inv <- solve(t(data$R*XZ) %*% diag(as.vector(yscore.pl*(1-yscore.pl))) %*% (XZ)/N)
  
  dgammaforjaco <- function(gammader) {
    p1_jaco <- exp(X1_A%*%gammader)/(1+exp(X1_A%*%gammader))
    p0_jaco <- 1-p1_jaco
    X_gamma_jaco <- X1_A%*%gammader
    
    # create dataset of phase I cases and phase I controls
    data_jaco <- cbind(data,p1_jaco,p0_jaco)
    
    #calculate X matrix
    Xmat_f0_jaco <- cbind(rep(1,nrow(data_jaco)),data$blood.lab.counts,data$MEAN.CELLULAR.HEMOGLOBIN.CONCENTRATION..mean,data$ALT..count, data$ALBUMIN..min,
                          data$n_is_Tumor,data$HEMOGLOBIN..last,data$X..LYMPHOCYTES..min,data$n_is_METS, data$PAT_AGE,data$ANION.GAP..last,
                          bs(data_jaco$p1_jaco,knots=prob_y,Boundary.knots = c(0,1)))
    Xmat_f1_jaco <- cbind(rep(1,nrow(data_jaco)),data$blood.lab.counts,data$MEAN.CELLULAR.HEMOGLOBIN.CONCENTRATION..mean,data$ALT..count, data$ALBUMIN..min,
                          data$n_is_Tumor,data$HEMOGLOBIN..last,data$X..LYMPHOCYTES..min,data$n_is_METS, data$PAT_AGE,data$ANION.GAP..last,
                          bs(data_jaco$p0_jaco,knots=(1-prob_y),Boundary.knots = c(0,1)))
    f0_hat_jaco <- exp(Xmat_f0_jaco%*%f0_fit$coefficients)/(1+exp(Xmat_f0_jaco%*%f0_fit$coefficients))
    f1_hat_jaco <- exp(Xmat_f1_jaco%*%f1_fit$coefficients)/(1+exp(Xmat_f1_jaco%*%f1_fit$coefficients))
    XZ_jaco <- as.matrix(cbind(rep(1,nrow(data)),as.numeric(X_gamma_jaco),data[,(dim_X+1):(dim_X+dim_Z)]))
    XZ_jaco[which(is.na(XZ_jaco))] <- 0
    yhat_jaco <- exp(XZ_jaco%*%coef)*(as.vector(f1_hat_jaco/f0_hat_jaco))/(1+exp(XZ_jaco%*%coef)*as.vector(f1_hat_jaco/f0_hat_jaco))
    
    #calculate score
    big_U1 <- apply(cbind(data_jaco$R,data_jaco$Y,XZ_jaco,yhat_jaco), 1, function(x) x[1]*x[3:(length(coef)+2)]*(x[2]-x[length(coef)+3]))
    mean_score <- rowMeans(big_U1)
    return(mean_score)
  }
  D_gamma <- jacobian(dgammaforjaco,gammader)
  
  # ystar.pl is equivalent to xi_i
  ystar.pl <- exp(XZ%*%coef)/(1+exp(XZ%*%coef))
  
  # Calculate D_tau1 and D_tau0
  D_tau1 <- t(data$R*XZ) %*% diag(as.vector(-yscore.pl*(1-yscore.pl)*(1-f1_hat))) %*% (Xmat_f1)/N
  D_tau0 <- t(data$R*XZ) %*% diag(as.vector(yscore.pl*(1-yscore.pl)*(1-f0_hat))) %*% (Xmat_f0)/N
  
  # Claculate influence function for beta
  h_beta <- neg_D_beta_inv %*% (U + D_tau1 %*% h_tau1 + D_tau0 %*% h_tau0 + D_gamma %*% h_gamma)
  var_coef <- diag(h_beta %*% t(h_beta)/N^2)
  
  # estimate of f(x,z)
  fhat.pl <- (data$R)/(as.vector(f1_hat)*ystar.pl+as.vector(f0_hat)*(1-ystar.pl))/N
  
  return(list(beta.pl=coef,yhat.pl=ystar.pl,fhat.pl=fhat.pl,h_beta=h_beta,h_tau0=h_tau0,h_tau1=h_tau1,h_gamma=h_gamma,
              f1_hat=f1_hat,f0_hat=f0_hat,Xmat_f1=Xmat_f1,Xmat_f0=Xmat_f0,var=var_coef,
              coef_tau1=coef_tau1,coef_tau0=coef_tau0,XZ=XZ,coef_gamma=coef_gamma,X1_A=X1_A))
}

# Two-Stage CC
Dat_twostage_XZ <- function(data,ada_lasso_fit) {
  # generate initial estimates of coefficients for model between Y and X using phase I samples
  HighD_X <- as.matrix(data[c(1:dim_X)])
  coef <- coef(ada_lasso_fit,s="lambda.min")
  coef_selected <- rownames(coef)[which(as.numeric(coef)!=0)]
  X_gamma <- cbind(rep(1,N),HighD_X)%*%coef
  
  prob_initial_Y1 <- predict(ada_lasso_fit,HighD_X,type="response",s="lambda.min")
  prob_initial_Y0 <- 1-prob_initial_Y1
  
  # score function for gamma
  X1 <- as.matrix(cbind(rep(1,nrow(data)),data[1:(dim_X)]))
  X1_A <- X1[,which(as.numeric(coef)!=0)]
  dim_X_A <- ncol(X1_A)
  S_gamma <- apply(cbind(data$Y,X1_A, prob_initial_Y1), 1, function(x) x[2:(dim_X_A+1)]*(x[1]-x[dim_X_A+2]))
  
  # Fisher information matrix for gamma
  I_gamma <- solve(t(X1_A) %*% diag(as.vector(prob_initial_Y1*(prob_initial_Y0))) %*% X1_A/N)
  
  # influence function for gamma
  h_gamma <-  I_gamma %*% S_gamma
  
  # fit the stage II model
  Y <- data$Y
  X_gamma <- as.numeric(X_gamma)
  Z <- data[,(dim_X+1):(dim_X+dim_Z)]
  tmp <- cbind(Y,X_gamma,Z)
  mod_stageII <- glm(Y~.,family="binomial",data=tmp)
  coef_stageII <- mod_stageII$coefficients
  # yscore.pl is equivalent to u_i
  yscore <- predict(mod_stageII,type="response")
  X_mat <- as.matrix(cbind(rep(1,nrow(data)),X_gamma,Z))
  
  # Calculate U for estimating equation
  U <- apply(cbind(data$Y,X_mat, yscore), 1, function(x) x[2:(length(coef_stageII)+1)]*(x[1]-x[length(coef_stageII)+2]))
  
  # Calculate D_beta_inv
  neg_D_beta_inv <- solve(t(X_mat) %*% diag(as.vector(yscore*(1-yscore))) %*% (X_mat)/N)
  
  gammader <- (coef)[which(as.numeric(coef)!=0)]
  dgammaforjaco <- function(gammader) {
    X_gamma_jaco <- X1_A%*%gammader
    X_mat_jaco <- cbind(rep(1,nrow(data)),X_gamma_jaco,Z)
    yhat_jaco <- exp(as.matrix(X_mat_jaco)%*%coef_stageII)/(1+exp(as.matrix(X_mat_jaco)%*%coef_stageII))
    
    #calculate score
    big_U1 <- apply(cbind(data$Y,X_mat_jaco,yhat_jaco), 1, function(x) x[2:(length(coef_stageII)+1)]*(x[1]-x[length(coef_stageII)+2]))
    mean_score <- rowMeans(big_U1)
    return(mean_score)
  }
  
  D_gamma <- jacobian(dgammaforjaco,gammader)
  
  
  # ystar is equivalent to xi_i
  ystar <- exp(X_mat%*%coef_stageII)/(1+exp(X_mat%*%coef_stageII))
  
  # Claculate influence function for beta
  h_beta <- neg_D_beta_inv %*% (U + D_gamma %*% h_gamma)
  var_coef <- diag(h_beta %*% t(h_beta)/N^2)
  
  return(list(beta=coef_stageII,yhat=ystar,h_beta=h_beta,var=var_coef,X_mat=X_mat,X1_A=X1_A,h_gamma=h_gamma,
              gamma=gammader))
}

valid_auc <- function(fpr, tpr){
  ord <- order(fpr)
  fpr2 <- fpr[ord]
  tpr2 <- tpr[ord]
  idx <- !duplicated(fpr2)
  return(MESS::auc(fpr2[idx], tpr2[idx]))
}

Dat_auc_PL <- function(data,plFit_list) {
  yhat.pl <- plFit_list$yhat.pl
  fhat.pl <- plFit_list$fhat.pl
  beta.pl <- plFit_list$beta.pl
  h_beta <- plFit_list$h_beta
  h_tau0 <- plFit_list$h_tau0
  h_tau1 <- plFit_list$h_tau1
  f1_hat <- as.vector(plFit_list$f1_hat)
  f0_hat <- as.vector(plFit_list$f0_hat)
  tau1.pl <- plFit_list$coef_tau1
  tau0.pl <- plFit_list$coef_tau0
  Xmat_f1 <- plFit_list$Xmat_f1
  Xmat_f0 <- plFit_list$Xmat_f0
  XZ <- plFit_list$XZ
  gamma <- plFit_list$coef_gamma
  h_gamma <- plFit_list$h_gamma
  X1_A <- plFit_list$X1_A
  
  # calculate the AUC
  denom_tpr.pl <- sum(yhat.pl*fhat.pl)
  denom_fpr.pl <- sum((1-yhat.pl)*fhat.pl)
  c <- c(seq(0, 1, by=0.01)*max(yhat.pl*data$R),1)
  num_tpr.pl <- numeric(length(c))
  num_fpr.pl <- numeric(length(c))
  for (j in 1:length(c)){
    ind <- (yhat.pl>=c[j])
    num_tpr.pl[j] <- sum((yhat.pl*fhat.pl)[ind==1])
    num_fpr.pl[j] <- sum(((1-yhat.pl)*fhat.pl)[ind==1])
  }
  tpr.pl <- num_tpr.pl/denom_tpr.pl
  fpr.pl <- num_fpr.pl/denom_fpr.pl
  auc <- valid_auc(fpr.pl, tpr.pl)
  
  # numerical method for dFPR/dbeta and dTPR/dbeta as a function of c
  # Estimate df/dbeta
  df_beta <- apply(cbind(XZ, -fhat.pl^2*N*(f1_hat-f0_hat)*yhat.pl*(1-yhat.pl)), 1, function(x) x[1:(dim_Z+2)]*x[(dim_Z+2)+1])
  
  dFPR_beta <- matrix(0,nrow=length(c),ncol=(dim_Z+2))
  dTPR_beta <- matrix(0,nrow=length(c),ncol=(dim_Z+2))
  yfit1 <- numeric(nrow(data))
  yfit2 <- numeric(nrow(data))
  for (i in 1:(dim_Z+2)) {
    tmp1 <- beta.pl
    tmp2 <- beta.pl
    tmp1[i] <- beta.pl[i] + 0.001
    tmp2[i] <- beta.pl[i] - 0.001
    yfit1 <- exp(XZ %*% tmp1)/(1+exp(XZ %*% tmp1))
    yfit2 <- exp(XZ %*% tmp2)/(1+exp(XZ %*% tmp2))
    fhat1 <- fhat.pl + df_beta[i,]*0.001
    fhat2 <- fhat.pl - df_beta[i,]*0.001
    num1_tpr <- 0
    num2_tpr <- 0
    num1_fpr <- 0
    num2_fpr <- 0
    denom1_tpr <- sum(yfit1*fhat1)
    denom2_tpr <- sum(yfit2*fhat2)
    denom1_fpr <- sum((1-yfit1)*fhat1)
    denom2_fpr <- sum((1-yfit2)*fhat2)
    
    for (j in 1:length(c)){
      ind1 <- (yfit1>=c[j])
      ind2 <- (yfit2>=c[j])
      num1_tpr <- sum((yfit1*fhat1)[ind1==1])
      num2_tpr <- sum((yfit2*fhat2)[ind2==1])
      num1_fpr <- sum(((1-yfit1)*fhat1)[ind1==1])
      num2_fpr <- sum(((1-yfit2)*fhat2)[ind2==1])
      dTPR_beta[j,i] <- (num1_tpr/denom1_tpr-num2_tpr/denom2_tpr)/(2*0.001)
      dFPR_beta[j,i] <- (num1_fpr/denom1_fpr-num2_fpr/denom2_fpr)/(2*0.001)
    }
  }
  numtau1 <- length(tau1.pl)
  epsilon_tau1 <- abs(tau1.pl)/numtau1
  df_tau1 <- apply(cbind(Xmat_f1, -fhat.pl^2*N*f1_hat*(1-f1_hat)*yhat.pl), 1, function(x) x[1:numtau1]*x[numtau1+1])
  dFPR_tau1 <- matrix(0, nrow=length(c), ncol=numtau1)
  dTPR_tau1 <- matrix(0, nrow=length(c), ncol=numtau1)
  
  for (i in 1:numtau1) {
    fhat1 <- fhat.pl + df_tau1[i,]*0.001
    fhat2 <- fhat.pl - df_tau1[i,]*0.001
    num1_tpr <- 0
    num2_tpr <- 0
    num1_fpr <- 0
    num2_fpr <- 0
    denom1_tpr <- sum(yhat.pl*fhat1)
    denom2_tpr <- sum(yhat.pl*fhat2)
    denom1_fpr <- sum((1-yhat.pl)*fhat1)
    denom2_fpr <- sum((1-yhat.pl)*fhat2)
    
    for (j in 1:length(c)){
      ind <- (yhat.pl>=c[j])
      num1_tpr <- sum((yhat.pl*fhat1)[ind==1])
      num2_tpr <- sum((yhat.pl*fhat2)[ind==1])
      num1_fpr <- sum(((1-yhat.pl)*fhat1)[ind==1])
      num2_fpr <- sum(((1-yhat.pl)*fhat2)[ind==1])
      
      dTPR_tau1[j,i] <- (num1_tpr/denom1_tpr-num2_tpr/denom2_tpr)/(2*0.001)
      dFPR_tau1[j,i] <- (num1_fpr/denom1_fpr-num2_fpr/denom2_fpr)/(2*0.001)
    }
  }
  
  numtau0 <- length(tau0.pl)
  epsilon_tau0 <- abs(tau0.pl)/numtau0
  df_tau0 <- apply(cbind(Xmat_f0, -fhat.pl^2*N*f0_hat*(1-f0_hat)*(1-yhat.pl)), 1, function(x) x[1:numtau0]*x[numtau0+1])
  dFPR_tau0 <- matrix(0, nrow=length(c), ncol=numtau0)
  dTPR_tau0 <- matrix(0, nrow=length(c), ncol=numtau0)
  for (i in 1:numtau1) {
    tmp1 <- tau0.pl
    tmp2 <- tau0.pl
    fhat1 <- fhat.pl + df_tau0[i,]*0.001
    fhat2 <- fhat.pl - df_tau0[i,]*0.001
    num1_tpr <- 0
    num2_tpr <- 0
    num1_fpr <- 0
    num2_fpr <- 0
    denom1_tpr <- sum(yhat.pl*fhat1)
    denom2_tpr <- sum(yhat.pl*fhat2)
    denom1_fpr <- sum((1-yhat.pl)*fhat1)
    denom2_fpr <- sum((1-yhat.pl)*fhat2)
    
    for (j in 1:length(c)){
      ind <- (yhat.pl>=c[j])
      num1_tpr <- sum((yhat.pl*fhat1)[ind==1])
      num2_tpr <- sum((yhat.pl*fhat2)[ind==1])
      num1_fpr <- sum(((1-yhat.pl)*fhat1)[ind==1])
      num2_fpr <- sum(((1-yhat.pl)*fhat2)[ind==1])
      
      dTPR_tau0[j,i] <- (num1_tpr/denom1_tpr-num2_tpr/denom2_tpr)/(2*0.001)
      dFPR_tau0[j,i] <- (num1_fpr/denom1_fpr-num2_fpr/denom2_fpr)/(2*0.001)
    }
  }
  
  # calculate the derivative of f with respect to gamma numerically
  df_dgamma <- matrix(0,nrow=nrow(data),ncol=length(gamma))
  for (i in 1:length(gamma)) {
    tmp1 <- gamma
    tmp2 <- gamma
    tmp1[i] <- gamma[i] + 0.001
    tmp2[i] <- gamma[i] - 0.001
    X_gamma1 <- X1_A%*%tmp1
    X_gamma2 <- X1_A%*%tmp2
    X_mat1 <- as.matrix(XZ)
    X_mat1[,2] <- X_gamma1
    X_mat2 <- as.matrix(XZ)
    X_mat2[,2] <- X_gamma2
    yhat1 <- exp(X_mat1 %*% beta.pl)/(1+exp(X_mat1 %*% beta.pl))
    yhat2 <- exp(X_mat2 %*% beta.pl)/(1+exp(X_mat2 %*% beta.pl))
    p1_1 <- exp(X1_A%*%tmp1)/(1+exp(X1_A%*%tmp1))
    p1_2 <- exp(X1_A%*%tmp2)/(1+exp(X1_A%*%tmp2))
    p0_1 <- 1-p1_1
    p0_2 <- 1-p1_2
    Xmat_f0_1 <- cbind(Xmat_f0[,1:(length(tau0.pl)-4)],bs(p1_1,knots=prob_y,Boundary.knots = c(0,1)))
    Xmat_f0_2 <- cbind(Xmat_f0[,1:(length(tau0.pl)-4)],bs(p1_2,knots=prob_y,Boundary.knots = c(0,1)))
    Xmat_f1_1 <- cbind(Xmat_f1[,1:(length(tau1.pl)-4)],bs(p0_1,knots=(1-prob_y),Boundary.knots = c(0,1)))
    Xmat_f1_2 <- cbind(Xmat_f1[,1:(length(tau1.pl)-4)],bs(p0_2,knots=(1-prob_y),Boundary.knots = c(0,1)))
    f0_hat_1 <- exp(Xmat_f0_1%*%tau0.pl)/(1+exp(Xmat_f0_1%*%tau0.pl))
    f0_hat_2 <- exp(Xmat_f0_2%*%tau0.pl)/(1+exp(Xmat_f0_2%*%tau0.pl))
    f1_hat_1 <- exp(Xmat_f1_1%*%tau1.pl)/(1+exp(Xmat_f1_1%*%tau1.pl))
    f1_hat_2 <- exp(Xmat_f1_2%*%tau1.pl)/(1+exp(Xmat_f1_2%*%tau1.pl))
    fhat1 <- data$R/(f1_hat_1*yhat1+f0_hat_1*(1-yhat1))/N
    fhat2 <- data$R/(f1_hat_2*yhat2+f0_hat_2*(1-yhat2))/N
    df_dgamma[,i] <- (fhat1-fhat2)/(2*0.001)
  }
  
  # numerical method for dFPR/dgamma and dTPR/dgamma as a function of c
  dFPR_gamma <- matrix(0,nrow=length(c),ncol=length(gamma))
  dTPR_gamma <- matrix(0,nrow=length(c),ncol=length(gamma))
  for (i in 1:length(gamma)) {
    tmp1 <- gamma
    tmp2 <- gamma
    tmp1[i] <- gamma[i] + 0.001
    tmp2[i] <- gamma[i] - 0.001
    X_gamma1 <- X1_A%*%tmp1
    X_gamma2 <- X1_A%*%tmp2
    X_mat1 <- XZ
    X_mat1[,2] <- X_gamma1
    X_mat2 <- XZ
    X_mat2[,2] <- X_gamma2
    yfit1 <- exp(X_mat1 %*% beta.pl)/(1+exp(X_mat1 %*% beta.pl))
    yfit2 <- exp(X_mat2 %*% beta.pl)/(1+exp(X_mat2 %*% beta.pl))
    fhat1 <- fhat.pl + df_dgamma[,i]*0.001
    fhat2 <- fhat.pl - df_dgamma[,i]*0.001
    num1_tpr <- 0
    num2_tpr <- 0
    num1_fpr <- 0
    num2_fpr <- 0
    denom1_tpr <- sum(yfit1*fhat1)
    denom2_tpr <- sum(yfit2*fhat2)
    denom1_fpr <- sum((1-yfit1)*fhat1)
    denom2_fpr <- sum((1-yfit2)*fhat2)
    
    for (j in 1:length(c)){
      ind1 <- (yfit1>=c[j])
      ind2 <- (yfit2>=c[j])
      num1_tpr <- sum((yfit1*fhat1)[ind1==1])
      num2_tpr <- sum((yfit2*fhat2)[ind2==1])
      num1_fpr <- sum(((1-yfit1)*fhat1)[ind1==1])
      num2_fpr <- sum(((1-yfit2)*fhat2)[ind2==1])
      dTPR_gamma[j,i] <- (num1_tpr/denom1_tpr-num2_tpr/denom2_tpr)/(2*0.001)
      dFPR_gamma[j,i] <- (num1_fpr/denom1_fpr-num2_fpr/denom2_fpr)/(2*0.001)
    }
  }
  
  # calculate influence function for TPR/FPR as a function of c
  mat_tpr_pl <- matrix(0, nrow=nrow(data), ncol=length(c))
  mat_fpr_pl <- matrix(0, nrow=nrow(data), ncol=length(c))
  for (j in 1:length(c)){
    ind <- (yhat.pl>=c[j])
    mat_tpr_pl[,j] <- N*fhat.pl*yhat.pl*(ind-tpr.pl[j])/denom_tpr.pl # some values greater than 1 here
    mat_fpr_pl[,j] <- N*fhat.pl*(1-yhat.pl)*(ind-fpr.pl[j])/denom_fpr.pl
  }
  
  # calculate the influence function for tpr/fpr
  h_tpr <- t(dTPR_beta %*% h_beta + dTPR_tau1 %*% h_tau1 + dTPR_tau0 %*% h_tau0 + dTPR_gamma %*% h_gamma) + mat_tpr_pl
  h_fpr <- t(dFPR_beta %*% h_beta + dFPR_tau1 %*% h_tau1 + dFPR_tau0 %*% h_tau0 + dFPR_gamma %*% h_gamma) + mat_fpr_pl
  
  # calculate influence function for auc
  h_auc <- tryCatch(0-apply(h_fpr, 1, function(x) valid_auc(tpr.pl, x))+apply(h_tpr, 1, function(x) valid_auc(fpr.pl, x)),
                    error=function(e){message <- "Cannot solve";message})
  
  var_auc_pl <- tryCatch(t(h_auc) %*% h_auc/N^2,error=function(e){message <- "Could not solve";message})
  return(list(auc=auc,var=var_auc_pl,h_auc=h_auc))
}

Dat_auc_twostage <- function(data,fit_list) {
  yhat <- fit_list$yhat
  beta <- fit_list$beta
  h_beta <- fit_list$h_beta
  X_mat <- fit_list$X_mat
  gamma <- fit_list$gamma
  h_gamma <- fit_list$h_gamma
  X1_A <- fit_list$X1_A
  
  # calculate the AUC
  denom_tpr <- sum(yhat)
  denom_fpr <- sum((1-yhat))
  c <- c(seq(0, 1, by=0.01)*max(yhat),1)
  num_tpr <- numeric(length(c))
  num_fpr <- numeric(length(c))
  for (j in 1:length(c)){
    ind <- (yhat>=c[j])
    num_tpr[j] <- sum((yhat)[ind==1])
    num_fpr[j] <- sum((1-yhat)[ind==1])
  }
  tpr <- num_tpr/denom_tpr
  fpr <- num_fpr/denom_fpr
  auc <- valid_auc(fpr, tpr)
  
  # numerical method for dFPR/dbeta and dTPR/dbeta as a function of c
  dFPR_beta <- matrix(0,nrow=length(c),ncol=ncol(X_mat))
  dTPR_beta <- matrix(0,nrow=length(c),ncol=ncol(X_mat))
  yfit1 <- numeric(nrow(X_mat))
  yfit2 <- numeric(nrow(X_mat))
  for (i in 1:ncol(X_mat)) {
    tmp1 <- beta
    tmp2 <- beta
    tmp1[i] <- beta[i] + 0.001
    tmp2[i] <- beta[i] - 0.001
    yfit1 <- exp(X_mat %*% tmp1)/(1+exp(X_mat %*% tmp1))
    yfit2 <- exp(X_mat %*% tmp2)/(1+exp(X_mat %*% tmp2))
    num1_tpr <- 0
    num2_tpr <- 0
    num1_fpr <- 0
    num2_fpr <- 0
    denom1_tpr <- sum(yfit1)
    denom2_tpr <- sum(yfit2)
    denom1_fpr <- sum(1-yfit1)
    denom2_fpr <- sum(1-yfit2)
    
    for (j in 1:length(c)){
      ind1 <- (yfit1>=c[j])
      ind2 <- (yfit2>=c[j])
      num1_tpr <- sum((yfit1)[ind1==1])
      num2_tpr <- sum((yfit2)[ind2==1])
      num1_fpr <- sum((1-yfit1)[ind1==1])
      num2_fpr <- sum((1-yfit2)[ind2==1])
      dTPR_beta[j,i] <- (num1_tpr/denom1_tpr-num2_tpr/denom2_tpr)/(2*0.001)
      dFPR_beta[j,i] <- (num1_fpr/denom1_fpr-num2_fpr/denom2_fpr)/(2*0.001)
    }
  }
  
  # numerical method for dFPR/dgamma and dTPR/dgamma as a function of c
  dFPR_gamma <- matrix(0,nrow=length(c),ncol=length(gamma))
  dTPR_gamma <- matrix(0,nrow=length(c),ncol=length(gamma))
  for (i in 1:length(gamma)) {
    tmp1 <- gamma
    tmp2 <- gamma
    tmp1[i] <- gamma[i] + 0.001
    tmp2[i] <- gamma[i] - 0.001
    X_gamma1 <- X1_A%*%tmp1
    X_gamma2 <- X1_A%*%tmp2
    X_mat1 <- X_mat
    X_mat1[,2] <- X_gamma1
    X_mat2 <- X_mat
    X_mat2[,2] <- X_gamma2
    yfit1 <- exp(X_mat1 %*% beta)/(1+exp(X_mat1 %*% beta))
    yfit2 <- exp(X_mat2 %*% beta)/(1+exp(X_mat2 %*% beta))
    num1_tpr <- 0
    num2_tpr <- 0
    num1_fpr <- 0
    num2_fpr <- 0
    denom1_tpr <- sum(yfit1)
    denom2_tpr <- sum(yfit2)
    denom1_fpr <- sum(1-yfit1)
    denom2_fpr <- sum(1-yfit2)
    
    for (j in 1:length(c)){
      ind1 <- (yfit1>=c[j])
      ind2 <- (yfit2>=c[j])
      num1_tpr <- sum((yfit1)[ind1==1])
      num2_tpr <- sum((yfit2)[ind2==1])
      num1_fpr <- sum((1-yfit1)[ind1==1])
      num2_fpr <- sum((1-yfit2)[ind2==1])
      dTPR_gamma[j,i] <- (num1_tpr/denom1_tpr-num2_tpr/denom2_tpr)/(2*0.001)
      dFPR_gamma[j,i] <- (num1_fpr/denom1_fpr-num2_fpr/denom2_fpr)/(2*0.001)
    }
  }
  #
  # calculate influence function for TPR/FPR as a function of c
  mat_tpr_pl <- matrix(0, nrow=nrow(data), ncol=length(c))
  mat_fpr_pl <- matrix(0, nrow=nrow(data), ncol=length(c))
  for (j in 1:length(c)){
    ind <- (yhat>=c[j])
    mat_tpr_pl[,j] <- yhat*(ind-tpr[j])/(denom_tpr/N) # some values greater than 1 here
    mat_fpr_pl[,j] <- (1-yhat)*(ind-fpr[j])/(denom_fpr/N)
  }
  
  # calculate the influence function for tpr/fpr
  h_tpr <- t(dTPR_beta %*% h_beta + dTPR_gamma %*%h_gamma) + mat_tpr_pl
  h_fpr <- t(dFPR_beta %*% h_beta + dTPR_gamma %*%h_gamma) + mat_fpr_pl
  
  
  # calculate influence function for auc
  h_auc <- tryCatch(0-apply(h_fpr, 1, function(x) valid_auc(tpr, x))+apply(h_tpr, 1, function(x) valid_auc(fpr, x)),
                    error=function(e){message <- "Cannot solve";message})
  var_auc <- tryCatch(t(h_auc) %*% h_auc/N^2,error=function(e){message <- "Could not solve";message})
  return(list(auc=auc,var=var_auc,h_auc=h_auc))
}

# Set Model parameters
N <- nrow(data_PL_fit)
m <- 5000
n_II <- sum(data_PL_fit$R)
prob_y <- mean(data_PL_fit$Y)
nfold <- 2

# standardize X matrix
X <- data_PL_fit[,c(1:dim_X)]
X_centered <- apply(X,2,function(x) x-mean(x))
Xs <- apply(X_centered,2,function(x) x/sqrt(sum(x^2)/nrow(X)))
dat_standardized <- data.frame(Xs,data_PL_fit[,(dim_X+1):ncol(data_PL_fit)])

# Fit lasso model between X and R
cv_lasso_XR <- cv.glmnet(
  as.matrix(Xs),
  dat_standardized$R,
  alpha = 1,
  family = "binomial"
)

coef_mat <- coef(cv_lasso_XR, s = "lambda.min")
nonzero_coef <- coef_mat[coef_mat[, 1] != 0, , drop = FALSE]
nonzero_df <- data.frame(
  variable = rownames(nonzero_coef),
  coefficient = as.numeric(nonzero_coef)
)
nonzero_df_sorted <- nonzero_df[
  order(abs(nonzero_df$coefficient), decreasing = TRUE),
]

# variable selection (include 4 PRO Variables)
dat_standardized_reduced <- dat_standardized[,-c(206,207,208,209,210,213,214,215)]
dim_Z <- 4
N <- nrow(dat_standardized_reduced)

# Stage I Model
P1OnlyFit_List_crossfit <- Dat_adaLasso_crossfit(dat_standardized)

#### Two-Stage I Model ###
## betas
plFit_list_crossfit <- Dat_PLFit_crossfit(dat_standardized_reduced,
                                          P1OnlyFit_List_crossfit$X_gamma,
                                          P1OnlyFit_List_crossfit$yhat.full,
                                          P1OnlyFit_List_crossfit$fit)

# calcualte ci for betas
plFit_list_crossfit$beta.pl[6] <- -1*plFit_list_crossfit$beta.pl[6]
ci_lower_crossfit <- plFit_list_crossfit$beta.pl+qnorm(0.025)*sqrt(plFit_list_crossfit$var)
ci_upper_crossfit <- plFit_list_crossfit$beta.pl+qnorm(0.975)*sqrt(plFit_list_crossfit$var)

## AUC
plauc_list_crossfit <- Dat_auc_PL(dat_standardized_reduced,plFit_list_crossfit)

# calcualte ci for aucs
auc_crossfit <- plauc_list_crossfit$auc
var_auc_crossfit <- plauc_list_crossfit$var
ci_lower_crossfit <- auc_crossfit+qnorm(0.025)*sqrt(var_auc_crossfit)
ci_upper_crossfit <- auc_crossfit+qnorm(0.975)*sqrt(var_auc_crossfit)

#### Two-Stage II Model ###
## betas
plFit_list_crossfit_II <- Dat_PLFit_crossfit_add(dat_standardized_reduced,
                                                 P1OnlyFit_List_crossfit$X_gamma,
                                                 P1OnlyFit_List_crossfit$yhat.full,
                                                 P1OnlyFit_List_crossfit$fit)

# calcualte ci for betas
plFit_list_crossfit_II$beta.pl[6] <- -1*plFit_list_crossfit_II$beta.pl[6]
ci_lower_crossfit_II <- plFit_list_crossfit_II$beta.pl+qnorm(0.025)*sqrt(plFit_list_crossfit_II$var)
ci_upper_crossfit_II <- plFit_list_crossfit_II$beta.pl+qnorm(0.975)*sqrt(plFit_list_crossfit_II$var)

## AUC
plauc_list_crossfit_II <- Dat_auc_PL(dat_standardized_reduced,plFit_list_crossfit_II)

# calcualte ci for aucs
auc_crossfit_II <- plauc_list_crossfit_II$auc
var_auc_crossfit_II <- plauc_list_crossfit_II$var
ci_lower_crossfit_II <- auc_crossfit_II+qnorm(0.025)*sqrt(var_auc_crossfit_II)
ci_upper_crossfit_II <- auc_crossfit_II+qnorm(0.975)*sqrt(var_auc_crossfit_II)

#### Complete Case ###
data_PL_fit2 <- dat_standardized_reduced[which(dat_standardized_reduced$R==1),]
HighD_X <- as.matrix(data_PL_fit2[c(1:(dim_X))])
N <- nrow(data_PL_fit2)
Y <- data_PL_fit2$Y
ridge <- cv.glmnet(HighD_X,Y,family="binomial",alpha=0)
best_ridge_coef <- coef(ridge, s = ridge$lambda.min)
best_ridge_coef <- as.numeric(best_ridge_coef)[-1]
count_vars <- c("X..BAND.NEUTROPHILS..count","X..MONOCYTES..count","X..NEUTROPHILS.MANUAL..count","X..LYMPHOCYTES.MANUAL..count",
                "X..LYMPHOCYTES..count","ALT..count","BILIRUBIN.DIRECT..count","CREATININE..count","HEMATOCRIT..count",
                "POTASSIUM..count","PROTEIN.TOTAL..count","URIC.ACID..count",'has.CA.125',
                'has.CANCER.ANTIGEN.15.3','has.CANCER.ANTIGEN.19.9','has.CARBOXYHEMOGLOBIN','has.CARCINOEMBRYONIC.ANT',
                'has.FERRITIN','has.FIBRINOGEN','has.LACTATE.DEHYDROGENASE','has.LACTATE','has.P.Axis','has.PR.Interval',
                'has.PT','has.PTT','has.Q.T.INTERVAL','has.QRS.DURATION','has.R.AXIS','has.T.AXIS',
                'has.THYROID.STIMULATING.HORMONE','has.VENTRICULAR.RATE','blood.lab.counts')
zero_penalty <- which(names(data_PL_fit2)%in%count_vars)
penalty <- 1/abs(best_ridge_coef)
penalty[zero_penalty] <- 0
adaptive_lasso <- cv.glmnet(x=HighD_X,y=Y,family="binomial",alpha=1,penalty.factor=penalty)

## betas
twostage_fit_list <- Dat_twostage_XZ(data_PL_fit2,adaptive_lasso)

# calcualte ci for betas
twostage_fit_list$beta[6] <- -1*twostage_fit_list$beta[6]
ci_lower <- twostage_fit_list$beta+qnorm(0.025)*sqrt(twostage_fit_list$var)
ci_upper <- twostage_fit_list$beta+qnorm(0.975)*sqrt(twostage_fit_list$var)

## AUC
twostage_auc_list <- Dat_auc_twostage(data_PL_fit2,twostage_fit_list)

# calcualte ci for aucs
auc_res <- twostage_auc_list$auc
auc_var <- twostage_auc_list$var
ci_lower <- auc_res+qnorm(0.025)*sqrt(auc_var)
ci_upper <- auc_res+qnorm(0.975)*sqrt(auc_var)
