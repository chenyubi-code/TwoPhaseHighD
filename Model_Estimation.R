Dat_adaLasso_crossfit <- function(data) {
  fold.ID <- sample(c(1:nfold),size=nrow(data), replace=TRUE)
  pred <- rep(NA,nrow(data))
  X_gamma <- NA
  
  # storage fit from each fold
  coef_mat <- matrix(NA, nrow = nfold, ncol = dim_X + 1)
  colnames(coef_mat) <- c("(Intercept)", colnames(data)[1:dim_X])
  
  # Cross-Fitting Loop
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
    adaptive_lasso <- cv.glmnet(x=X_train,y=Y_train,family="binomial",alpha=1,penalty.factor=1/abs(best_ridge_coef))
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

Dat_PLFit <- function(data,X_gamma,yhat.full,fit) {
  coef <- fit
  X_gamma <- X_gamma
  prob_initial_Y1 <- yhat.full
  prob_initial_Y0 <- 1-prob_initial_Y1
  data <- cbind(data,prob_initial_Y1,prob_initial_Y0)
  n <- ncol(data)
  names(data)[n-1] <- "prob_initial_Y1"
  names(data)[n]   <- "prob_initial_Y0"
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

Dat_twostage_XZ <- function(data,X_gamma) {
  
  # fit the stage II model
  Y <- data$Y
  X_gamma <- as.numeric(X_gamma)
  Z <- data[,(dim_X+1):(dim_X+dim_Z)]
  tmp <- cbind(Y,X_gamma,Z)
  mod_stageII <- glm(Y~.,family="binomial",data=tmp)
  coef_stageII <- mod_stageII$coefficients
  X_mat <- as.matrix(cbind(rep(1,nrow(data)),X_gamma,Z))
  
  # ystar is equivalent to xi_i
  ystar <- exp(X_mat%*%coef_stageII)/(1+exp(X_mat%*%coef_stageII))
  
  return(list(beta=coef_stageII,yhat=ystar))
}