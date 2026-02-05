valid_auc <- function(fpr, tpr){
  ord <- order(fpr)
  fpr2 <- fpr[ord]
  tpr2 <- tpr[ord]
  
  # sort and remove duplicate FPR values
  idx <- !duplicated(fpr2)
  
  return(MESS::auc(fpr2[idx], tpr2[idx]))
}

# For our proposed method
Dat_auc_PL_var <- function(data,plFit_list) {
  yhat.pl <- plFit_list$yhat.pl
  fhat.pl <- plFit_list$fhat.pl
  beta.pl <- plFit_list$beta.pl
  h_beta <- plFit_list$h_beta
  h_tau0 <- plFit_list$h_tau0
  h_tau1 <- plFit_list$h_tau1
  h_gamma <- plFit_list$h_gamma
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
    Xmat_f0_1 <- cbind(rep(1,nrow(data)),bs(p1_1,knots=prob_y,Boundary.knots = c(0,1)))
    Xmat_f0_2 <- cbind(rep(1,nrow(data)),bs(p1_2,knots=prob_y,Boundary.knots = c(0,1)))
    Xmat_f1_1 <- cbind(rep(1,nrow(data)),bs(p0_1,knots=(1-prob_y),Boundary.knots = c(0,1)))
    Xmat_f1_2 <- cbind(rep(1,nrow(data)),bs(p0_2,knots=(1-prob_y),Boundary.knots = c(0,1)))
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

# For two-stage model (benchmark)
Dat_fullauc <- function(yhat) {
  # calculate AUC 
  c <- c(seq(0, 1, by=0.01)*max(yhat),1)
  denom_tpr <- sum(yhat)
  denom_fpr <- sum(1-yhat)
  num_tpr <- numeric(length(c))
  num_fpr <- numeric(length(c))
  for (j in 1:length(c)) {
    ind <- (yhat>=c[j])
    num_tpr[j] <- sum(yhat[ind==1])
    num_fpr[j] <- sum((1-yhat)[ind==1])
  }
  tpr <- num_tpr/denom_tpr
  fpr <- num_fpr/denom_fpr
  auc <- valid_auc(fpr,tpr)
  tpr_05 <- sum(yhat[yhat>=0.5])/denom_tpr
  fpr_05 <- sum((1-yhat)[yhat>=0.5])/denom_fpr
  
  return(list(auc=auc,tpr_05=tpr_05,fpr_05=fpr_05))
}