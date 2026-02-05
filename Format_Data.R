func_A1gen <- function(rho,p){
  A1=matrix(0,p,p)
  for(i in 1:p){
    for(j in 1:p){
      A1[i,j]<-rho^(abs(i-j))
    }
  }
  A1
}

data_gen <- function(N,prob_y,nonzero_X,Z_logOR,rho,gamma_0) {
  # generate full data
  Cov <- func_A1gen(rho, (dim_X+dim_Z))
  XZ <- mvrnorm(N, rep(0,(dim_X+dim_Z)), Cov)
  
  # organize vector of coefficients for model used to generate the outcome, Y
  X_logOR <- NULL
  for (i in 1:nonzero_X) {
    X_logOR <- c(X_logOR,nonzeroX_logOR[i],rep(log(1),dim_X/nonzero_X-1))
  }
  gamma_logOR <- c(X_logOR,Z_logOR)
  gamma0=gamma_0
  exp_val = exp(gamma0 + XZ%*%as.matrix(gamma_logOR))
  epsilon_beta <- abs(c(gamma0,nonzeroX_logOR,Z_logOR))/numbeta
  
  # generate the outcome, Y
  prob <- exp_val/(1+exp_val)
  Y <- rbinom(N,1,prob)
  
  data_P1 <- data.frame(cbind(XZ,Y))
  return(data_P1)
}

Dat_gensrs <- function(prob_nII,n) {
  R <- rbinom(N,1,prob_nII)
  return(R)
}

Dat_genlinear <- function(X,tao_0) {
  Xmat <- X[,1:q]
  prob <- exp(tao_0 + 5 * rowSums(Xmat) / sqrt(q))/(1+exp(tao_0 +5 * rowSums(Xmat) / sqrt(q)))
  R <- rbinom(nrow(Xmat),1,prob)
  return(R)
}

Dat_gensingle_index <- function(X,tao_0) {
  Xmat <- X[,1:q]
  prob <- exp(tao_0 + 3 * rowSums(Xmat) / sqrt(q) + (3 / 2) * rowSums(Xmat ^ 2)/q)/(1 + exp(tao_0 + 3 * rowSums(Xmat) / sqrt(q) + (3 / 2) * rowSums(Xmat ^ 2)/q))
  R <- rbinom(nrow(Xmat),1,prob)
  return(R)
}

Dat_fullFit <- function(data,gamma_logOR) {
  tmp <- cbind(data$Y,data[,which(gamma_logOR!=0)])
  names(tmp)[1] <- "Y"
  mod <- glm(Y~.,data=tmp,family=binomial)
  yhat.full <- predict(mod,tmp,type="response")
  beta.full <- mod$coefficients
  return(list(yhat.full=yhat.full, beta.full=beta.full))
}
