##------------------------------------------------
## FUNS
##------------------------------------------------
#  IMPULSE RESPONSE 
# be aware of matrix orientation (reg koeficienty pro jednu promenou jsou sloupce!!)
# recursive computation ALA CHRIS SIMS
irf_function <- function(beta=beta,sigma=sigma,nstep=11,n=n,p=p){
  neq <- n #  pomocna promena
  response <- array(0,dim = c(neq,n,nstep))
  response[,,1] <- t(chol(sigma)) # choleski ordering
  # response[,,1] <- diag(n) # nerekurzivne 
  B <- array(0,dim <-  c(neq,n,p))
  for (pp in 1:p) {
    B[,,pp]<-t(beta)[1:n,(n*pp-n+1):(n*pp)]
  }
  
  # rekurzivni vypocet IRF (prepsat do C++) 
  for (it in 2:nstep) {
    for (ilag in 1:min(p,(it-1))) {
      response[,,it]<-response[,,it]+B[,,ilag]%*% response[,,it-ilag]
    }
  }
  return(response)
}
##------------------------------------------------

ewma <- function(x, v = 0.98, y0) {
  ## exponentially-weighted moving average as in Cieslak-Povala
  ## v = 0.95 quarterly data
  tau <- numeric(length(x))
  if (missing(y0)) {
    tau[1] <- x[1]
  } else {
    tau[1] <- y0
  }
  for (t in 2:length(x))
    tau[t] <- tau[t-1] + (1-v)*(x[t] - tau[t-1])
  tau
}


##------------------------------------------------
#  IMPULSE RESPONSE VAR(1) approach
# vyzaduje expm package!!
# jednodusi na pochopeni
# pozor na orientaci (reg koeficinty pro jednu promenou radky)
# irf_function_BIG a irf_function musi davat stejny vysledek!!!

irf_function_BIG <- function(A=A,sigma=sigma,nstep=11,n=n,p=p){
  neq <- n
  response <- array(0,dim = c(neq,n,nstep))
  JJ<- cbind(diag(n), matrix(0,n,n*(p-1)))
  for (it in 0:(nstep-1)) {
    response[,,(it+1)]<-JJ %*% (A %^% it) %*% t(JJ)%*% t(chol(sigma)) # choleski ordering
    # response[,,(it+1)]<-JJ %*% (A %^% it) %*% t(JJ)%*% diag(n) # nenresursive ordering 
  }
  # # # # IRF vzhledem k prvni promene
  # ts.plot(t(response[,1,]),col='red')
  return(response)
}
##------------------------------------------------

##------------------------------------------------
# simple OLS VAR(p)
VAR_model<- function(Ycyc, p, exogen)  {
  Ycyc <- as.matrix(Ycyc)
  TT <- nrow(Ycyc)
  n <- ncol(Ycyc)
  k <- n * p
  # constructing the matrix of regressors
  x <- matrix(0, nrow = TT, ncol = k)
  for (i in 1:p) {
    x[, ((i-1)*n+1):(i*n)] <- as.matrix(dplyr::lag(Ycyc,n=i))
  }
  x <- x[(p+1):TT, ]
  y <- Ycyc[(p+1):TT, ]
  T <- nrow(y)
  if(missing(exogen)){
    beta <- pracma::mldivide((crossprod(x,x) ), (t(x) %*% y ), pinv = FALSE)
  }else{
    exogen <- as.matrix(exogen)
    exogen <-  exogen[(p+1):TT,] 
    x <- cbind(x,exogen)
    beta <- pracma::mldivide((crossprod(x,x) ), (t(x) %*% y ), pinv = FALSE)
    
  }
  
  eps <- y - x %*% beta
  sigma <- ((crossprod(eps,eps))) / (T-p-p*n)
  VAR_LIST <- list(beta=beta, sigma=sigma,eps=eps)
  return(VAR_LIST)
  
}
##------------------------------------------------


##------------------------------------------------
# Forecast error variance decomposition
# 
fevd_fnc <- function(A, p=p,sigma,nstep=11,n=n)  {
  neq <- n
  # IRF vypocet
  PSI <- array(0,dim = c(neq,n,nstep))
  JJ <- cbind(diag(n), matrix(0,n,n*(p-1)))
  for (it in 0:(nstep-1)) {
    PSI[,,(it+1)] <- JJ %*% (A %^% it) %*% t(JJ)
  }
  
  # 
  VD <- array(0,dim = c(nstep,n,n)) # varience demposition
  SE <- matrix(0,ncol=n,nrow=nstep) # standard error
  
  
  # loop for the shocks
  for (ii in 1:n) {
    # jmenovatel
    MSE <- array(0,dim = c(n,n,nstep))
    MSE[,,1] <- sigma
    # MSE[,,1] <- response[,,1] %*% (response[,,1]) # alternativne z vysledku IRF
    for (nn in 2:nstep) {
      MSE[,,nn] <- MSE[,,nn-1]+ PSI[,,nn]%*% sigma%*% t(PSI[,,nn])
      # MSE[,,nn] <- MSE[,,nn-1]+ response[,,nn]%*% t(response[,,nn]) # alternativne z vysledku IRF
      
    }
    
    # citatel
    MSE_shock <- array(0,dim = c(n,n,nstep))
    B <- t(chol(sigma))
    MSE_shock[,,1]= (B[,ii])%*%t(B[,ii])
    # MSE_shock[,,1]= crossprod(t(response[,ii,1]),(response[,ii,1])) # alternativne z vysledku IRF
    
    for (nn in 2:nstep) {
      MSE_shock[,,nn] <- MSE_shock[,,nn-1]+ PSI[,,nn]%*% MSE_shock[,,1]%*% t(PSI[,,nn])
      # MSE_shock[,,nn] <- MSE_shock[,,nn-1]+  response[,ii,nn] %*%t(response[,ii,nn])  # alternativne z vysledku IRF

    }
    # Compute the Forecast Error Covariance Decomposition
    FECD <- MSE_shock[1:n,1:n,]/MSE[1:n,1:n,]
    
    # Select only the variance (diagonal) terms
    
    for (nn in 1:nstep) {
      for (kk in 1:n) {
        VD[nn,ii,kk]  <-  100*FECD[kk,kk,nn]
        SE[nn,]  <-  sqrt(t(diag(MSE[1:n,1:n,nn])))
        
      }
      
    }
    
  }
  
  FEVD_LIST <- list(VD=VD, SE=SE)
  return(FEVD_LIST)
  
}
##------------------------------------------------
##------------------------------------------------

# historical decomposition
VAR_historical_decomp <- function(residuals_res,A, p=p,sigma,n=n){
  ## Retrieve and initialize variables 
  invA    <- t(chol(sigma))           # inverse of the A matrix
  eps     <- pinv(invA) %*% t(residuals_res)          # structural errors 
  nobs    <- nrow(residuals_res)                  # number of observations
  
  ## Compute historical decompositions
  
  # Contribution of each shock
  invA_big <- matrix(0,n*p,n)
  invA_big[1:n,] <- invA
  
  Icomp <- cbind(diag(n), matrix(0,n,(p-1)*n))
  HDshock_big <- array(0, dim=c(p*n,nobs+1,n))
  HDshock <- array(0, dim=c(n,(nobs+1),n))
  
  for (j in 1:n){  # for each variable
    eps_big <- matrix(0,n,(nobs+1)) # matrix of shocks conformable with companion
    eps_big[j,2:ncol(eps_big)] <- eps[j,]
    for (i in 2:(nobs+1)){
      HDshock_big[,i,j] <- invA_big %*% eps_big[,i] + A %*% HDshock_big[,(i-1),j]
      HDshock[,i,j] <-  Icomp %*% HDshock_big[,i,j]
    } 
    
  } 
  
  #### Initial value contribution
  # X treba dodat X (matice regressoru)
  # HDinit_big = matrix(0,n*p,nobs+1)
  # HDinit =  matrix(0,n,nobs+1)
  # HDinit_big[,1] = t([1,])
  # HDinit_big[,1] = Icomp %*% HDinit_big[,1]
  # for (i in 2:(nobs+1)) {
  #   HDinit_big[,i,j] <- invA_big %*% eps_big[,i] + A %*% HDinit_big[,(i-1)]
  #   HDinit[,i,j] <-  Icomp %*% HDinit_big[,i]
  #   
  # }
  
  
  HD <- array(0, dim=c((nobs+p),n,n))   # [nobs x shock x var]
  
  for (i in 1:n){
    
    for (j in 1:n){
      HD[,j,i] <- c(rep(NA,p), HDshock[i,(2:dim(HDshock)[2]),j])
    }
  }
  
  return(HD)
  
}
##------------------------------------------------




##------------------------------------------------
# simple BVAR model
# minnesota style prior (prior)
# conjugate normal-inverse-Wishart model
BVAR <- function(Ycyc, p, b,   draw=1,PSI, exogen,alpha = 2)  {
  # data matrix manipulations
  TT <- nrow(Ycyc)
  n <- ncol(Ycyc)
  k <- n * p
  e_dim <- 0
  if(!missing(exogen)){
    e_dim <- dim(as.matrix(exogen))[2]
    
  }
  
  # constructing the matrix of regressors
  x <- matrix(0, nrow = TT, ncol = k)
  for (i in 1:p) {
    x[, ((i-1)*n+1):(i*n)] <- as.matrix(dplyr::lag(Ycyc,n=i))
  }
  # 
  # 
  
  # if missing observation
  if(missing(exogen)){
    x <- x[(p+1):TT, ]
    y <- Ycyc[(p+1):TT, ]
    y0 <- ifelse(p==1,mean(y[(1:p), ]),rowMeans(y[(1:p), ]))
    T <- nrow(y)
  }else{
    x <- x[(p+1):TT, ]
    exogen <- as.matrix(exogen)
    exogen <-  exogen[(p+1):(TT),]
    x <- cbind(x,exogen)
    
    y <- Ycyc[(p+1):TT, ]
    y0 <- ifelse(p==1,mean(y[(1:p), ]),rowMeans(y[(1:p), ]))
    T <- nrow(y)
    #   
  }
  
  
  # MINNESOTA PRIOR
  #  initiate matrix
  s_endo <- apply(y, 2, stats::var)
  
  #  AR 1 process
  for (index in 1:dim(y)[2]) {
    y_pomocna <- y[,index]
    s_endo[index] <-  (summary(lm(y_pomocna ~ -1+lag(y_pomocna)))$sigma)^2
    
  }
  
  
  s_endo <- (diag(s_endo))
  
  if(missing(PSI)){
    PSI <- diag(s_endo)
  }else{
    PSI <- PSI
  }
  
  
  # starting values for the minimization
  lambda  <-  .2 # std of MN prior
  # alpha <- 2 # lag-decaying parameter of the MN prior
  d <- n + 2 # df for the covariance
  
  # priors
  
  omega <- rep(0, k)
  for (i in 1:p) {
    omega[((i-1)*n+1):(i*n)] <- (d - n - 1) * (lambda^2) * (1/(i^alpha)) / PSI
  }
  
  
  
  
  Omega_inv <- diag(1 / omega)
  
  Omega_inv_big<- diag(k+e_dim)*0
  Omega_inv_big[1:k,1:k] <- Omega_inv

  if(!missing(exogen)){
    Omega_inv_big[((k+1):(k+e_dim)),((k+1):(k+e_dim))] <- 100    
  }
  
  Omega_inv_big[1:k,1:k] <- Omega_inv
  
  Omega_inv <- Omega_inv_big
  
  
  
  ### make uniformative----------------
  # PSI <- diag(s_endo)*0
  # Omega_inv <- 0*Omega_inv
  # d <- 0
  # -------------------------------
  
  
  
  
  # dummy observation
  # SOC prior
  # miu = 3 # std of SOC
  # ydnoc = (1/miu) * diag(y0) * 0
  # xdnoc = matrix(rep(0, n * p), nrow = n, ncol = p)
  # xdnoc = t(xdnoc)
  # y <- rbind(y, ydnoc)
  # x <- rbind(x, xdnoc)
  # T <- T + n
  # 
  
  
  
  # posterior mode of the VAR coefficients
  betahat <- pracma::mldivide((crossprod(x,x) + Omega_inv), (t(x) %*% y + Omega_inv %*% b), pinv = FALSE)
  
  
  # VAR residuals
  epshat <- y - x %*% betahat
  sigmahat <- ((crossprod(epshat,epshat)) + diag(PSI) + t(betahat - b) %*% Omega_inv %*% (betahat - b)) / (T + d + n + 1)
  
  
  
  stationary <- 0
  if(draw == 1){
    
    # nutim to byt stacionarni
    while (stationary == 0) {
      Sinv <- pracma::pinv((sigmahat * (T + d + n + 1)))
      eta <- MASS::mvrnorm((T+d),mu=rep(0,n),Sigma=Sinv)
      sigma <- pracma::mldivide(crossprod(eta,eta), diag(n)) #  pinv(crossprod(eta,eta)) %*% diag(n)
      cholSIGMA <- RobustCalibration::Chol_Eigen(sigma) 
      
      # beta <- betahat + t(MASS::mvrnorm(n, mu=rep(0,k),Sigma=mldivide(crossprod(x,x)+Omega_inv,diag(k), pinv = FALSE)))%*%cholSIGMA
      beta <- betahat + t(MASS::mvrnorm(n, mu=rep(0,k+e_dim),Sigma=mldivide(crossprod(x,x)+Omega_inv,diag(k+e_dim))))%*%cholSIGMA
      
      
      # kontroluju jenom stacionaritu endogenich promenych
      AA <- matrix(0, n * p, n * p)
      AA[1:n, 1:(n * p)] <- t(beta[(1:(n*p)),])
      if (p>1){
        AA[(n + 1):(n * p), 1:(n * (p - 1))] <- diag(rep(1, n * (p - 1)))
      }
      stationary <- sum(abs(eigen(AA)$values) > 1) == 0
      
    }
    
  } else {
    beta <- betahat
    sigma <- sigmahat 
  }
  
  
  
  # final residuals
  eps <- y - x %*% beta
  
  
  BVAR_LIST <- list(beta=beta, sigma=sigma,eps=eps,epshat=epshat)
  return(BVAR_LIST)
  
}
##------------------------------------------------

