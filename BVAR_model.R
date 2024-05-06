# ---------------------------------------------
# clean, set directory
# rm(list = ls())
# graphics.off()

currentDir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(currentDir)
# -------------------------------------------------------
# libraries
# .libPaths(new = "N://R_LIBRARY4")
libs <- c("stats", "pracma",'mfilter',"tidyverse",'zoo',"KFAS","expm","readxl",'lubridate')
lapply(libs, require, character.only = T)


##------------------------------------------------
# state space BVAR model
##------------------------------------------------
# can handle missing observation a conditional forecast
##------------------------------------------------
## FUN
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






# Y <-  cbind(Ygap,Pie.qoq,W.qoq,Ir)
Y <-  detr_dataset

Ndraws <-  10000
burnin <- Ndraws/2 # Number of burn-in draws
store <- Ndraws - burnin
p <- 2  # NUmber of lags in the VAR for the cycle;
pred_ahead <- 20






T <- nrow(Y)
n <- ncol(Y)







  # MINESSOTA PRIOR, ALE CENTRUJEME TO NA NULU I NA PRVNIM LAGU
b0 <- matrix(0, n * p, n)
diag(b0[1:n, 1:n]) <- rep(1, n)*0

# exogeni promene
# e_n=1
# b0 <- matrix(0, e_n+n * p, n)
# diag(b0[1:n, 1:n]) <- rep(1, n)*0




# slavek bussiness cycle priors
diag(b0[1:n, 1:n]) <- rep(1, n)*1.4
diag(b0[(n+1):(n*2), 1:n]) <- rep(1, n)*-0.5



# varience minesota prioru nechavam pocitat ve var funkci
# 
# PSI=c(1,1,1,1,1)
#   
# PSI=rep(1,dim(Y)[2])

  
# 
# COVARIENCE full sample
# v zasade je to jedno, protoze ve VAR modelu si ted pocitam podle minnesoty
Psi <- diag(cov(na.omit(Y)))
  



# observation equation
C  <-  matrix(0, n,n*p)
C[1:n,1:n] <- diag(rep(1, n))
  
  

# transition equation
A <- matrix(0, n * p, n * p)
A[(n + 1):(n * p), 1:(n * p - n)] <- diag(rep(1, n * (p - 1)))
  

# error terms in observation equation,  v podstate nula
R <- diag(rep(1e-12, n))
  

# covariance matrix
Q <- matrix(0, nrow = n*p, ncol = n*p)
for (ii in 1:n) {
    Q[ii,ii] <- Psi[ii]
  }
  

  
#  expected values of the initial states
S0 <- rep(0, n * p)

# covariance matrix of the nondiffuse part of the initial state vector
P0 <- diag(c(rep(Psi, p)))

  
  
  
  
# --------------------------------
# data holders
P_acc <- rep(NA, store)
States <- array(NA, dim = c(T,  n*p, store))
LogLik <- rep(NA, store)
Theta <- rep(NA, store)
AA <- array(NA, dim = c( n*p,  n*p, store))
QQ <- array(NA, dim = c( n*p,  n*p, store))
CC <- array(NA, dim = c(n,  n*p, store))
RR <- array(NA, dim = c(n, n, store))
pred_array <- array(NA, dim = c(pred_ahead, n, store))
# -----------

  
# BVAR model with missing observation
  
  
for (jm in 1:Ndraws) {
# # sampluji states
# simulation smoother
sim <- simulateSSM(SSModel(as.matrix(Y) ~ -1 + SSMcustom(Z = C, T = A,  Q = Q, a1 = S0, P1 = P0), H = R), type="states")  # simulation smoother
    
# Ycyc <- sim[, 1:n,1]
Ycyc <- Y

# nahrazuji jenom missing values
if(any(is.na(Ycyc))){
  Ycyc_pomocna=as.numeric((as.matrix(Ycyc)))
  Ycyc_pomocna[which(is.na(as.numeric((as.matrix(Ycyc)))))]=as.numeric(sim[, 1:n,1])[which(is.na(Ycyc))]
  Ycyc=matrix(Ycyc_pomocna,nrow=dim(Ycyc)[1],ncol=dim(Ycyc)[2])
  
}



# mean zero.... at nemusim resit kontantu
colmean_holder <- colMeans(Ycyc)
Ycyc <- sweep(Ycyc, 2, colMeans(Ycyc))
# ----------------
    

# BVAR
result <- BVAR(Ycyc, p, b0,   1)
    
beta <- result$beta
sigma <- result$sigma
A[1:n, 1:( n * p)] <- t(beta)
Q[1:n,  1: n] <- sigma


    
    
# covariance matrix of the nondiffuse part of the initial state vector
# klasicka dlouhodoba unconditional variance
vecP0full <- pinv(diag((n*p)^2) - kronecker(A, A)) %*% as.vector(Q)
P0full <- matrix(vecP0full, nrow = ( n * p), ncol = ( n * p))
P0[1:nrow(P0), (1):ncol(P0)] <- P0full[1:nrow(P0), ( 1):ncol(P0)]
    
    
# predikce
pred <- predict(SSModel(as.matrix(Y) ~ -1 + SSMcustom(Z = C, T = A,  Q = Q, a1 = S0, P1 = P0), H = R), n.ahead=pred_ahead,interval="none",level=0.9)
pred <-  matrix(unlist(pred), ncol=n)
    


if (jm > burnin) {
      # draws_a[, jm - burnin] <- a
      # draws_sigma[, jm - burnin] <- u_sigma
      
      
      States[,,jm - burnin] <-  as.matrix(sim[,,1])
      # LogLik[jm] <- loglik
      AA[,,jm - burnin] <- A
      QQ[,,jm - burnin] <- Q
      CC[,,jm - burnin] <- C
      RR[,,jm - burnin] <- R
      pred_array[,,jm - burnin] <- pred
      # 
      
      
      
    }
    
    if(jm %% 1000==0) {
      # Print on the screen some message
      cat(paste0("iteration: ", jm, "\n"))
    }
    
    
}
  

predikce_median  <-  bind_cols(
  apply(pred_array[,1,], 1, median),
  apply(pred_array[,2,], 1, median),
  apply(pred_array[,3,], 1, median),
  apply(pred_array[,4,], 1, median),
  apply(pred_array[,5,], 1, median),
  apply(pred_array[,6,], 1, median),
  apply(pred_array[,7,], 1, median),
)
colnames(predikce_median)   <-  colnames(Y)


prazdna=matrix(NA,ncol=dim(Ycyc)[2],nrow=dim(Ycyc)[1])
colnames(prazdna)   <-  colnames(Y)

prazdna[dim(Ycyc)[1],]=as.matrix(Y[dim(Ycyc)[1],])

predikce_graf=rbind(prazdna,predikce_median)

ts.plot(predikce_graf[,2],ylim=c(-5,10))
lines(Y[,2],col='red')


# median beta/sigma
A_median <-  apply(AA,c(1,2), median)
sigma_median <-  apply(QQ,c(1,2), median)[1:n,1:n]




irf_vysledek <-  irf_function_BIG(A=A,sigma=sigma,nstep=11,n=n,p=p)
irf_vysledek <-  irf_function_BIG(A=A_median,sigma=sigma_median,nstep=11,n=n,p=p)
# # # # IRF vzhlem k prvni promene
ts.plot(t(irf_vysledek[,1,]),col='red')

irf_list <-  list()
irf_list2 <-  list()

for (ii in 1:store) {
  res=irf_function_BIG(A=AA[,,ii],sigma=QQ[1:n,1:n,ii],nstep=11,n=n,p=p)
  irf_list[[ii]]=res
  irf_list2[[ii]]=t(res[1,1,])
  
}
ts.plot(matrix(unlist(irf_list2),ncol=5000,nrow=11))


OLS_VAR <-  VAR_model((Ycyc),2)
irf_vysledek=irf_function(OLS_VAR$beta,OLS_VAR$sigma,nstep=11,n=n,p=p)
lines((irf_vysledek[1,1,]),col='red')
# # # # IRF vzhlem k prvni promene
ts.plot(t(irf_vysledek[,1,]),col='red')

a <-  fevd_fnc(rbind(t(OLS_VAR$beta),cbind(diag(4),diag(4)*0)),OLS_VAR$sigma,n=n,p=p)
a$VD[,,1]
a=fevd_fnc(A_median,sigma_median,n=n,p=p)
ts.plot(a$VD[,,1])
BVAR_res <-  BVAR(Y,2,b0)

  
  

irf_vysledek=irf_function(BVAR_res$beta,BVAR_res$sigma,nstep=11,n=n,p=p)
ts.plot(t(irf_vysledek[,1,]),col='red')




HDD=VAR_historical_decomp(residuals_res=OLS_VAR$eps,A=rbind(t(OLS_VAR$beta),cbind(diag(4),diag(4)*0)), p=p,sigma=OLS_VAR$sigma,n=n)

ts.plot(HDD[,,1])

