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
source('BVAR_model_fun.R')


# create/select dataset
Y <-  detr_dataset

Ndraws <-  10000
burnin <- Ndraws/2 # Number of burn-in draws
store <- Ndraws - burnin
p <- 2  # NUmber of lags
pred_ahead <- 20






T <- nrow(Y)
n <- ncol(Y)




  # MINESSOTA PRIOR, centered on 0
b0 <- matrix(0, n * p, n)
diag(b0[1:n, 1:n]) <- rep(1, n)*0

# exogeni promene
# e_n=1
# b0 <- matrix(0, e_n+n * p, n)
# diag(b0[1:n, 1:n]) <- rep(1, n)*0




# bussiness cycle priors
#diag(b0[1:n, 1:n]) <- rep(1, n)*1.4
#diag(b0[(n+1):(n*2), 1:n]) <- rep(1, n)*-0.5




  
# 
# COVARIENCE full sample
# just for first draw
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
# # sampling states
# simulation smoother
sim <- simulateSSM(SSModel(as.matrix(Y) ~ -1 + SSMcustom(Z = C, T = A,  Q = Q, a1 = S0, P1 = P0), H = R), type="states")  # simulation smoother
    
# Ycyc <- sim[, 1:n,1]
Ycyc <- Y

# replace missing values
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
# unconditional variance
vecP0full <- pinv(diag((n*p)^2) - kronecker(A, A)) %*% as.vector(Q)
P0full <- matrix(vecP0full, nrow = ( n * p), ncol = ( n * p))
P0[1:nrow(P0), (1):ncol(P0)] <- P0full[1:nrow(P0), ( 1):ncol(P0)]
    
    
# prediction
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
      # iteration
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
# # # # IRF relative to first variable
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

