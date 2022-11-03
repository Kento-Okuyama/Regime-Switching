set.seed(7)

library(bssm)
library(car)  

stateProb <- list()

##################################
# collect only drop out clusters #
##################################

for (t in ts){
  if (gmm[[t]]$parameters$mean[18,1] > gmm[[t]]$parameters$mean[18,2]) 
    doClust <- 1
  else 
    doClust <- 2
  stateProb[[t]] <- gmm[[t]]$z[,doClust]
}

###########################
# rearrange the stateProb #
###########################

x <- array(NA, c(N, Nt))
for (t in ts){
  count <- 1
  for (i in 1:N){
    if (sum(Codes[i]==gmmCodes[[t]])==1){
      x[i,t] <- stateProb[[t]][count]
      count <- count + 1
    }
  }
}
stateProb <- x

# set the initial dropout probability to be very small
stateProb[,1] <- abs(rnorm(N, mean=0, sd=0.1))

######################################
# impute NA of dropout probabilities #
######################################

for (i in 1:N){
  for (t in 1:Nt){
    if (is.na(stateProb[i,t])==TRUE)
      stateProb[i,t] <- stateProb[i,t-1]
  }
}

# intra-individual covariates
X <- yw[,,2:18]

# n covariates
nVar <- length(X[1,1,])

############################################
# impute NA of intra-individual covariates #
############################################

for (var in 1:nVar){
  for (i in 1:N){
    for (t in 2:Nt){
      if (is.na(X[i,t,var])==TRUE)
        X[i,t,var] <- X[i,t-1,var]
      }
    }
}
X[is.na(X)==TRUE] <- 0

##############################
# parameters to be estimated #
##############################
# gradient/gibbs sampler?
alpha <- array(rnorm(2, mean=0, sd=0.01), c(2,1))
beta <- array(rnorm(2, mean=0.99, sd=0.001), c(2,1))
gamma <- array(rnorm(2*nVar, mean=0, sd=0.01), c(2,nVar))
d <- rnorm(2, mean=0, sd=1)
Lmd <- array(rnorm(2, mean=0, sd=1), c(2,1))
A <- array(rnorm(2*nVar*2, mean=0, sd=0.01), c(2,nVar))

#############
# variables #
#############
K <- array(0, c(2,1))
eta <- eta_ <- array(0, c(N,Nt,2)) # elements of initial eta are fixed to be a vector of zeros
P <- P_ <- array(NA, c(N,Nt,2)); P[,1,] <- 100 # diagonal elements of initial P are set to some arbitrarily large constants 
EF <- array(NA, c(N,Nt))
y <- yw[,,19]
v <- array(0, c(N,Nt))
f <- Lf <- array(0, c(N,Nt))

Logistic <- function(x)
  exp(x) / (1 + exp(x))

#######
# EKF #
#######

for (i in 1:N){
  for (t in 2:Nt){
    eta_[i,t,] <- alpha + beta * eta[i,t-1,] + gamma %*% X[i,t,] # (A.1)
    P_[i,t,] <- beta^2 * P[i,t-1,] # (A.2)
    v[i,t] <- y[i,t] - (Logistic(d[1] + Lmd[1] * eta_[i,t,1] + A[1,] %*% X[i,t,]) * stateProb[i,t] + Logistic(d[2] + Lmd[2] * eta_[i,t,2] + A[2,] %*% X[i,t,]) * (1-stateProb[i,t])) # (A.3)
    EF[i,t] <- (Lmd[1]^2 * P_[i,t,1] + rnorm(1, mean=0, sd=1)) * stateProb[i,t] + (Lmd[2]^2 * P_[i,t,2] + rnorm(1, mean=0, sd=1)) * (1-stateProb[i,t]) # (A.4)
    # f[i,t] <- (2*pi)**(-1/2) * abs(EF[i,t])**(-1/2) * exp(-1/2 * v[i,t]**2 * EF[i,t]**(-1)) # likelihood
    Lf[i,t] <- -1/2 * (log(2*pi) + log(abs(EF[i,t])) + v[i,t]**2 * EF[i,t]**(-1)) # log-likelihood
    K <- P_[i,t,] * Lmd * EF[i,t]**(-1) # Kalman gain function
    eta[i,t,] <- eta_[i,t,] + K * v[i,t] # (A.5)
    P[i,t,] <- P_[i,t,] - K * Lmd * P_[i,t,] # (A.6)
  }
}

for (i in 1:N){
  for (t in 1:Nt){
    if (is.nan(Lf[i,t])==TRUE || is.infinite(Lf[i,t])==TRUE)
      Lf[i,t] <- Lf[i,t-1]
  }
}
sum(Lf)
