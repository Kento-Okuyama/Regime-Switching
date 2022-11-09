set.seed(7)

# library(bssm)
library(car)  

stateProb <- list()
for (t in ts){
  if (gmm[[t]]$parameters$mean[dim(gmm[[t]]$parameters$mean)[1],1] > gmm[[t]]$parameters$mean[dim(gmm[[t]]$parameters$mean)[1],2]) 
    doClust <- 1
  else 
    doClust <- 2
  stateProb[[t]] <- gmm[[t]]$z[,doClust]
}

x <- array(NA, c(N, Nt))
for (t in ts){
  count <- 1
  for (i in 1:N){
    if (sum(Codes[i] == gmmCodes[[t]])==1){
      x[i,t] <- stateProb[[t]][count]
      count <- count + 1
    }
  }
}

stateProb <- x
stateProb[,1] <- abs(rnorm(N, mean=0, sd=1e-5))
for (i in 1:N){
  for (t in 1:Nt){
    if (is.na(stateProb[i,t])==TRUE)
      stateProb[i,t] <- stateProb[i,t-1]
  }
}

# write.csv(stateProb, file='stateProb.csv')
# write.csv(y, file='y.csv')

#######################
# interpret stateProb #
#######################
# ID: 9, 58, 63, 74, 81, 82 did not drop out, and we classified them correctly 
# -> small True Negative

X <- yw[,,2:18]
nVar <- length(X[1,1,])
for (var in 1:nVar){
  for (i in 1:N){
    for (t in 2:Nt){
      if (is.na(X[i,t,var])==TRUE)
        X[i,t,var] <- X[i,t-1,var]
      }
    }
}
X[is.na(X)==TRUE] <- 0

##############
# Parameters #
##############

alpha <- array(rnorm(2, mean=0, sd=0.01), c(2,1))
beta <- array(rnorm(2, mean=0.99, sd=0.001), c(2,1))
gamma <- array(rnorm(2*nVar, mean=0, sd=0.01), c(2,nVar))
d <- rnorm(2, mean=0, sd=1)
Lmd <- array(rnorm(2, mean=0, sd=1), c(2,1))
A <- array(rnorm(2*nVar, mean=0, sd=0.01), c(2,nVar))

K <- array(0, c(2,1))
eta <- eta_ <- array(0, c(N,Nt,2))
P <- P_ <- array(NA, c(N,Nt,2)); P[,1,] <- 100
EF <- array(NA, c(N,Nt,2))
y <- yw[,,19]
v <- array(0, c(N,Nt))
f <- Lf <- array(0, c(N,Nt,2))

Logistic <- function(x)
  exp(x) / (1 + exp(x))

for (i in 1:N){
  for (t in 2:Nt){
    eta_[i,t,] <- alpha + beta * eta[i,t-1,] + gamma %*% X[i,t,]
    
    P_[i,t,] <- beta^2 * P[i,t-1,] 
    
    v[i,t] <- y[i,t] - (Logistic(d[1] + Lmd[1] * eta_[i,t,1] + A[1,] %*% X[i,t,]) * stateProb[i,t] + Logistic(d[2] + Lmd[2] * eta_[i,t,2] + A[2,] %*% X[i,t,]) * (1-stateProb[i,t]))
    
    EF[i,t,] <- Lmd^2 * P_[i,t,] + abs(rnorm(2, mean=0, sd=1))
    
    f[i,t,] <- (2*pi)**(-1/2) * EF[i,t,]**(-1/2) * exp(-1/2 * v[i,t]**2 * EF[i,t,]**(-1)) 
    
    Lf[i,t,] <- -1/2 * (log(2*pi) + EF[i,t,] + v[i,t]**2 * EF[i,t,]**(-1))  
    
    K <- P_[i,t,] * Lmd * EF[i,t,]**(-1)
    
    eta[i,t,] <- eta_[i,t,] + K * v[i,t] 
    
    P[i,t,] <- P_[i,t,] - K * Lmd * P_[i,t,]
  }
}

LL <- f[,,1] * stateProb + f[,,2] * (1-stateProb)
logLL <- Lf[,,1] * stateProb + Lf[,,2] * (1-stateProb)
for (i in 1:N){
  for (t in 1:Nt){
    if (is.nan(logLL[i,t])==TRUE)
      logLL[i,t] <- logLL[i,t-1]
  }
}
sumLogLL <- sum(logLL)

# derivatives
# step 1
dv <- array(0, c(2,1))
dF <- array(0, c(2,1))

dd <- array(0, c(2,1))
dLmd <- array(0, c(2,1))
dA <- array(0, c(2,nVar))

for (i in 1:N){
  for (t in 2:Nt){
    
    stateProbs <- array(c(stateProb[i,t], 1-stateProb[i,t]), c(2,1))
    X_it <- array(X[i,t,], c(nVar, 1))
    
    dv <- dv - v[i,t] / EF[i,t,] * Logistic(d + Lmd * eta_[i,t,] + A %*% X[i,t,]) * (1-Logistic(d + Lmd * eta_[i,t,] + A %*% X[i,t,])) * stateProbs
    
    dF <- dF - 1/2 * (EF[i,t,]**(-1) - v[i,t]**2 / EF[i,t,]**2) * stateProbs
    
    dd <- dd + dv
    # dd <- dd - v[i,t] / EF[i,t,] * Logistic(d + Lmd * eta_[i,t,] + A %*% X[i,t,]) * (1-Logistic(d + Lmd * eta_[i,t,] + A %*% X[i,t,])) * stateProbs
    
    dLmd <- dLmd + dd * eta_[i,t,] + dF * 2 * Lmd * P_[i,t,]
    # dLmd <- dLmd - v[i,t] / EF[i,t,] * Logistic(d + Lmd * eta_[i,t,] + A %*% X[i,t,]) * (1-Logistic(d + Lmd * eta_[i,t,] + A %*% X[i,t,])) * eta_[i,t,] * stateProbs -
    #   (EF[i,t,]**(-1) - v[i,t]**2 / EF[i,t,]**2) * Lmd * P_[i,t,] * stateProbs
    
    dA <- dA + dd %*% t(X_it)
    # dA <- dA - v[i,t] / EF[i,t,] * Logistic(d + Lmd * eta_[i,t,] + A %*% X[i,t,]) * (1-Logistic(d + Lmd * eta_[i,t,] + A %*% X[i,t,])) * stateProbs %*% t(X_it)
  }
}

lr <- 0.01
d <- d + lr * sign(dd)
Lmd <- Lmd + lr * sign(dLmd)
A <- A + lr * sign(dA)

dv <- array(0, c(2,1))
dF <- array(0, c(2,1))

# step 2
dalpha <- array(0, c(2,1))
dbeta <- array(0, c(2,1))
dgamma <- array(0, c(2,nVar))

for (i in 1:N){
  for (t in 2:Nt){
    
    stateProbs <- array(c(stateProb[i,t], 1-stateProb[i,t]), c(2,1))
    X_it <- array(X[i,t,], c(nVar, 1))

    dv <- dv - v[i,t] / EF[i,t,] * Logistic(d + Lmd * eta_[i,t,] + A %*% X[i,t,]) * (1-Logistic(d + Lmd * eta_[i,t,] + A %*% X[i,t,])) * stateProbs
    
    dF <- dF - 1/2 * (EF[i,t,]**(-1) - v[i,t]**2 / EF[i,t,]**2) * stateProbs
    
    dalpha <- dalpha + dd * Lmd
    # dalpha <- dalpha - v[i,t] / EF[i,t,] * Logistic(d + Lmd * eta_[i,t,] + A %*% X[i,t,]) * (1-Logistic(d + Lmd * eta_[i,t,] + A %*% X[i,t,])) * Lmd * stateProbs
    
    dbeta <- dbeta + dalpha * eta[i,t-1,] + dF * 2 * Lmd**2 * beta * P[i,t-1,] 
    # dbeta <- dbeta - [i,t] / EF[i,t,] * Logistic(d + Lmd * eta_[i,t,] + A %*% X[i,t,]) * (1-Logistic(d + Lmd * eta_[i,t,] + A %*% X[i,t,])) * eta[i,t-1,] * stateProbs -
    #   (EF[i,t,]**(-1) - v[i,t]**2 / EF[i,t,]**2) * Lmd**2 * beta * P[i,t-1,] * stateProbs
    
    dgamma <- dgamma + dalpha %*% t(X_it)
    # dgamma <- dgamma - v[i,t] / EF[i,t,] * Logistic(d + Lmd * eta_[i,t,] + A %*% X[i,t,]) * (1-Logistic(d + Lmd * eta_[i,t,] + A %*% X[i,t,])) * Lmd * stateProbs %*% t(X_it)
  }
}

lr <- 0.01
alpha <- alpha + lr * sign(dalpha)
beta <- beta + lr * sign(dbeta)
gamma <- gamma + lr * sign(dgamma)
