# implementing Newton Rhapson algorithm with ordered probit model

# simulate data
N <- 100
T <- 200

beta <- c(1,-1,-0.4)
x <- array(rnorm(N*T*3,0,1),dim = c(N,T,3))
a <- rnorm(N,0,1)
meanx_i <- rowMeans(x[,,1])
alpha <- sqrt(T)*meanx_i + a
w <- matrix(alpha,nrow=N,ncol=T,byrow=FALSE) + beta[1]*x[,,1] + beta[2]*x[,,2] + beta[3]*x[,,3]  
epsilon <- matrix(rnorm(N*T,0,1),nrow=N)
y <- (w + epsilon > 0) + (w + epsilon > 2) + (w + epsilon > 3)

# ordered probit
k <- length(beta)
beta <- rep(0,k)
tau <- c(1,2)
l <- length(tau)
gamma <- c(beta,tau)
alpha <- rep(0,N)
eps <- 10^(-3)*rep(1,length(gamma))
dif <- 100*rep(1,length(gamma))
gammaold <- NULL
alphaold <- NULL

while(sum(dif > eps)>0){
  sum1 <- matrix(0,length(gamma),length(gamma))
  sum2 <- matrix(0,length(gamma))
  
  cdfl <- list()
  pdfl <- list()
  tau2 <- c(-Inf,0,tau,Inf) # set one to 0 for identification
  w <- x[,,1]*beta[1] + x[,,2]*beta[2] + x[,,3]*beta[3] + matrix(alpha,nrow=N,ncol=T,byrow=FALSE)
  for(l in 1:length(tau2)){
    cdfl[[l]] <- pnorm(tau2[l] - w)
    pdfl[[l]] <- dnorm(tau2[l] - w)    
  }

  # derivatives to individual specific alpha and beta
  alph <- matrix(0,N,T)
  alph2 <- matrix(0,N,T)
  gbeta <- matrix(0,k,1)
  hbeta <- matrix(0,k,k)
  hbetaalpha <- matrix(0,N,k)
  ta <- list(matrix(0,N,T))
  ta2 <- list(matrix(0,N,T))
  gtau <- matrix(0,l,1)
  htau <- matrix(0,l,l)
  taalph <- list(matrix(0,N,T))
  htabeta <- list(matrix(0,k,1))
  ta12 <- matrix(0,N,T)
  for(i in 1:N){
    for(t in 1:T){
      if(y[i,t]==0){
        alph[i,t] <- (1/(cdfl[[2]][i,t]-cdfl[[1]][i,t]))*(pdfl[[1]][i,t]-pdfl[[2]][i,t])
        alph2[i,t] <- -(1/(cdfl[[2]][i,t]-cdfl[[1]][i,t])^2)*((pdfl[[1]][i,t]-pdfl[[2]][i,t])^2) + 
        (1/(cdfl[[2]][i,t]-cdfl[[1]][i,t]))*(-pdfl[[2]][i,t]*(tau2[2]-w[i,t]))
        gbeta <- gbeta + alph[i,t]*x[i,t,]
        hbeta <- hbeta + alph2[i,t]*(x[i,t,]%*%t(x[i,t,]))
        hbetaalpha[i,] <- hbetaalpha[i,] + alph2[i,t]*x[i,t,]
        #ta[[1]][i,t] <- (1/(cdfl[[2]][i,t]-cdfl[[1]][i,t]))*(pdfl[[2]][i,t]) 
        #ta2[[1]][i,t] <- -(1/(cdfl[[2]][i,t]-cdfl[[1]][i,t])^2)*(pdfl[[2]][i,t]^2) - (1/(cdfl[[2]][i,t]-cdfl[[1]][i,t]))*(pdfl[[2]][i,t]*(tau2[2]-w[i,t]))
        #taalph[[1]][i,t] <- -(1/(cdfl[[2]][i,t]-cdfl[[1]][i,t])^2)*(pdfl[[1]][i,t]-pdfl[[2]][i,t])*pdfl[[2]][i,t] +
        #  (1/(cdfl[[2]][i,t]-cdfl[[1]][i,t]))*pdfl[[2]][i,t]*(tau2[2]-w[i,t]) 
        #htabeta[[1]] <- htabeta[[1]] + taalph[[1]][i,t]*x[i,t,]
      }
      else if(y[i,t]==1){
        alph[i,t] <- (1/(cdfl[[3]][i,t]-cdfl[[2]][i,t]))*(pdfl[[2]][i,t]-pdfl[[3]][i,t])
        alph2[i,t] <- -(1/(cdfl[[3]][i,t]-cdfl[[2]][i,t])^2)*((pdfl[[2]][i,t]-pdfl[[3]][i,t])^2) + 
        (1/(cdfl[[3]][i,t]-cdfl[[2]][i,t]))*(pdfl[[2]][i,t]*(tau2[2]-w[i,t])-pdfl[[3]][i,t]*(tau2[3]-w[i,t]))
        gbeta <- gbeta + alph[i,t]*x[i,t,]
        hbeta <- hbeta + alph2[i,t]*(x[i,t,]%*%t(x[i,t,]))
        hbetaalpha[i,] <- hbetaalpha[i,] + alph2[i,t]*x[i,t,]
        #ta[[1]][i,t] <- - (1/(cdfl[[3]][i,t]-cdfl[[2]][i,t]))*(pdfl[[2]][i,t])
        #ta2[[1]][i,t] <- -(1/(cdfl[[3]][i,t]-cdfl[[2]][i,t])^2)*(pdfl[[2]][i,t]^2) + (1/(cdfl[[3]][i,t]-cdfl[[2]][i,t]))*(pdfl[[2]][i,t]*(tau2[2]-w[i,t]))
        #taalph[[1]][i,t] <- (1/(cdfl[[3]][i,t]-cdfl[[2]][i,t])^2)*(pdfl[[2]][i,t]-pdfl[[3]][i,t])*pdfl[[2]][i,t] -
        #  (1/(cdfl[[3]][i,t]-cdfl[[2]][i,t]))*pdfl[[2]][i,t]*(tau2[2]-w[i,t])
        #htabeta[[1]] <- htabeta[[1]] + taalph[[1]][i,t]*x[i,t,]
        ta[[1]][i,t] <- (1/(cdfl[[3]][i,t]-cdfl[[2]][i,t]))*(pdfl[[3]][i,t]) 
        ta2[[1]][i,t] <- -(1/(cdfl[[3]][i,t]-cdfl[[2]][i,t])^2)*(pdfl[[3]][i,t]^2) - (1/(cdfl[[3]][i,t]-cdfl[[2]][i,t]))*(pdfl[[3]][i,t]*(tau2[3]-w[i,t]))
        taalph[[1]][i,t] <- -(1/(cdfl[[3]][i,t]-cdfl[[2]][i,t])^2)*(pdfl[[2]][i,t]-pdfl[[3]][i,t])*pdfl[[3]][i,t] +
          (1/(cdfl[[3]][i,t]-cdfl[[2]][i,t]))*pdfl[[3]][i,t]*(tau2[3]-w[i,t]) 
        htabeta[[1]] <- htabeta[[1]] + taalph[[1]][i,t]*x[i,t,]
        #ta12[i,t] <- (1/(cdfl[[3]][i,t]-cdfl[[2]][i,t])^2)*pdfl[[2]][i,t]*pdfl[[3]][i,t]
      }
      else{
        alph[i,t] <- (1/(cdfl[[4]][i,t]-cdfl[[3]][i,t]))*(pdfl[[3]][i,t]-pdfl[[4]][i,t])
        alph2[i,t] <- -(1/(cdfl[[4]][i,t]-cdfl[[3]][i,t])^2)*((pdfl[[3]][i,t]-pdfl[[4]][i,t])^2) + 
        (1/(cdfl[[4]][i,t]-cdfl[[3]][i,t]))*(pdfl[[3]][i,t]*(tau2[3]-w[i,t]))
        gbeta <- gbeta + alph[i,t]*x[i,t,]
        hbeta <- hbeta + alph2[i,t]*(x[i,t,]%*%t(x[i,t,]))
        hbetaalpha[i,] <- hbetaalpha[i,] + alph2[i,t]*x[i,t,]
        ta[[1]][i,t] <- - (1/(cdfl[[4]][i,t]-cdfl[[3]][i,t]))*(pdfl[[3]][i,t])
        ta2[[1]][i,t] <- -(1/(cdfl[[4]][i,t]-cdfl[[3]][i,t])^2)*(pdfl[[3]][i,t]^2) + (1/(cdfl[[4]][i,t]-cdfl[[3]][i,t]))*(pdfl[[3]][i,t]*(tau2[3]-w[i,t]))
        taalph[[1]][i,t] <- (1/(cdfl[[4]][i,t]-cdfl[[3]][i,t])^2)*(pdfl[[3]][i,t]-pdfl[[4]][i,t])*pdfl[[3]][i,t] -
          (1/(cdfl[[4]][i,t]-cdfl[[3]][i,t]))*pdfl[[3]][i,t]*(tau2[3]-w[i,t])
        htabeta[[1]] <- htabeta[[1]] + taalph[[1]][i,t]*x[i,t,]
      }
    }
  }
  galpha <- rowSums(alph)
  hii <- rowSums(alph2)

  # derivatives to tau
  #gtau <- c(sum(ta[[1]]),sum(ta[[2]]))
  gtau <- sum(ta[[1]])
  #htau <- matrix(c(sum(ta2[[1]]),sum(ta12),sum(ta12),sum(ta2[[2]])),nrow = 2,byrow=TRUE)
  htau <- sum(ta2[[1]])
    
  # cross derivatives
  #htaualpha <- cbind(rowSums(taalph[[1]]),rowSums(taalph[[2]]))
  htaualpha <- rowSums(taalph[[1]])
  hgammai <- cbind(hbetaalpha,htaualpha)
  #htaubeta <- cbind(htabeta[[1]],htabeta[[2]])
  htaubeta <- htabeta[[1]]
    
  ggamma <- c(gbeta,gtau)
  Hgamma <- matrix(0,length(gamma),length(gamma))
  Hgamma[1:k,1:k] <- hbeta
  Hgamma[(k+1),1:k] <- htaubeta
  Hgamma[1:k,(k+1)] <- t(htaubeta)
  Hgamma[(k+1),(k+1)] <- htau
  
  # Newton Rhapson step
  for(i in 1:N){
    sum1 <- sum1 + (1/hii[i])*hgammai[i,]%*%t(hgammai[i,])
    sum2 <- sum2 + (galpha[i]/hii[i])*hgammai[i,]
  }
  deltagamma <- -solve(Hgamma - sum1)%*%(ggamma - sum2)
  deltaalpha <- -(1/hii)*(galpha + hgammai%*%deltagamma)
  
  gammaold <- cbind(gammaold,gamma)
  alphaold <- cbind(alphaold,alpha)
  gamma <- gamma + deltagamma
  alpha <- alpha + deltaalpha
  beta <- gamma[1:k]
  tau <- gamma[k+1]
  dif <- abs(gammaold[,ncol(gammaold)]-gamma)
}
