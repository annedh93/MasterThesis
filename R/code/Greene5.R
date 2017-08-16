# implementing Newton Rhapson algorithm with ordered probit model

ptm <- proc.time()
# simulate data
N <- 100
T <- 200

beta <- c(1,-1) # change this for different number of parameters
k <- length(beta)
x <- array(rnorm(N*T*k,0,1),dim = c(N,T,k))
a <- rnorm(N,0,1)
meanx_i <- rowMeans(x[,,1])
alpha <- sqrt(T)*meanx_i + a
w <- matrix(alpha,nrow=N,ncol=T,byrow=FALSE) + eval(parse(text=paste(paste0("beta","[",1:k,"]*","x","[,,",1:k,"]",collapse="+"))))
epsilon <- matrix(rnorm(N*T,0,1),nrow=N)
y <- (w + epsilon > 0) + (w + epsilon > 2) + (w + epsilon > 3) # change this for different (number of) limits

# ordered probit
beta <- rep(0,k)
tau <- c(0.1,2) # make sure initialisation is separated enough # change this for different (number of) limits
p <- length(tau)
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
  w <- eval(parse(text=paste(paste0("beta","[",1:k,"]*","x","[,,",1:k,"]",collapse="+")))) + matrix(alpha,nrow=N,ncol=T,byrow=FALSE)
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
  ta <- array(0,dim=c(N,T,p))
  ta2 <- array(0,dim=c(N,T,p))
  gtau <- matrix(0,p,1)
  htau <- matrix(0,p,p)
  taalph <- array(0,dim=c(N,T,p))
  htabeta <- matrix(0,k,p)
  if(p > 1){
    ta12 <- array(0,dim=c(N,T,(p-1)))    
  }
  for(i in 1:N){
    for(t in 1:T){
      c <- y[i,t] + 2
      q <- y[i,t]
      if(y[i,t]==0){
        alph[i,t] <- (1/(cdfl[[c]][i,t]-cdfl[[c-1]][i,t]))*(pdfl[[c-1]][i,t]-pdfl[[c]][i,t])
        alph2[i,t] <- -(1/(cdfl[[2]][i,t]-cdfl[[1]][i,t])^2)*((pdfl[[1]][i,t]-pdfl[[2]][i,t])^2) + 
          (1/(cdfl[[c]][i,t]-cdfl[[c-1]][i,t]))*(-pdfl[[c]][i,t]*(tau2[c]-w[i,t]))
        gbeta <- gbeta + alph[i,t]*x[i,t,]
        hbeta <- hbeta + alph2[i,t]*(x[i,t,]%*%t(x[i,t,]))
        hbetaalpha[i,] <- hbetaalpha[i,] + alph2[i,t]*x[i,t,]
      }
      else if(y[i,t]==1){
        alph[i,t] <- (1/(cdfl[[c]][i,t]-cdfl[[c-1]][i,t]))*(pdfl[[c-1]][i,t]-pdfl[[c]][i,t])
        alph2[i,t] <- -(1/(cdfl[[c]][i,t]-cdfl[[c-1]][i,t])^2)*((pdfl[[c-1]][i,t]-pdfl[[c]][i,t])^2) + 
          (1/(cdfl[[c]][i,t]-cdfl[[c-1]][i,t]))*(pdfl[[c-1]][i,t]*(tau2[c-1]-w[i,t]) - pdfl[[c]][i,t]*(tau2[c]-w[i,t]))
        gbeta <- gbeta + alph[i,t]*x[i,t,]
        hbeta <- hbeta + alph2[i,t]*(x[i,t,]%*%t(x[i,t,]))
        hbetaalpha[i,] <- hbetaalpha[i,] + alph2[i,t]*x[i,t,]
        # first limit
        ta[i,t,q] <- (1/(cdfl[[c]][i,t]-cdfl[[c-1]][i,t]))*(pdfl[[c]][i,t]) 
        ta2[i,t,q] <- -(1/(cdfl[[c]][i,t]-cdfl[[c-1]][i,t])^2)*(pdfl[[c]][i,t]^2) - (1/(cdfl[[c]][i,t]-cdfl[[c-1]][i,t]))*(pdfl[[c]][i,t]*(tau2[c]-w[i,t]))
        taalph[i,t,q] <- -(1/(cdfl[[c]][i,t]-cdfl[[c-1]][i,t])^2)*(pdfl[[c-1]][i,t]-pdfl[[c]][i,t])*pdfl[[c]][i,t] +
          (1/(cdfl[[c]][i,t]-cdfl[[c-1]][i,t]))*pdfl[[c]][i,t]*(tau2[c]-w[i,t]) 
        htabeta[,q] <- htabeta[,q] + taalph[i,t,q]*x[i,t,]
      }
      else if( (y[i,t]>1) & (y[i,t]<(p+1)) ){
        alph[i,t] <- (1/(cdfl[[c]][i,t]-cdfl[[c-1]][i,t]))*(pdfl[[c-1]][i,t]-pdfl[[c]][i,t])
        alph2[i,t] <- -(1/(cdfl[[c]][i,t]-cdfl[[c-1]][i,t])^2)*((pdfl[[c-1]][i,t]-pdfl[[c]][i,t])^2) + 
          (1/(cdfl[[c]][i,t]-cdfl[[c-1]][i,t]))*(pdfl[[c-1]][i,t]*(tau2[c-1]-w[i,t]) - pdfl[[c]][i,t]*(tau2[c]-w[i,t]))
        gbeta <- gbeta + alph[i,t]*x[i,t,]
        hbeta <- hbeta + alph2[i,t]*(x[i,t,]%*%t(x[i,t,]))
        hbetaalpha[i,] <- hbetaalpha[i,] + alph2[i,t]*x[i,t,]
        # first limit
        ta[i,t,(q-1)] <- - (1/(cdfl[[c]][i,t]-cdfl[[c-1]][i,t]))*(pdfl[[c-1]][i,t])
        ta2[i,t,(q-1)] <- -(1/(cdfl[[c]][i,t]-cdfl[[c-1]][i,t])^2)*(pdfl[[c-1]][i,t]^2) + (1/(cdfl[[c]][i,t]-cdfl[[c-1]][i,t]))*(pdfl[[c-1]][i,t]*(tau2[c-1]-w[i,t]))
        taalph[i,t,(q-1)] <- (1/(cdfl[[c]][i,t]-cdfl[[c-1]][i,t])^2)*(pdfl[[c-1]][i,t]-pdfl[[c]][i,t])*pdfl[[c-1]][i,t] -
          (1/(cdfl[[c]][i,t]-cdfl[[c-1]][i,t]))*pdfl[[c-1]][i,t]*(tau2[c-1]-w[i,t])
        htabeta[,(q-1)] <- htabeta[,(q-1)] + taalph[i,t,(q-1)]*x[i,t,]
        # second limit
        ta[i,t,q] <- (1/(cdfl[[c]][i,t]-cdfl[[c-1]][i,t]))*(pdfl[[c]][i,t]) 
        ta2[i,t,q] <- -(1/(cdfl[[c]][i,t]-cdfl[[c-1]][i,t])^2)*(pdfl[[c]][i,t]^2) - (1/(cdfl[[c]][i,t]-cdfl[[c-1]][i,t]))*(pdfl[[c]][i,t]*(tau2[c]-w[i,t]))
        taalph[i,t,q] <- -(1/(cdfl[[c]][i,t]-cdfl[[c-1]][i,t])^2)*(pdfl[[c-1]][i,t]-pdfl[[c]][i,t])*pdfl[[c]][i,t] +
          (1/(cdfl[[c]][i,t]-cdfl[[c-1]][i,t]))*pdfl[[c]][i,t]*(tau2[c]-w[i,t]) 
        htabeta[,q] <- htabeta[,q] + taalph[i,t,q]*x[i,t,]
        # double derivative limits
        ta12[i,t,(q-1)] <- (1/(cdfl[[c]][i,t]-cdfl[[c-1]][i,t])^2)*pdfl[[c-1]][i,t]*pdfl[[c]][i,t]
      }
      else{
        alph[i,t] <- (1/(cdfl[[c]][i,t]-cdfl[[c-1]][i,t]))*(pdfl[[c-1]][i,t]-pdfl[[c]][i,t])
        alph2[i,t] <- -(1/(cdfl[[c]][i,t]-cdfl[[c-1]][i,t])^2)*((pdfl[[c-1]][i,t]-pdfl[[c]][i,t])^2) + 
          (1/(cdfl[[c]][i,t]-cdfl[[c-1]][i,t]))*(pdfl[[c-1]][i,t]*(tau2[c-1]-w[i,t]))
        gbeta <- gbeta + alph[i,t]*x[i,t,]
        hbeta <- hbeta + alph2[i,t]*(x[i,t,]%*%t(x[i,t,]))
        hbetaalpha[i,] <- hbetaalpha[i,] + alph2[i,t]*x[i,t,]
        # second limit
        ta[i,t,(q-1)] <- - (1/(cdfl[[c]][i,t]-cdfl[[c-1]][i,t]))*(pdfl[[c-1]][i,t])
        ta2[i,t,(q-1)] <- -(1/(cdfl[[c]][i,t]-cdfl[[c-1]][i,t])^2)*(pdfl[[c-1]][i,t]^2) + (1/(cdfl[[c]][i,t]-cdfl[[c-1]][i,t]))*(pdfl[[c-1]][i,t]*(tau2[c-1]-w[i,t]))
        taalph[i,t,(q-1)] <- (1/(cdfl[[c]][i,t]-cdfl[[c-1]][i,t])^2)*(pdfl[[c-1]][i,t]-pdfl[[c]][i,t])*pdfl[[c-1]][i,t] -
          (1/(cdfl[[c]][i,t]-cdfl[[c-1]][i,t]))*pdfl[[c-1]][i,t]*(tau2[c-1]-w[i,t])
        htabeta[,(q-1)] <- htabeta[,(q-1)] + taalph[i,t,(q-1)]*x[i,t,]
      }
    }
  }
  galpha <- rowSums(alph)
  hii <- rowSums(alph2)
  
  # derivatives to tau
  gtau <- c(sum(ta[,,1]),sum(ta[,,2]))
  gtau <- eval(parse(text = paste0("c(",paste0("sum(ta[,,",1:p,"])",collapse=","),")")))
  htau <- matrix(0,p,p)
  for(l in 1:(p-1)){
    htau[l:(l+1),l:(l+1)] <- matrix(c(sum(ta2[,,l]),sum(ta12[,,l]),sum(ta12[,,l]),sum(ta2[,,(l+1)])))
  }

  # cross derivatives
  htaualpha <- eval(parse(text = paste0("cbind(",paste0("rowSums(taalph[,,",1:p,"])",collapse=","),")")))
  hgammai <- cbind(hbetaalpha,htaualpha)
  htaubeta <-  eval(parse(text = paste0("cbind(",paste0("htabeta[,",1:p,"]",collapse=","),")")))
  
  ggamma <- c(gbeta,gtau)
  Hgamma <- matrix(0,length(gamma),length(gamma))
  Hgamma[1:k,1:k] <- hbeta
  Hgamma[(k+1):(k+p),1:k] <- htaubeta
  Hgamma[1:k,(k+1):(k+p)] <- t(htaubeta)
  Hgamma[(k+1):(k+p),(k+1):(k+p)] <- htau
  
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
  rownames(gamma) <- c(paste0("beta",1:k),paste0("tau",1:p))
  alpha <- alpha + deltaalpha
  beta <- gamma[1:k]
  tau <- gamma[(k+1):(k+p)]
  dif <- abs(gammaold[,ncol(gammaold)]-gamma)
  show(proc.time() - ptm)
}

