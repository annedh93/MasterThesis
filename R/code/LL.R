steplength <- function(ro,c,stepbeg,par,x,delta,der){
  step <- stepbeg
  while(LL(par + step*delta,x,y) <= (LL(par,x,y) + c*step*t(der)%*%delta)){
    step <- ro*step
  }
  return(step)
}

LL <- function(par,x,y){
  beta <- par[1:k]
  tau <- par[(k+1):(k+p)]
  alpha <- par[(k+p+1):length(par)]
  w <- eval(parse(text=paste(paste0("beta","[",1:k,"]*","x","[,,",1:k,"]",collapse="+")))) + 
    matrix(alpha,nrow=N,ncol=T,byrow=FALSE)
  
  for(l in 1:length(tau2)){
    cdfl[[l]] <- pnorm(tau2[l] - w)
  }
  
  logl <- matrix(0,N,T)
  for(j in 1:4){
    ll <- (y==(j-1))*log(cdfl[[j+1]]-cdfl[[j]])
    ll[(y!=(j-1))] <- 0
    logl <- logl + ll
  }
  return(sum(logl))
}

LLna <- function(par,x,y,N,T,k){

  completecases <- matrix(NA,N,(k+1))
  completecases[,1] <- complete.cases(y)
  for(p in 1:k){
    completecases[,(p+1)] <- complete.cases(x[,,p])
  }
  completecases[(rowSums(completecases)<(k+1)),] <- rep(FALSE,(k+1)) 
  completecases <- completecases[,1]
  
  y <- y - rowMeans(y)
  y <- y[completecases,]
  xdat <- array(NA,dim=c(N,T,k))
  for(p in 1:k){
    xdat[,,p] <- x[,,p] - rowMeans(x[,,p])
  }
  x <- xdat[completecases,,]
  
  beta <- par
  logl <- log( (1/sqrt(2*pi))*exp(-(1/2)*((y-eval(parse(text=paste(paste0("beta","[",1:k,"]*","x","[,,",1:k,"]",collapse="+")))))^2)))

  return(-sum(logl))
}

NR <- function(par,xdat,ydat,N,T,k){
  completecases <- matrix(NA,N,(k+1))
  completecases[,1] <- complete.cases(ydat)
  for(p in 1:k){
    completecases[,(p+1)] <- complete.cases(xdat[,,p])
  }
  completecases[(rowSums(completecases)<(k+1)),] <- rep(FALSE,(k+1)) 
  completecases <- completecases[,1]
  
  y <- ydat - rowMeans(ydat)
  y <- y[completecases,]
  x <- array(NA,dim=c(N,T,k))
  for(p in 1:k){
    x[,,p] <- xdat[,,p] - rowMeans(xdat[,,p])
  }
  x <- x[completecases,,]
  N <- nrow(y)
  T <- 5

  beta <- par
  betaold <- rep(10^6,k)
  while(sum(abs(beta-betaold))>(10^-3)){
    grad <- rep(0,k)
    Hes <- matrix(0,k,k)
    for(i in 1:N){
      for(t in 1:T){
        grad <- grad + y[i,t]*x[i,t,] - x[i,t,]%*%t(x[i,t,])%*%beta 
        Hes <- Hes - x[i,t,]%*%t(x[i,t,])
      }
    }
    betaold <- beta
    beta <- beta - solve(Hes)%*%grad
  }
  return(beta)
}