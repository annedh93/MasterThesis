# implementing Newton Rhapson algorithm with ordered probit model

ptm <- proc.time()
# simulate data
N <- 100
T <- 5

beta <- c(1) # change this for different number of parameters
o <- length(beta)
x <- array(rnorm(N*T*o,0,1),dim = c(N,T,o))
a <- rnorm(N,0,1)
meanx_i <- rowMeans(x[,,1])
alpha <- meanx_i + a
w <- matrix(alpha,nrow=N,ncol=T,byrow=FALSE) + eval(parse(text=paste(paste0("beta","[",1:o,"]*","x","[,,",1:o,"]",collapse="+"))))
epsilon <- matrix(rnorm(N*T,0,1),nrow=N)
y <- (w + epsilon > 0) + (w + epsilon > 2) + (w + epsilon > 3) # change this for different (number of) limits

# ordered probit
beta1 <- rep(0,o)
beta2 <- rep(0,o) 
tau <- c(1,2.5) # make sure initialisation is separated enough # change this for different (number of) limits
p <- length(tau)
gamma <- c(beta1,beta2,tau)
eps <- 10^(-3)*rep(1,length(gamma))
dif <- 100*rep(1,length(gamma))
gammaold <- NULL
rowmeans <- eval(parse(text = paste0("cbind(",paste0("rowMeans(x[,,",1:o,"])",collapse =","),  ")")))

while(sum(dif > eps)>0){

  cdfl <- list()
  pdfl <- list()
  tau2 <- c(-Inf,0,tau,Inf) # set one to 0 for identification
  w <- eval(parse(text=paste(paste0("beta1","[",1:o,"]*","x","[,,",1:o,"]",collapse="+")))) + 
    matrix(eval(parse(text=paste(paste0("beta2","[",1:o,"]*","rowMeans(x","[,,",1:o,"])",collapse="+")))),nrow=N,ncol=T,byrow=FALSE)
  for(l in 1:length(tau2)){
    cdfl[[l]] <- pnorm(tau2[l] - w)
    pdfl[[l]] <- dnorm(tau2[l] - w)    
  }
  
  # derivatives to individual specific alpha and beta
  k <- length(gamma) - length(tau)
  gbeta <- matrix(0,k,1)
  hbeta <- matrix(0,k,k)
  ta <- array(0,dim=c(N,T,p))
  ta2 <- array(0,dim=c(N,T,p))
  gtau <- matrix(0,p,1)
  htau <- matrix(0,p,p)
  htabeta <- matrix(0,k,p)
  if(p > 1){
    ta12 <- array(0,dim=c(N,T,(p-1)))    
  }
  for(i in 1:N){
    for(t in 1:T){
      var <- c(x[i,t,],rowmeans[i,])
      if(y[i,t]==0){
        c <- y[i,t] + 2
        
        alph <- (1/(cdfl[[c]][i,t]-cdfl[[c-1]][i,t]))*(pdfl[[c-1]][i,t]-pdfl[[c]][i,t])
        alph2 <- -(1/(cdfl[[2]][i,t]-cdfl[[1]][i,t])^2)*((pdfl[[1]][i,t]-pdfl[[2]][i,t])^2) + 
          (1/(cdfl[[c]][i,t]-cdfl[[c-1]][i,t]))*(-pdfl[[c]][i,t]*(tau2[c]-w[i,t]))
        gbeta <- gbeta + alph*var
        hbeta <- hbeta + alph2*(var%*%t(var))
      }
      else if(y[i,t]==1){
        c <- y[i,t] + 2
        q <- y[i,t]
        
        alph <- (1/(cdfl[[c]][i,t]-cdfl[[c-1]][i,t]))*(pdfl[[c-1]][i,t]-pdfl[[c]][i,t])
        alph2 <- -(1/(cdfl[[c]][i,t]-cdfl[[c-1]][i,t])^2)*((pdfl[[c-1]][i,t]-pdfl[[c]][i,t])^2) + 
          (1/(cdfl[[c]][i,t]-cdfl[[c-1]][i,t]))*(pdfl[[c-1]][i,t]*(tau2[c-1]-w[i,t]) - pdfl[[c]][i,t]*(tau2[c]-w[i,t]))
        gbeta <- gbeta + alph*var
        hbeta <- hbeta + alph2*(var%*%t(var))
        # first limit
        ta[i,t,q] <- (1/(cdfl[[c]][i,t]-cdfl[[c-1]][i,t]))*(pdfl[[c]][i,t]) 
        ta2[i,t,q] <- -(1/(cdfl[[c]][i,t]-cdfl[[c-1]][i,t])^2)*(pdfl[[c]][i,t]^2) - (1/(cdfl[[c]][i,t]-cdfl[[c-1]][i,t]))*(pdfl[[c]][i,t]*(tau2[c]-w[i,t]))
        taalph <- -(1/(cdfl[[c]][i,t]-cdfl[[c-1]][i,t])^2)*(pdfl[[c-1]][i,t]-pdfl[[c]][i,t])*pdfl[[c]][i,t] +
          (1/(cdfl[[c]][i,t]-cdfl[[c-1]][i,t]))*pdfl[[c]][i,t]*(tau2[c]-w[i,t]) 
        htabeta[,q] <- htabeta[,q] + taalph*var
      }
      else if( (y[i,t]>1) & (y[i,t]<(p+1)) ){
        c <- y[i,t] + 2
        q <- y[i,t]
        
        alph <- (1/(cdfl[[c]][i,t]-cdfl[[c-1]][i,t]))*(pdfl[[c-1]][i,t]-pdfl[[c]][i,t])
        alph2 <- -(1/(cdfl[[c]][i,t]-cdfl[[c-1]][i,t])^2)*((pdfl[[c-1]][i,t]-pdfl[[c]][i,t])^2) + 
          (1/(cdfl[[c]][i,t]-cdfl[[c-1]][i,t]))*(pdfl[[c-1]][i,t]*(tau2[c-1]-w[i,t]) - pdfl[[c]][i,t]*(tau2[c]-w[i,t]))
        gbeta <- gbeta + alph*var
        hbeta <- hbeta + alph2*(var%*%t(var))
        # first limit
        ta[i,t,(q-1)] <- - (1/(cdfl[[c]][i,t]-cdfl[[c-1]][i,t]))*(pdfl[[c-1]][i,t])
        ta2[i,t,(q-1)] <- -(1/(cdfl[[c]][i,t]-cdfl[[c-1]][i,t])^2)*(pdfl[[c-1]][i,t]^2) + (1/(cdfl[[c]][i,t]-cdfl[[c-1]][i,t]))*(pdfl[[c-1]][i,t]*(tau2[c-1]-w[i,t]))
        taalph <- (1/(cdfl[[c]][i,t]-cdfl[[c-1]][i,t])^2)*(pdfl[[c-1]][i,t]-pdfl[[c]][i,t])*pdfl[[c-1]][i,t] -
          (1/(cdfl[[c]][i,t]-cdfl[[c-1]][i,t]))*pdfl[[c-1]][i,t]*(tau2[c-1]-w[i,t])
        htabeta[,(q-1)] <- htabeta[,(q-1)] + taalph*var
        # second limit
        ta[i,t,q] <- (1/(cdfl[[c]][i,t]-cdfl[[c-1]][i,t]))*(pdfl[[c]][i,t]) 
        ta2[i,t,q] <- -(1/(cdfl[[c]][i,t]-cdfl[[c-1]][i,t])^2)*(pdfl[[c]][i,t]^2) - (1/(cdfl[[c]][i,t]-cdfl[[c-1]][i,t]))*(pdfl[[c]][i,t]*(tau2[c]-w[i,t]))
        taalph <- -(1/(cdfl[[c]][i,t]-cdfl[[c-1]][i,t])^2)*(pdfl[[c-1]][i,t]-pdfl[[c]][i,t])*pdfl[[c]][i,t] +
          (1/(cdfl[[c]][i,t]-cdfl[[c-1]][i,t]))*pdfl[[c]][i,t]*(tau2[c]-w[i,t]) 
        htabeta[,q] <- htabeta[,q] + taalph*var
        # double derivative limits
        ta12[i,t,(q-1)] <- (1/(cdfl[[c]][i,t]-cdfl[[c-1]][i,t])^2)*pdfl[[c-1]][i,t]*pdfl[[c]][i,t]
      }
      else{
        c <- y[i,t] + 2
        q <- y[i,t]
        
        alph <- (1/(cdfl[[c]][i,t]-cdfl[[c-1]][i,t]))*(pdfl[[c-1]][i,t]-pdfl[[c]][i,t])
        alph2 <- -(1/(cdfl[[c]][i,t]-cdfl[[c-1]][i,t])^2)*((pdfl[[c-1]][i,t]-pdfl[[c]][i,t])^2) + 
          (1/(cdfl[[c]][i,t]-cdfl[[c-1]][i,t]))*(pdfl[[c-1]][i,t]*(tau2[c-1]-w[i,t]))
        gbeta <- gbeta + alph*var
        hbeta <- hbeta + alph2*(var%*%t(var))
        # second limit
        ta[i,t,(q-1)] <- - (1/(cdfl[[c]][i,t]-cdfl[[c-1]][i,t]))*(pdfl[[c-1]][i,t])
        ta2[i,t,(q-1)] <- -(1/(cdfl[[c]][i,t]-cdfl[[c-1]][i,t])^2)*(pdfl[[c-1]][i,t]^2) + (1/(cdfl[[c]][i,t]-cdfl[[c-1]][i,t]))*(pdfl[[c-1]][i,t]*(tau2[c-1]-w[i,t]))
        taalph <- (1/(cdfl[[c]][i,t]-cdfl[[c-1]][i,t])^2)*(pdfl[[c-1]][i,t]-pdfl[[c]][i,t])*pdfl[[c-1]][i,t] -
          (1/(cdfl[[c]][i,t]-cdfl[[c-1]][i,t]))*pdfl[[c-1]][i,t]*(tau2[c-1]-w[i,t])
        htabeta[,(q-1)] <- htabeta[,(q-1)] + taalph*var
      }
    }
  }

  # derivatives to tau
  gtau <- eval(parse(text = paste0("c(",paste0("sum(ta[,,",1:p,"])",collapse=","),")")))
  htau <- matrix(0,p,p)
  for(l in 1:(p-1)){
    htau[l:(l+1),l:(l+1)] <- matrix(c(sum(ta2[,,l]),sum(ta12[,,l]),sum(ta12[,,l]),sum(ta2[,,(l+1)])))
  }
  
  # cross derivatives
  htaubeta <-  eval(parse(text = paste0("cbind(",paste0("htabeta[,",1:p,"]",collapse=","),")")))
  
  ggamma <- c(gbeta,gtau)
  Hgamma <- matrix(0,length(gamma),length(gamma))
  Hgamma[1:k,1:k] <- hbeta
  Hgamma[(k+1):(k+p),1:k] <- htaubeta
  Hgamma[1:k,(k+1):(k+p)] <- t(htaubeta)
  Hgamma[(k+1):(k+p),(k+1):(k+p)] <- htau
  
  # Newton Rhapson step
  deltagamma <- -solve(Hgamma)%*%ggamma
  gammaold <- cbind(gammaold,gamma)
  gamma <- gamma + 0.7*deltagamma
  rownames(gamma) <- c(paste0("beta1",1:o),paste0("beta2",1:o),paste0("tau",1:p))
  beta1 <- gamma[1:o]
  beta2 <- gamma[(o+1):(o*2)]
  tau <- gamma[(k+1):(k+p)]
  dif <- abs(gammaold[,ncol(gammaold)]-gamma)
  show(proc.time() - ptm)
}
gamma
