# with real data
N <- nrow(y)
T <- ncol(y)
p <- dim(x)[3]
X <- x
c <- 4

# dichotomized ordered logit model
betak <- matrix(0,p,c)
gbeta <- array(0,dim=c(p,N,c))

for(k in 2:c){
  dich <- (y >= k)
  beta <- rep(0,p)
  betaold <- NULL
  
  dif <- rep(100,p)
  eps <- rep(10^-3,p)
  
  while(sum(dif>eps)>0){
    Hbeta <- matrix(0,p,p)
    
    for(i in 1:N){
      if(sum(is.na(dich[i,]))!=T){
        ai <- sum(dich[i,],na.rm=TRUE)
        Ti <- sum(!is.na(dich[i,]))
        index <- which(!is.na(dich[i,]))
        
        # creating set Bi
        comb <- factorial(Ti)/(factorial(ai)*factorial(Ti-ai)) 
        Bi <- matrix(0,T,comb)
        if(ai > 0){
          ind <- t(combinations(Ti,ai,1:Ti))
          for(l in 1:comb){
            Bi[index[ind[,l]],l] <- 1 
          }
        }
        
        num1 <- rep(0,Ti)
        num2 <- rep(0,p)
        denom <- 0
        sum <- matrix(0,p,p)
        for(l in 1:comb){
          num1 <- num1 + Bi[index,l]%*%exp(t(Bi[index,l])%*%X[i,index,]%*%beta)
          num2 <- num2 + t(exp(t(Bi[index,l])%*%X[i,index,]%*%beta)%*%t(Bi[index,l])%*%X[i,index,])
          denom <- denom + exp(t(Bi[index,l])%*%X[i,index,]%*%beta)
        }
        for(m in 1:comb){
          sum <- sum + as.vector(exp(t(Bi[index,l])%*%X[i,index,]%*%beta)/denom)*(t(X[i,index,])%*%Bi[index,l] - 
                                                                                    num2/as.vector(denom))%*%t(t(X[i,index,])%*%Bi[index,l] - num2/as.vector(denom))
        }
        gbeta[,i,k] <- t(X[i,index,])%*%(as.matrix(dich[i,index]) - num1/as.vector(denom))
        Hbeta <- Hbeta - sum
      }
    }
    betaold <- cbind(betaold,beta)
    deltabeta <- -solve(Hbeta)%*%rowSums(gbeta[,,k])
    beta <- beta + deltabeta
    dif <- abs(betaold[,ncol(betaold)]-beta)
    show(k)
  }
  betak[,k] <- beta

}

mom1 <- matrix(0,p*(c-1),p*(c-1))
mom2 <- matrix(0,p*(c-1),p*(c-1))
for(i in 1:N){
  for(a in 1:(c-1)){
    for(b in 1:(c-1)){
      mom1[(a*p-p+1):(a*p),(a*p-p+1):(a*p)] <- mom1[(a*p-p+1):(a*p),(a*p-p+1):(a*p)] + gbeta[,i,(a+1)]%*%t(gbeta[,i,(a+1)])
      mom2[(a*p-p+1):(a*p),(b*p-p+1):(b*p)] <- mom2[(a*p-p+1):(a*p),(b*p-p+1):(b*p)] + gbeta[,i,(a+1)]%*%t(gbeta[,i,(b+1)])
    }
  }
}

mom1 <- mom1/3
omega <- solve((1/N)*mom1)%*%((1/N)*mom2)%*%solve((1/N)*mom1)
I <- matrix(rep(diag(1,p),(c-1)),p*(c-1),p,byrow=TRUE)
betaf <- solve(t(I)%*%solve(omega)%*%I)%*%t(I)%*%solve(omega)%*%c(betak[,2:4])
rownames(betaf) <- regressors
var <- solve(t(I)%*%solve(omega)%*%I)
