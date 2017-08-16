# with real data
N <- nrow(y)
T <- ncol(y)
p <- dim(x)[3]
X <- x
c <- 5

# dichotomized ordered logit model
gbeta <- array(0,dim=c(p,N,c))
Hbeta <- array(0,dim=c(p,p,c))
beta <- rep(0,p)
betaold <- NULL

dif <- rep(100,p)
eps <- rep(10^-3,p)

while(sum(dif>eps)>0){
  for(k in 2:c){
    dich <- (y >= k)
    
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
        Hbeta[,,k] <- Hbeta[,,k] - sum
      }
    }
  }
  betaold <- cbind(betaold,beta)
  deltabeta <- -solve(apply(Hbeta,c(1,2),sum))%*%rowSums(gbeta)
  beta <- beta + deltabeta
  dif <- abs(betaold[,ncol(betaold)]-beta)
  show(sum(dif>eps))
}

mom <- matrix(0,p,p)
for(i in 1:N){
  for(a in 1:c){
    for(b in 1:c){
      mom <- mom + (1/N)*gbeta[,i,a]%*%t(gbeta[,i,b])
      
    }
  }
}

H <- apply((1/N)*Hbeta,c(1,2),sum)
var <- solve(H)%*%mom%*%solve(H)
