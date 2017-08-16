# experimental setup 
install.packages("gtools")
library(gtools)

N <- dim(x)[1]
T <- dim(x)[2]
R <- 100

betaor <- c(0.1,0.02,0.1) # change this for different number of parameters
p <- length(betaor)
betaf <- matrix(0,p,R)
for(r in 1:R){
beta <- betaor
u <- matrix(rnorm(N*T,0,1),N,T)
#X <- array(rnorm(N*T*p,0,1),dim=c(N,T,p))
X <- x
a <- sqrt(T)*rowMeans(X[,,1],na.rm=TRUE) + sqrt(T)*rowMeans(u)
epsilon <- matrix(rlogis(N*T,0,1),N,T)
ystar <- matrix(a,nrow=N,ncol=T,byrow=FALSE) + eval(parse(text=paste(paste0("beta","[",1:p,"]*","X","[,,",1:p,"]",collapse="+")))) + epsilon
y <- (ystar > 150) + (ystar > 170) + (ystar > 200)
c <- length(table(y))-1
#sample1 <- sample(N,ceiling((3/4)*N))
#sample2 <- sample(N,ceiling((4/5)*N))
#X[-sample1,1,] <- NA
#X[-sample2,1:2,] <- NA
#y[-sample1,1] <- NA
#y[-sample2,1:2] <- NA


# with real data
N <- nrow(y)
T <- ncol(y)
p <- dim(x)[3]
X <- x
c <- 4

# dichotomized ordered logit model
betak <- matrix(0,p,c)
gbeta <- array(0,dim=c(p,N,c))

for(k in 1:c){
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

mom1 <- matrix(0,p*c,p*c)
mom2 <- matrix(0,p*c,p*c)
for(i in 1:N){
  for(a in 1:c){
    for(b in 1:c){
      mom1[(a*p-p+1):(a*p),(a*p-p+1):(a*p)] <- mom1[(a*p-p+1):(a*p),(a*p-p+1):(a*p)] + gbeta[,i,a]%*%t(gbeta[,i,a])
      mom2[(a*p-p+1):(a*p),(b*p-p+1):(b*p)] <- mom2[(a*p-p+1):(a*p),(b*p-p+1):(b*p)] + gbeta[,i,a]%*%t(gbeta[,i,b])
    }
  }
}

omega <- solve((1/N)*mom1)%*%((1/N)*mom2)%*%solve((1/N)*mom1)
H <- matrix(rep(diag(1,p),c),p*c,p,byrow=TRUE)
betaf <- solve(t(H)%*%solve(omega)%*%H)%*%t(H)%*%solve(omega)%*%c(betak)
rownames(betaf) <- regressors
var <- solve(t(H)%*%solve(omega)%*%H)
}