# implementing Newton Rhapson algorithm with ordered probit model

# simulate data
N <- 100
T <- 8

beta <- c(1,-1)
x <- array(rnorm(N*T*2,0,1),dim = c(N,T,2))
a <- rnorm(N,0,1)
meanx_i <- rowMeans(x1)
alpha <- sqrt(T)*meanx_i + a
w <- matrix(alpha,nrow=N,ncol=T,byrow=FALSE) + beta[1]*x[,,1] + beta[2]*x[,,2]  
epsilon <- matrix(rnorm(N*T,0,1),nrow=N)
y <- (w + epsilon > 0) + (w + epsilon > 3)
ind0 <- which(y==0,arr.ind=TRUE)
ind1 <- which(y==1,arr.ind=TRUE)
ind2 <- which(y==2,arr.ind=TRUE)

# ordered probit
beta <- c(1,0)
tau <- c(0,1)
gamma <- c(beta,tau)
alpha <- rep(0,N)
eps <- 10^(-3)*rep(1,length(gamma))
dif <- 100*rep(1,length(gamma))
gammaold <- NULL

while(sum(dif > eps)>0){
  sum1 <- matrix(0,length(gamma),length(gamma))
  sum2 <- matrix(0,length(gamma))
  
  cdfl <- list()
  pdfl <- list()
  tau2 <- c(-Inf,tau,Inf)
  w <- x[,,1]*beta[1] + x[,,2]*beta[2] + matrix(alpha,nrow=N,ncol=T,byrow=FALSE)
  for(l in 1:length(tau2)){
    cdfl[[l]] <- pnorm(tau2[l] - w)
    pdfl[[l]] <- dnorm(tau2[l] - w)    
  }
  tau2 <- c(0,tau,0)
  
  alph <- matrix(0,N,T)
  alph2 <- matrix(0,N,T)
  ta <- list(matrix(0,N,T),matrix(0,N,T))
  ta2 <- list(matrix(0,N,T),matrix(0,N,T))
  taalph <- list(matrix(0,N,T),matrix(0,N,T))
  ta12 <- matrix(0,N,T)
  # derivatives to individual specific alpha and beta
  for(c in 0:2){
    ind <- c+2
    # calculate derivatives to alpha
    alph[y==c] <- ((1/(cdfl[[ind]]-cdfl[[ind-1]]))*(pdfl[[ind-1]]-pdfl[[ind]]))[y==c]
    alph2[y==c] <- ((1/(cdfl[[ind]]-cdfl[[ind-1]]))*(pdfl[[ind-1]]*(tau2[ind-1]-w)-pdfl[[ind]]*(tau2[ind]-w)) - 
                   (1/(cdfl[[ind]]-cdfl[[ind-1]])^2)*((pdfl[[ind-1]]-pdfl[[ind]])^2))[y==c] 
  }
  
  for(c in 0:1){
    ind <- c+2
    # calculate derivatives to tau
    ta[[c+1]][y==c] <- ((1/(cdfl[[ind]]-cdfl[[ind-1]]))*(pdfl[[ind]]))[y==c] 
    ta[[c+1]][y==(c+1)] <- (-(1/(cdfl[[ind+1]]-cdfl[[ind]]))*(pdfl[[ind]]))[y==(c+1)]
    ta2[[c+1]][y==c] <- (-(1/(cdfl[[ind]]-cdfl[[ind-1]]))*(pdfl[[ind]]*(tau2[ind]-w))-(1/(cdfl[[ind]]-cdfl[[ind-1]])^2)*(pdfl[[ind]]^2))[y==c] 
    ta2[[c+1]][y==(c+1)] <- (-(1/(cdfl[[ind+1]]-cdfl[[ind]])^2)*(pdfl[[ind]]^2) + (1/(cdfl[[ind+1]]-cdfl[[ind]]))*(pdfl[[ind]]*(tau2[ind]-w)))[y==(c+1)]
    taalph[[c+1]][y==c] <- (-(1/(cdfl[[ind]]-cdfl[[ind-1]])^2)*(pdfl[[ind-1]]-pdfl[[ind]])*pdfl[[ind]] +
                        (1/(cdfl[[ind]]-cdfl[[ind-1]]))*pdfl[[ind]]*(tau2[ind]-w))[y==c] 
    taalph[[c+1]][y==(c+1)] <- ((1/(cdfl[[ind+1]]-cdfl[[ind]])^2)*(pdfl[[ind]]-pdfl[[ind+1]][i,t])*pdfl[[ind]][i,t] -
                        (1/(cdfl[[ind+1]]-cdfl[[ind]]))*pdfl[[ind]]*(tau2[ind]-w))[y==(c+1)]
  }
  ta12 <- ((1/(cdfl[[3]]-cdfl[[2]])^2)*pdfl[[2]]*pdfl[[3]])[y==1] 
  
  gbeta <- matrix(0,2,1)
  hbeta <- matrix(0,2,2)
  hbetaalpha <- matrix(0,N,2)
  htaubeta <- matrix(0,2,2)
  for(i in 1:N){
    for(t in 1:T){
      gbeta <- gbeta + alph[i,t]*x[i,t,]
      hbeta <- hbeta + alph2[i,t]*x[i,t,]%*%t(x[i,t,])
      hbetaalpha[i,] <- hbetaalpha[i,] + alph2[i,t]*x[i,t,]
      htaubeta[,1] <- htaubeta[,1] + taalph[[1]][i,t]*x[i,t,] 
      htaubeta[,2] <- htaubeta[,2] + taalph[[2]][i,t]*x[i,t,] 
    }
  }

  galpha <- rowSums(alph)
  hii <- rowSums(alph2)

  # derivatives to tau
  gtau <- c(sum(ta[[1]]),sum(ta[[2]]))
  htau <- matrix(c(sum(ta2[[1]]),sum(ta12),sum(ta12),sum(ta2[[2]])),nrow = 2,byrow=TRUE)
  
  # cross derivatives
  htaualpha <- cbind(rowSums(taalph[[1]]),rowSums(taalph[[2]]))
  hgammai <- cbind(hbetaalpha,htaualpha)
  
  
  ggamma <- c(gbeta,gtau)
  Hgamma <- matrix(0,length(gamma),length(gamma))
  Hgamma[1:2,1:2] <- hbeta
  Hgamma[3:4,1:2] <- htaubeta
  Hgamma[1:2,3:4] <- t(htaubeta)
  Hgamma[3:4,3:4] <- htau
  
  # Newton Rhapson step
  for(i in 1:N){
    sum1 <- sum1 + (1/hii[i])*hgammai[i,]%*%t(hgammai[i,])
    sum2 <- sum2 + (galpha[i]/hii[i])*hgammai[i,]
  }
  deltagamma <- -solve(Hgamma - sum1)%*%(ggamma - sum2)
  deltaalpha <- -(1/hii)*(galpha + hgammai%*%deltagamma)
  
  gammaold <- cbind(gammaold,gamma)
  gamma <- gamma + deltagamma
  alpha <- alpha + deltaalpha
  beta <- gamma[1]
  tau <- gamma[2:3]
  dif <- abs(gammaold[,ncol(gammaold)]-gamma)
}
