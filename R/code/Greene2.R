# implementing Newton Rhapson algorithm with ordered probit model

# simulate data
N <- 100
T <- 300

beta <- 1
x <- matrix(rnorm(N*T,0,1),nrow=N)
a <- rnorm(N,0,1)
meanx_i <- rowMeans(x)
alpha <- sqrt(T)*meanx_i + a
w <- matrix(alpha,nrow=N,ncol=T,byrow=FALSE) + beta*x 
epsilon <- matrix(rnorm(N*T,0,1),nrow=N)
y <- (w + epsilon > 0) + (w + epsilon > 3)

# ordered probit
beta <- 1.2
tau <- c(0,2)
gamma <- c(beta,tau)
alpha <- rep(0,N)
eps <- 10^(-3)*rep(1,3)
dif <- 100*rep(1,3)
gammaold <- NULL

while(sum(dif > eps)>0){
  sum1 <- matrix(0,3,3)
  sum2 <- matrix(0,3)
  
  cdfl <- list()
  pdfl <- list()
  tau2 <- c(-Inf,tau,Inf)
  w <- x*beta + matrix(alpha,nrow=N,ncol=T,byrow=FALSE)
  for(l in 1:length(tau2)){
    cdfl[[l]] <- pnorm(tau2[l] - w)
    pdfl[[l]] <- dnorm(tau2[l] - w)    
  }

  # derivatives to individual specific alpha and beta
  alph <- matrix(0,N,T)
  alph2 <- matrix(0,N,T)
  ta <- list(matrix(0,N,T),matrix(0,N,T))
  ta2 <- list(matrix(0,N,T),matrix(0,N,T))
  taalph <- list(matrix(0,N,T),matrix(0,N,T))
  ta12 <- matrix(0,N,T)
  for(i in 1:N){
    for(t in 1:T){
      if(y[i,t]==0){
        alph[i,t] <- (1/(cdfl[[2]][i,t]-cdfl[[1]][i,t]))*(pdfl[[1]][i,t]-pdfl[[2]][i,t])
        alph2 <- -(1/(cdfl[[2]]-cdfl[[1]])^2)*((pdfl[[1]]-pdfl[[2]])^2) + 
        (1/(cdfl[[2]]-cdfl[[1]]))*(-pdfl[[2]]*(tau2[2]-w))
        ta[[1]][i,t] <- (1/(cdfl[[2]][i,t]-cdfl[[1]][i,t]))*(pdfl[[2]][i,t]) 
        ta2[[1]][i,t] <- -(1/(cdfl[[2]][i,t]-cdfl[[1]][i,t])^2)*(pdfl[[2]][i,t]^2) - (1/(cdfl[[2]][i,t]-cdfl[[1]][i,t]))*(pdfl[[2]][i,t]*(tau2[2]-w[i,t]))
        taalph[[1]][i,t] <- -(1/(cdfl[[2]][i,t]-cdfl[[1]][i,t])^2)*(pdfl[[1]][i,t]-pdfl[[2]][i,t])*pdfl[[2]][i,t] +
          (1/(cdfl[[2]][i,t]-cdfl[[1]][i,t]))*pdfl[[2]][i,t]*(tau2[2]-w[i,t]) 
      }
      else if(y[i,t]==1){
        alph[i,t] <- (1/(cdfl[[3]][i,t]-cdfl[[2]][i,t]))*(pdfl[[2]][i,t]-pdfl[[3]][i,t])
        alph2 <- -(1/(cdfl[[3]]-cdfl[[2]])^2)*((pdfl[[2]]-pdfl[[3]])^2) + 
        (1/(cdfl[[3]]-cdfl[[2]]))*(pdfl[[2]]*(tau2[2]-w)-pdfl[[3]]*(tau2[3]-w))
        ta[[1]][i,t] <- - (1/(cdfl[[3]][i,t]-cdfl[[2]][i,t]))*(pdfl[[2]][i,t])
        ta2[[1]][i,t] <- -(1/(cdfl[[3]][i,t]-cdfl[[2]][i,t])^2)*(pdfl[[2]][i,t]^2) + (1/(cdfl[[3]][i,t]-cdfl[[2]][i,t]))*(pdfl[[2]][i,t]*(tau2[2]-w[i,t]))
        taalph[[1]][i,t] <- (1/(cdfl[[3]][i,t]-cdfl[[2]][i,t])^2)*(pdfl[[2]][i,t]-pdfl[[3]][i,t])*pdfl[[2]][i,t] -
          (1/(cdfl[[3]][i,t]-cdfl[[2]][i,t]))*pdfl[[2]][i,t]*(tau2[2]-w[i,t])
        ta[[2]][i,t] <- (1/(cdfl[[3]][i,t]-cdfl[[2]][i,t]))*(pdfl[[3]][i,t]) 
        ta2[[2]][i,t] <- -(1/(cdfl[[3]][i,t]-cdfl[[2]][i,t])^2)*(pdfl[[3]][i,t]^2) - (1/(cdfl[[3]][i,t]-cdfl[[2]][i,t]))*(pdfl[[3]][i,t]*(tau2[3]-w[i,t]))
        taalph[[2]][i,t] <- -(1/(cdfl[[3]][i,t]-cdfl[[2]][i,t])^2)*(pdfl[[2]][i,t]-pdfl[[3]][i,t])*pdfl[[3]][i,t] +
          (1/(cdfl[[3]][i,t]-cdfl[[2]][i,t]))*pdfl[[3]][i,t]*(tau2[3]-w[i,t]) 
        ta12[i,t] <- (1/(cdfl[[3]][i,t]-cdfl[[2]][i,t])^2)*pdfl[[2]][i,t]*pdfl[[3]][i,t]
      }
      else{
        alph[i,t] <- (1/(cdfl[[4]][i,t]-cdfl[[3]][i,t]))*(pdfl[[3]][i,t]-pdfl[[4]][i,t])
        alph2 <- -(1/(cdfl[[4]]-cdfl[[3]])^2)*((pdfl[[3]]-pdfl[[4]])^2) + 
        (1/(cdfl[[4]]-cdfl[[3]]))*(pdfl[[3]]*(tau2[3]-w))
        ta[[2]][i,t] <- - (1/(cdfl[[4]][i,t]-cdfl[[3]][i,t]))*(pdfl[[3]][i,t])
        ta2[[2]][i,t] <- -(1/(cdfl[[4]][i,t]-cdfl[[3]][i,t])^2)*(pdfl[[3]][i,t]^2) + (1/(cdfl[[4]][i,t]-cdfl[[3]][i,t]))*(pdfl[[3]][i,t]*(tau2[3]-w[i,t]))
        taalph[[2]][i,t] <- (1/(cdfl[[4]][i,t]-cdfl[[3]][i,t])^2)*(pdfl[[3]][i,t]-pdfl[[4]][i,t])*pdfl[[3]][i,t] -
          (1/(cdfl[[4]][i,t]-cdfl[[3]][i,t]))*pdfl[[3]][i,t]*(tau2[3]-w[i,t])
      }
    }
  }
  galpha <- rowSums(alph)
  hii <- rowSums(alph2)
  gbeta <- sum(alph*x)
  hbeta <- sum(alph2*(x^2))
  
  # derivatives to tau
  gtau <- c(sum(ta[[1]]),sum(ta[[2]]))
  htau <- matrix(c(sum(ta2[[1]]),sum(ta12),sum(ta12),sum(ta2[[2]])),nrow = 2,byrow=TRUE)

  # cross derivatives
  hbetaalpha <- rowSums(alph2*x)
  htaualpha <- cbind(rowSums(taalph[[1]]),rowSums(taalph[[2]]))
  htaubeta <- c(sum(taalph[[1]]*x),sum(taalph[[2]]*x))
  
  hgammai <- cbind(hbetaalpha,htaualpha)
  
  
  ggamma <- c(gbeta,gtau)
  Hgamma <- matrix(0,length(gamma),length(gamma))
  Hgamma[,1] <- c(hbeta,htaubeta)
  Hgamma[1,] <- c(hbeta,htaubeta)
  Hgamma[2:length(gamma),2:length(gamma)] <- htau
  
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
