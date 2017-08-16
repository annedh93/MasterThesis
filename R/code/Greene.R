# implementing Newton Rhapson algorithm with ordered probit model

# simulate data
N <- 100
T <- 3

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
tau <- c(0,1)
alpha <- rep(0,N)
eps <- 10^(-3)*rep(1,3)
dif <- 100*rep(1,3)

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
  phi <- matrix(0,N,T)
  ta <- list()
  for(l in 1:(length(tau2)-1)){
    alph <- alph + (y==(l-1))*(1/(cdfl[[l+1]]-cdfl[[l]]))*(pdfl[[l]]-pdfl[[l+1]])
    if(l==1){
      phi <- phi - (y==(l-1))*(1/(cdfl[[l+1]]-cdfl[[l]])^2)*((pdfl[[l]]-pdfl[[l+1]])^2) + 
      (y==(l-1))*(1/(cdfl[[l+1]]-cdfl[[l]]))*(-pdfl[[l+1]]*(tau2[l+1]-w))
    }
    else if(l==(length(tau2)-1)){
      phi <- phi - (y==(l-1))*(1/(cdfl[[l+1]]-cdfl[[l]])^2)*((pdfl[[l]]-pdfl[[l+1]])^2) + 
      (y==(l-1))*(1/(cdfl[[l+1]]-cdfl[[l]]))*(pdfl[[l]]*(tau2[l]-w))
    }
    else{
      phi <- phi - (y==(l-1))*(1/(cdfl[[l+1]]-cdfl[[l]])^2)*((pdfl[[l]]-pdfl[[l+1]])^2) + 
      (y==(l-1))*(1/(cdfl[[l+1]]-cdfl[[l]]))*(pdfl[[l]]*(tau2[l]-w)-pdfl[[l+1]]*(tau2[l+1]-w))
    }
  }
  galpha <- rowSums(alph)
  hii <- rowSums(phi)
  
  # derivatives to tau
  for(k in 1:(length(tau))){
    ta[[k]] <- (y==(k-1))*(1/(cdfl[[k+1]]-cdfl[[k]]))*(pdfl[[k+1]]) - (y==k)*(1/(cdfl[[k+2]]-cdfl[[k+1]]))*(pdfl[[k+1]])
  }
  
  
  
  delta <- (y==0)*(1/(cdfl[[2]]))*(-pdfl[[2]]) + (y==1)*(1/(cdfl[[3]]-cdfl[[2]]))*(pdfl[[2]]-pdfl[[3]]) + (y==2)*(1/(1-cdfl[[3]]))*(pdfl[[3]])
  phi2 <- -(y==0)*(1/(cdfl[[2]])^2)*(pdfl[[2]])^2 - (y==0)*(1/(cdfl[[2]]))* (tau2[2]-w)*pdfl[[2]] -
          (y==1)*(1/(cdfl[[3]]-cdfl[[2]])^2)*(-pdfl[[3]] + pdfl[[2]])^2 + (y==1)*(1/(cdfl[[3]]-cdfl[[2]]))*(-(tau2[3]-w)*pdfl[[3]] + (tau2[2]-w)*pdfl[[2]]) -
          (y==2)*(1/(1-cdfl[[3]])^2)*(pdfl[[3]])^2 + (y==2)*(1/(1-cdfl[[3]]))*(tau2[3]-w)*pdfl[[3]]
  galpha <- rowSums(delta) 
  hii <- rowSums(phi)
  
  # derivatives to beta and tau
  gbeta <- sum(delta*x)
  gtau0 <- sum((y==0)*(1/(cdfl[[1]]))*(pdfl[[1]]) - (y==1)*(1/(cdfl[[2]]-cdfl[[1]]))*(pdfl[[1]])) 
  gtau1 <- sum((y==1)*(1/(cdfl[[2]]-cdfl[[1]]))*(pdfl[[2]]) - (y==2)*(1/(1-cdfl[[2]]))*(pdfl[[2]]))
  ggamma <- c(gbeta,gtau0,gtau1)
  
  hbeta <- sum(phi*(x^2))
  htau0 <- sum(-(y==0)*(1/(cdfl[[1]])^2)*(pdfl[[1]])^2 - (y==0)*(1/(cdfl[[1]]))*(tau[1]-w)*pdfl[[1]] +
            (y==1)*(1/(cdfl[[2]]-cdfl[[1]])^2)*(pdfl[[1]])^2 + (y==1)*(1/(cdfl[[2]]-cdfl[[1]]))*((tau[1]-w)*pdfl[[1]]))
  htau1 <- sum(-(y==1)*(1/(cdfl[[2]]-cdfl[[1]])^2)*(pdfl[[2]])^2 - (y==1)*(1/(cdfl[[2]]-cdfl[[1]]))*((tau[2]-w)*pdfl[[2]]) +
            (y==2)*(1/(1-cdfl[[2]])^2)*(pdfl[[2]])^2 + (y==2)*(1/(1-cdfl[[2]]))*((tau[2]-w)*pdfl[[2]]))
  
  tau0alpha <- -(y==0)*(1/(cdfl[[1]])^2)*(-pdfl[[1]])*pdfl[[1]] + (y==0)*(1/(cdfl[[1]]))*((tau[1]-w)*pdfl[[1]]) +
                 (y==1)*(1/(cdfl[[2]]-cdfl[[1]])^2)*(pdfl[[1]]-pdfl[[2]])*pdfl[[1]] - (y==1)*(1/(cdfl[[2]]-cdfl[[1]]))*((tau[1]-w)*pdfl[[1]])
  tau1alpha <- -(y==1)*(1/(cdfl[[2]]-cdfl[[1]])^2)*(pdfl[[1]]-pdfl[[2]])*pdfl[[2]] + (y==1)*(1/(cdfl[[2]]-cdfl[[1]]))*((tau[2]-w)*pdfl[[2]]) +
                 (y==2)*(1/(1-cdfl[[2]])^2)*(pdfl[[2]])*pdfl[[2]] - (y==2)*(1/(1-cdfl[[2]]))*((tau[2]-w)*pdfl[[2]])
  tau0tau1 <- -(y==1)*(1/(cdfl[[2]]-cdfl[[1]])^2)*(pdfl[[1]])*(pdfl[[2]])
  
  hbetaalpha <- rowSums(phi*x)
  htau0alpha <- rowSums(tau0alpha)
  htau1alpha <- rowSums(tau1alpha)
  hgammai <- cbind(hbetaalpha,htau0alpha,htau1alpha)
  
  htau0beta <- sum(tau0alpha*x)
  htau1beta <- sum(tau1alpha*x)
  htau0tau1 <- sum(tau0tau1)
  
  Hgamma <- cbind(c(hbeta,htau0beta,htau1beta),c(htau0beta,htau0,htau0tau1),c(htau1beta,htau0tau1,htau1))
  
  # Newton Rhapson step
  for(i in 1:N){
    sum1 <- sum1 + (1/hii[i])*hgammai[i,]%*%t(hgammai[i,])
    sum2 <- sum2 + (galpha[i]/hii[i])*hgammai[i,]
  }
  deltagamma <- -solve(Hgamma - sum1)%*%(ggamma - sum2)
  deltaalpha <- -(1/hii)*(galpha + hgammai%*%deltagamma)
  
  gammaold <- gamma
  gamma <- gamma + deltagamma
  alpha <- alpha + deltaalpha
  beta <- gamma[1]
  tau <- gamma[2:3]
  dif <- abs(gammaold-gamma)
}