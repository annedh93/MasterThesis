N <- 1000
k <- 5
X <- matrix(rnorm(k*N,0,1),N,k)
beta <- c(1,0.2,-1,-2,3)
y <- X%*%beta + rnorm(N,0,1)

b <- rep(0,k)
LLOLS <- function(b){
  return(-sum(log(dnorm(y-X%*%b))))
}
par <- optim(b,LLOLS,method="BFGS",hessian=TRUE)
bML <- par$par
HessianML <- par$hessian

dif <- rep(10,k)
eps <- rep(10^-3,k)
while(sum(dif>eps)>0){
  gbeta <- -t(X)%*%y + t(X)%*%X%*%b
  Hbeta <- t(X)%*%X
  delta <- -solve(Hbeta)%*%gbeta
  bold <- b
  b <- b + delta
  dif <- abs(bold-b)
}

sd <- sqrt(diag(solve(Hbeta))) # no - because this is Hessian of - log likelihood
z <- b/sd
pvalue <- pnorm(abs(z),lower.tail=FALSE)
star <- rep("",length(b))
star[pvalue < 0.05] <- "*"
star[pvalue < 0.01] <- "**"
star[pvalue < 0.001] <- "***"
show(as.data.frame(cbind(b,sd,z,pvalue,star)))

# OLS
OLS <- lm(y~X)
