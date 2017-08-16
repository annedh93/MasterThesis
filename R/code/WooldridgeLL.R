# implementing Newton Rhapson algorithm with ordered probit model

ptm <- proc.time()
# simulate data
N <- 100
T <- 5

beta <- c(1,2,0.4) # change this for different number of parameters
o <- length(beta)
x <- array(rnorm(N*T*o,0,1),dim = c(N,T,o))
a <- rnorm(N,0,1)
meanx_i <- rowMeans(x[,,1])
alpha <- sqrt(T)*meanx_i 
w <- matrix(alpha,nrow=N,ncol=T,byrow=FALSE) + eval(parse(text=paste(paste0("beta","[",1:o,"]*","x","[,,",1:o,"]",collapse="+"))))
epsilon <- matrix(rnorm(N*T,0,1),nrow=N)
y <- (w + epsilon > 0) + (w + epsilon > 1.5) + (w + epsilon > 2) + (w + epsilon > 3) # change this for different (number of) limits

# ordered probit
k <- 2*o
x <-  
beta <- rep(0,k)
tau <- c(1.5,2.5,3) # make sure initialisation is separated enough # change this for different (number of) limits
p <- length(tau)
gamma <- c(beta,tau)

par <- optim(gamma,LL,"BFGS",x=x,y=y)
parameters <- par$par
show(parameters)

LL <- function(gamma,x,y){
  beta <- gamma[1:k]
  tau <- gamma[(k+1):(k+p)]
  cdfl <- list()
  pdfl <- list()
  tau2 <- c(-Inf,0,tau,Inf) # set one to 0 for identification
  w <- eval(parse(text=paste(paste0("beta","[",1:k,"]*","x","[,,",1:k,"]",collapse="+"))))
  
  for(l in 1:length(tau2)){
    cdfl[[l]] <- pnorm(tau2[l] - w)
    pdfl[[l]] <- dnorm(tau2[l] - w)    
  }
  
  logl <- matrix(0,N,T)
  for(j in 1:4){
    ll <- (y==(j-1))*log(cdfl[[j+1]]-cdfl[[j]])
    ll[(y!=(j-1))] <- 0
    logl <- logl + ll
  }
  return(-sum(logl))
}

