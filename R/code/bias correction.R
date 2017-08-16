# Implementing the bias correction of Arellano Hahn

v <- matrix(0,N,T);
u <- array(0,dim = c(N,T,(k+p)));
valpha <- matrix(0,N,T);
ualpha <- array(0,dim=c(N,T,(k+p)));
utheta <- array(0,dim=c(N,T,(k+p),(k+p)));
v2alpha <- matrix(0,N,T);
u2alpha <- array(0,dim=c(N,T,(k+p)));
sumuv <- matrix(0,N,(k+p));
sumvv <- matrix(0,N,1);
sumvalpha <- matrix(0,N,1);
sumualpha <- matrix(0,N,(k+p));
sumvualpha <- matrix(0,N,(k+p));
sumvvalpha <- matrix(0,N,1);
sumv2alpha <- matrix(0,N,1);
sumu2alpha <- matrix(0,N,(k+p));
sumutheta <- array(0,dim=c(N,(k+p),(k+p)));
part1b <- matrix(0,N,1);
part2b <- matrix(0,N,(k+p));
part3b <- matrix(0,N,(k+p));

I <- array(0,dim=c(N,(k+p),(k+p)));
b <- matrix(0,N,(k+p));
numerator = 0;
U = matrix(0,N,T);
for(i in 1:N){
  for(t in 1:T){
    
    cdfl <- list()
    pdfl <- list()
    tau2 <- c(-Inf,0,tau,Inf) # set one to 0 for identification
    w <- eval(parse(text=paste(paste0("beta","[",1:k,"]*","x","[,,",1:k,"]",collapse="+")))) + matrix(alpha,nrow=N,ncol=T,byrow=FALSE)
    for(l in 1:length(tau2)){
      cdfl[[l]] <- pnorm(tau2[l] - w)
      pdfl[[l]] <- dnorm(tau2[l] - w)    
    }
    c <- y[i,t] + 2
    q <- y[i,t]
    
    if(y[i,t]==0){
      v[i,t] <- (1/(cdfl[[c]][i,t]-cdfl[[c-1]][i,t]))*(pdfl[[c-1]][i,t]-pdfl[[c]][i,t])
      u[i,t,(1:k)] <- v[i,t]*x[i,t,]
      valpha[i,t] <- -(1/(cdfl[[2]][i,t]-cdfl[[1]][i,t])^2)*((pdfl[[1]][i,t]-pdfl[[2]][i,t])^2) + 
        (1/(cdfl[[c]][i,t]-cdfl[[c-1]][i,t]))*(-pdfl[[c]][i,t]*(tau2[c]-w[i,t]))
      ualpha[i,t,(1:k)] <- valpha[i,t]*x[i,t,]
      utheta[i,t,(1:k),(1:k)] <- valpha[i,t]*x[i,t,]%*%t(x[i,t,])
      v2alpha[i,t] <- 2*(1/(cdfl[[c]][i,t]-cdfl[[c-1]][i,t])^3)*((pdfl[[c-1]][i,t]-pdfl[[c]][i,t])^3) - 
        3*(1/(cdfl[[c]][i,t]-cdfl[[c-1]][i,t])^2)*(pdfl[[c-1]][i,t]-pdfl[[c]][i,t])*(- pdfl[[c]][i,t]*(tau2[c]-w[i,t])) +
        (1/(cdfl[[c]][i,t]-cdfl[[c-1]][i,t]))*(- pdfl[[c]][i,t]*(tau2[c]-w[i,t])^2 + pdfl[[c]][i,t])
      u2alpha[i,t,(1:k)] <- v2alpha[i,t]*x[i,t,]
    }
    else if(y[i,t]==1){
      v[i,t] <- (1/(cdfl[[c]][i,t]-cdfl[[c-1]][i,t]))*(pdfl[[c-1]][i,t]-pdfl[[c]][i,t])
      u[i,t,1:k] <- v[i,t]*x[i,t,]
      u[i,t,(k+q)] <- (1/(cdfl[[c]][i,t]-cdfl[[c-1]][i,t]))*(pdfl[[c]][i,t]) 
      valpha[i,t] <- -(1/(cdfl[[c]][i,t]-cdfl[[c-1]][i,t])^2)*((pdfl[[c-1]][i,t]-pdfl[[c]][i,t])^2) + 
        (1/(cdfl[[c]][i,t]-cdfl[[c-1]][i,t]))*(pdfl[[c-1]][i,t]*(tau2[c-1]-w[i,t]) - pdfl[[c]][i,t]*(tau2[c]-w[i,t]))
      ualpha[i,t,(1:k)] <- valpha[i,t]*x[i,t,]
      ualpha[i,t,(k+q)] <- -(1/(cdfl[[c]][i,t]-cdfl[[c-1]][i,t])^2)*(pdfl[[c-1]][i,t]-pdfl[[c]][i,t])*pdfl[[c]][i,t] +
        (1/(cdfl[[c]][i,t]-cdfl[[c-1]][i,t]))*pdfl[[c]][i,t]*(tau2[c]-w[i,t]) 
      utheta[i,t,(1:k),(1:k)] <- valpha[i,t]*x[i,t,]%*%t(x[i,t,])
      utheta[i,t,(k+q),(k+q)] <- -(1/(cdfl[[c]][i,t]-cdfl[[c-1]][i,t])^2)*(pdfl[[c]][i,t]^2) - (1/(cdfl[[c]][i,t]-cdfl[[c-1]][i,t]))*(pdfl[[c]][i,t]*(tau2[c]-w[i,t]))
      utheta[i,t,(1:k),(k+q)] <- ualpha[i,t,(k+q)]*x[i,t,]
      utheta[i,t,(k+q),(1:k)] <- ualpha[i,t,(k+q)]*x[i,t,]
      v2alpha[i,t] <- 2*(1/(cdfl[[c]][i,t]-cdfl[[c-1]][i,t])^3)*((pdfl[[c-1]][i,t]-pdfl[[c]][i,t])^3) - 
        3*(1/(cdfl[[c]][i,t]-cdfl[[c-1]][i,t])^2)*(pdfl[[c-1]][i,t]-pdfl[[c]][i,t])*(pdfl[[c-1]][i,t]*(tau2[c-1]-w[i,t]) - pdfl[[c]][i,t]*(tau2[c]-w[i,t])) +
        (1/(cdfl[[c]][i,t]-cdfl[[c-1]][i,t]))*(pdfl[[c-1]][i,t]*(tau2[c-1]-w[i,t])^2 - pdfl[[c-1]][i,t] - pdfl[[c]][i,t]*(tau2[c]-w[i,t])^2 + pdfl[[c]][i,t])
      u2alpha[i,t,(1:k)] <- v2alpha[i,t]*x[i,t,]
      u2alpha[i,t,(k+q)] <- 2*(1/(cdfl[[c]][i,t]-cdfl[[c-1]][i,t])^3)*((pdfl[[c-1]][i,t]-pdfl[[c]][i,t])^2)*pdfl[[c]][i,t] -
        (1/(cdfl[[c]][i,t]-cdfl[[c-1]][i,t])^2)*(pdfl[[c-1]][i,t]*(tau2[c-1]-w[i,t]) - pdfl[[c]][i,t]*(tau2[c]-w[i,t]))*pdfl[[c]][i,t] -
        2*(1/(cdfl[[c]][i,t]-cdfl[[c-1]][i,t])^2)*(pdfl[[c-1]][i,t]-pdfl[[c]][i,t])*pdfl[[c]][i,t]*(tau2[c]-w[i,t]) +
        (1/(cdfl[[c]][i,t]-cdfl[[c-1]][i,t]))*pdfl[[c]][i,t]*((tau2[c]-w[i,t])^2) -
        (1/(cdfl[[c]][i,t]-cdfl[[c-1]][i,t]))*pdfl[[c]][i,t]
      
    }
    else if( (y[i,t]>1) & (y[i,t]<(p+1)) ){
      v[i,t] <- (1/(cdfl[[c]][i,t]-cdfl[[c-1]][i,t]))*(pdfl[[c-1]][i,t]-pdfl[[c]][i,t])
      u[i,t,(1:k)] <- v[i,t]*x[i,t,]
      u[i,t,(k+q-1)] <- - (1/(cdfl[[c]][i,t]-cdfl[[c-1]][i,t]))*(pdfl[[c-1]][i,t])
      u[i,t,(k+q)] <- (1/(cdfl[[c]][i,t]-cdfl[[c-1]][i,t]))*(pdfl[[c]][i,t]) 
      valpha[i,t] <- -(1/(cdfl[[c]][i,t]-cdfl[[c-1]][i,t])^2)*((pdfl[[c-1]][i,t]-pdfl[[c]][i,t])^2) + 
        (1/(cdfl[[c]][i,t]-cdfl[[c-1]][i,t]))*(pdfl[[c-1]][i,t]*(tau2[c-1]-w[i,t]) - pdfl[[c]][i,t]*(tau2[c]-w[i,t]))
      ualpha[i,t,(1:k)] <- valpha[i,t]*x[i,t,]
      ualpha[i,t,(k+q-1)] <- (1/(cdfl[[c]][i,t]-cdfl[[c-1]][i,t])^2)*(pdfl[[c-1]][i,t]-pdfl[[c]][i,t])*pdfl[[c-1]][i,t] -
        (1/(cdfl[[c]][i,t]-cdfl[[c-1]][i,t]))*pdfl[[c-1]][i,t]*(tau2[c-1]-w[i,t])
      ualpha[i,t,(k+q)] <- -(1/(cdfl[[c]][i,t]-cdfl[[c-1]][i,t])^2)*(pdfl[[c-1]][i,t]-pdfl[[c]][i,t])*pdfl[[c]][i,t] +
        (1/(cdfl[[c]][i,t]-cdfl[[c-1]][i,t]))*pdfl[[c]][i,t]*(tau2[c]-w[i,t]) 
      utheta[i,t,(1:k),(1:k)] <- valpha[i,t]*x[i,t,]%*%t(x[i,t,])
      utheta[i,t,(k+q-1),(k+q-1)] <- -(1/(cdfl[[c]][i,t]-cdfl[[c-1]][i,t])^2)*(pdfl[[c-1]][i,t]^2) + (1/(cdfl[[c]][i,t]-cdfl[[c-1]][i,t]))*(pdfl[[c-1]][i,t]*(tau2[c-1]-w[i,t]))
      utheta[i,t,(1:k),(k+q-1)] <- ualpha[i,t,(k+q-1)]*x[i,t,]
      utheta[i,t,(k+q-1),(1:k)] <- ualpha[i,t,(k+q-1)]*x[i,t,]
      utheta[i,t,(k+q),(k+q)] <- -(1/(cdfl[[c]][i,t]-cdfl[[c-1]][i,t])^2)*(pdfl[[c]][i,t]^2) - (1/(cdfl[[c]][i,t]-cdfl[[c-1]][i,t]))*(pdfl[[c]][i,t]*(tau2[c]-w[i,t]))       
      utheta[i,t,(1:k),(k+q)] <- ualpha[i,t,(k+q)]*x[i,t,]
      utheta[i,t,(k+q),(1:k)] <- ualpha[i,t,(k+q)]*x[i,t,]
      utheta[i,t,(k+q-1),(k+q)] <- (1/(cdfl[[c]][i,t]-cdfl[[c-1]][i,t])^2)*pdfl[[c-1]][i,t]*pdfl[[c]][i,t]
      utheta[i,t,(k+q),(k+q-1)] <- (1/(cdfl[[c]][i,t]-cdfl[[c-1]][i,t])^2)*pdfl[[c-1]][i,t]*pdfl[[c]][i,t]
      v2alpha[i,t] <- 2*(1/(cdfl[[c]][i,t]-cdfl[[c-1]][i,t])^3)*((pdfl[[c-1]][i,t]-pdfl[[c]][i,t])^3) - 
        3*(1/(cdfl[[c]][i,t]-cdfl[[c-1]][i,t])^2)*(pdfl[[c-1]][i,t]-pdfl[[c]][i,t])*(pdfl[[c-1]][i,t]*(tau2[c-1]-w[i,t]) - pdfl[[c]][i,t]*(tau2[c]-w[i,t])) +
        (1/(cdfl[[c]][i,t]-cdfl[[c-1]][i,t]))*(pdfl[[c-1]][i,t]*(tau2[c-1]-w[i,t])^2 - pdfl[[c-1]][i,t] - pdfl[[c]][i,t]*(tau2[c]-w[i,t])^2 + pdfl[[c]][i,t])
      u2alpha[i,t,(1:k)] <- v2alpha[i,t]*x[i,t,]
      u2alpha[i,t,(k+q-1)] <- - 2*(1/(cdfl[[c]][i,t]-cdfl[[c-1]][i,t])^3)*((pdfl[[c-1]][i,t]-pdfl[[c]][i,t])^2)*pdfl[[c-1]][i,t] +
        (1/(cdfl[[c]][i,t]-cdfl[[c-1]][i,t])^2)*(pdfl[[c-1]][i,t]*(tau2[c-1]-w[i,t]) - pdfl[[c]][i,t]*(tau2[c]-w[i,t]))*pdfl[[c-1]][i,t] +
        2*(1/(cdfl[[c]][i,t]-cdfl[[c-1]][i,t])^2)*(pdfl[[c-1]][i,t]-pdfl[[c]][i,t])*pdfl[[c-1]][i,t]*(tau2[c-1]-w[i,t]) -
        (1/(cdfl[[c]][i,t]-cdfl[[c-1]][i,t]))*pdfl[[c-1]][i,t]*((tau2[c-1]-w[i,t])^2) +
        (1/(cdfl[[c]][i,t]-cdfl[[c-1]][i,t]))*pdfl[[c-1]][i,t]
      u2alpha[i,t,(k+q)] <- 2*(1/(cdfl[[c]][i,t]-cdfl[[c-1]][i,t])^3)*((pdfl[[c-1]][i,t]-pdfl[[c]][i,t])^2)*pdfl[[c]][i,t] -
        (1/(cdfl[[c]][i,t]-cdfl[[c-1]][i,t])^2)*(pdfl[[c-1]][i,t]*(tau2[c-1]-w[i,t]) - pdfl[[c]][i,t]*(tau2[c]-w[i,t]))*pdfl[[c]][i,t] -
        2*(1/(cdfl[[c]][i,t]-cdfl[[c-1]][i,t])^2)*(pdfl[[c-1]][i,t]-pdfl[[c]][i,t])*pdfl[[c]][i,t]*(tau2[c]-w[i,t]) +
        (1/(cdfl[[c]][i,t]-cdfl[[c-1]][i,t]))*pdfl[[c]][i,t]*((tau2[c]-w[i,t])^2) -
        (1/(cdfl[[c]][i,t]-cdfl[[c-1]][i,t]))*pdfl[[c]][i,t]
      
    }
    else{
      v[i,t] <- (1/(cdfl[[c]][i,t]-cdfl[[c-1]][i,t]))*(pdfl[[c-1]][i,t]-pdfl[[c]][i,t])
      u[i,t,(1:k)] <- v[i,t]*x[i,t,]
      u[i,t,(k+q-1)] <- - (1/(cdfl[[c]][i,t]-cdfl[[c-1]][i,t]))*(pdfl[[c-1]][i,t])
      valpha[i,t] <- -(1/(cdfl[[c]][i,t]-cdfl[[c-1]][i,t])^2)*((pdfl[[c-1]][i,t]-pdfl[[c]][i,t])^2) + 
        (1/(cdfl[[c]][i,t]-cdfl[[c-1]][i,t]))*(pdfl[[c-1]][i,t]*(tau2[c-1]-w[i,t]))
      ualpha[i,t,(1:k)] <- valpha[i,t]*x[i,t,]
      ualpha[i,t,(k+q-1)] <- (1/(cdfl[[c]][i,t]-cdfl[[c-1]][i,t])^2)*(pdfl[[c-1]][i,t]-pdfl[[c]][i,t])*pdfl[[c-1]][i,t] -
        (1/(cdfl[[c]][i,t]-cdfl[[c-1]][i,t]))*pdfl[[c-1]][i,t]*(tau2[c-1]-w[i,t])
      utheta[i,t,(1:k),(1:k)] <- valpha[i,t]*x[i,t,]%*%t(x[i,t,])
      utheta[i,t,(k+q-1),(k+q-1)] <- -(1/(cdfl[[c]][i,t]-cdfl[[c-1]][i,t])^2)*(pdfl[[c-1]][i,t]^2) + (1/(cdfl[[c]][i,t]-cdfl[[c-1]][i,t]))*(pdfl[[c-1]][i,t]*(tau2[c-1]-w[i,t]))
      utheta[i,t,(1:k),(k+q-1)] <- ualpha[i,t,(k+q-1)]*x[i,t,]
      utheta[i,t,(k+q-1),(1:k)] <- ualpha[i,t,(k+q-1)]*x[i,t,]
      v2alpha[i,t] <- 2*(1/(cdfl[[c]][i,t]-cdfl[[c-1]][i,t])^3)*((pdfl[[c-1]][i,t]-pdfl[[c]][i,t])^3) - 
        3*(1/(cdfl[[c]][i,t]-cdfl[[c-1]][i,t])^2)*(pdfl[[c-1]][i,t]-pdfl[[c]][i,t])*(pdfl[[c-1]][i,t]*(tau2[c-1]-w[i,t])) +
        (1/(cdfl[[c]][i,t]-cdfl[[c-1]][i,t]))*(pdfl[[c-1]][i,t]*(tau2[c-1]-w[i,t])^2 - pdfl[[c-1]][i,t])
      u2alpha[i,t,(1:k)] <- v2alpha[i,t]*x[i,t,]
      u2alpha[i,t,(k+q-1)] <- - 2*(1/(cdfl[[c]][i,t]-cdfl[[c-1]][i,t])^3)*((pdfl[[c-1]][i,t]-pdfl[[c]][i,t])^2)*pdfl[[c-1]][i,t] +
        (1/(cdfl[[c]][i,t]-cdfl[[c-1]][i,t])^2)*(pdfl[[c-1]][i,t]*(tau2[c-1]-w[i,t]))*pdfl[[c-1]][i,t] +
        2*(1/(cdfl[[c]][i,t]-cdfl[[c-1]][i,t])^2)*(pdfl[[c-1]][i,t]-pdfl[[c]][i,t])*pdfl[[c-1]][i,t]*(tau2[c-1]-w[i,t]) -
        (1/(cdfl[[c]][i,t]-cdfl[[c-1]][i,t]))*pdfl[[c-1]][i,t]*((tau2[c-1]-w[i,t])^2) +
        (1/(cdfl[[c]][i,t]-cdfl[[c-1]][i,t]))*pdfl[[c-1]][i,t]
    }

    sumuv[i,] <- sumuv[i,] + ((v[i,t]*u[i,t,]))/T
    sumvv[i] <- sumvv[i] + (v[i,t]^2)/T
    sumvalpha[i] <- sumvalpha[i] + (valpha[i,t])/T
    sumualpha[i,] <- sumualpha[i,] + (ualpha[i,t,])/T

    sumvualpha[i,] <- sumvualpha[i,] + (v[i,t]*ualpha[i,t,])/T
    sumvvalpha[i] <- sumvvalpha[i] + (v[i,t]*valpha[i,t])/T
    sumv2alpha[i] <- sumv2alpha[i] + (v2alpha[i,t])/T
    sumu2alpha[i,] <- sumu2alpha[i,] + (u2alpha[i,t,])/T
    sumutheta[i,,] <- sumutheta[i,,] + (utheta[i,t,,])/T
  }
  
  part1b[i] <- -sumvv[i,]/sumvalpha[i]
  part2b[i,] <- -(1/-sumvv[i])*(sumvualpha[i,] - sumvvalpha[i]*(sumualpha[i,]/sumvalpha[i]))
  part3b[i,] <- -(1/(2*sumvalpha[i]))*(sumu2alpha[i,] - sumv2alpha[i]*(sumualpha[i,]/sumvalpha[i]))
  
  I[i,,] <- -(sumutheta[i,,] - (1/ sumvalpha[i])*sumualpha[i,]%*%t(sumualpha[i,]))
  b[i,] <- part1b[i]*(part2b[i,]+part3b[i,])

}

B <- solve(colMeans(I))%*%colMeans(b)
bias <- (1/T)*B

gammaold = gamma
gamma = gamma - bias