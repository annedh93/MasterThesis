steplength <- function(ro,c,stepbeg,par,x,delta,der){
  step <- stepbeg
  while(LL(par + step*delta,x,y) < (LL(par,x,y) + c*step*t(der)%*%delta)){
    step <- ro*step
    show(step)
  }
  return(step)
}

LL <- function(par,x,y){
  beta <- par[1:k]
  tau <- par[(k+1):(k+p)]
  alpha <- par[(k+p+1):length(par)]
  w <- eval(parse(text=paste(paste0("beta","[",1:k,"]*","x","[,,",1:k,"]",collapse="+")))) + 
    matrix(alpha,nrow=N,ncol=T,byrow=FALSE)
  
  for(l in 1:length(tau2)){
    cdfl[[l]] <- pnorm(tau2[l] - w)
  }
  
  logl <- matrix(0,N,T)
  for(j in 1:4){
    ll <- (y==(j-1))*log(cdfl[[j+1]]-cdfl[[j]])
    ll[(y!=(j-1))] <- 0
    logl <- logl + ll
  }
  return(sum(logl,na.rm=TRUE))
}

LLna <- function(par,x,y,N,T,k){
  ind <- which(rowSums(is.na(y))<(T-1))
  ydat <- y[ind,] - rowMeans(y[ind,],na.rm=TRUE)
  xdat <- array(NA,dim=c(length(ind),T,k))
  for(p in 1:k){
    xdat[,,p] <- x[ind,,p] - rowMeans(x[ind,,p],na.rm=TRUE)
  }

  beta <- par
  logl <- log( (1/sqrt(2*pi))*exp(-(1/2)*((ydat-eval(parse(text=paste(paste0("beta","[",1:k,"]*","xdat","[,,",1:k,"]",collapse="+")))))^2)))

  return(-sum(logl,na.rm=TRUE))
}

NR <- function(par,x,y,N,T,k){
  ind <- which(rowSums(is.na(y))<(T-1))
  ydat <- y[ind,] - rowMeans(y[ind,],na.rm=TRUE)
  xdat <- array(NA,dim=c(length(ind),T,k))
  for(p in 1:k){
    xdat[,,p] <- x[ind,,p] - rowMeans(x[ind,,p],na.rm=TRUE)
  }
 
  N <- nrow(ydat)

  beta <- par
  betaold <- rep(10^6,k)
  while(sum(abs(beta-betaold))>(10^-3)){
    grad <- rep(0,k)
    Hes <- matrix(0,k,k)
    for(i in 1:N){
      for(t in 1:T){
        plus <- ydat[i,t]*xdat[i,t,] - xdat[i,t,]%*%t(xdat[i,t,])%*%beta 
        if(!is.na(sum(plus))){grad <- grad + plus}
        plus2 <- - xdat[i,t,]%*%t(xdat[i,t,])
        if(!is.na(sum(plus2))){Hes <- Hes + plus2}
      }
    }
    betaold <- beta
    beta <- beta - solve(Hes)%*%grad
  }
  return(beta)
}

as.dummy <- function(x,omit,name){
  tab <- table(x)
  N <- length(x)
  G <- matrix(0,N,(length(tab)-1))
  k <- 1
  for(i in 1:length(tab)){
    if(i != omit){
      G[which(x==names(tab[i])),k] <- 1
      k <- k + 1
    }
  }
  colnames(G) <- paste0(name,"d",1:(length(tab)-1))
  return(G)
}

as.dummy2 <- function(x,name){
  tab <- table(x)
  N <- length(x)
  G <- matrix(0,N,length(tab))
  for(i in 1:length(tab)){
    G[which(x==names(tab[i])),i] <- 1
  }
  colnames(G) <- paste0(name,"d",1:length(tab))
  return(G)
}

routing <- function(samp){
# select random sample from data
numberid <- as.numeric(as.factor(samp[,"id"]))
samp <- cbind(numberid,samp)
samp <- samp[order(samp[,"numberid"]),]
uniq <- unique(numberid)

for(i in 1:length(uniq)){
  coup <- FALSE
  wavesp <- NULL
  
  obs <- samp[which(samp[,"numberid"]%in%uniq[i]),"mergeid"][1]
  indsamp <- which(samp[,"mergeid"]%in%obs)
  inddat <- which(data[,"mergeid"]%in%obs)
  waves <- samp[indsamp,"wave"]
  
  if(sum(data[inddat,"coupleid"]%in%"")!=length(inddat)){
    coup <- TRUE
    couple <- data[inddat[which(!(data[inddat,"coupleid"]%in%""))][1],"coupleid"]
    indcouple <- which(data[,"coupleid"]%in%couple)
    min <- match(inddat,indcouple)
    min <- min[complete.cases(min)]
    indpartner <- indcouple[-min]
    wavesp <- data[indpartner,"wave"]
  }
  
  for(t in 1:length(waves)){
    
    # marital status
    if(data[inddat[t],"mstatc"]==4 & !is.na(data[inddat[t],"mstatc"]) & is.na(samp[indsamp[t],"mstat"]) & (t>1)){
      samp[indsamp[t],"mstat"] <- samp[(indsamp[t-1]),"mstat"]
      data[inddat[t],"mstat"] <- data[(inddat[t-1]),"mstat"]
    }
    if(is.na(samp[indsamp[t],"mstat"]) & waves[t]%in%wavesp & coup){
      for(tp in 1:length(wavesp)){
        if(data[indpartner[tp],"mstatc"]==4 & !is.na(data[indpartner[tp],"mstatc"]) & is.na(data[indpartner[tp],"mstat"]) & (tp>1)){
          data[indpartner[tp],"mstat"] <- data[(indpartner[tp-1]),"mstat"]
        } 
      }
      tp <- which(wavesp%in%waves[t])
      samp[indsamp[t],"mstat"] <- data[indpartner[tp],"mstat"] 
    }
    
    # children
    if(is.na(samp[indsamp[t],"child"]) & waves[t]%in%wavesp & coup){
      tp <- which(wavesp%in%waves[t])
      samp[indsamp[t],"child"] <- data[indpartner[tp],"child"]
    }
    if(is.na(samp[indsamp[t],"child"])){
      samp[indsamp[t],"child"] <- data[inddat[which(!is.na(data[inddat,"child"]))][1],"child"]
    }
    
    # education
    if(is.na(samp[indsamp[t],"isced"])){
      if(t > 1){
        samp[indsamp[t],"isced"] <- data[(inddat[t]-1),"isced"]
      }
      if(is.na(samp[indsamp[t],"isced"]) & sum(!is.na(samp[indsamp,"isced"]))>0){
        samp[indsamp[t],"isced"] <- samp[indsamp[which(!is.na(samp[indsamp,"isced"]))][1],"isced"]
      }
    }
    
  }
  # marital status
  if(sum(is.na(samp[indsamp,"mstat"]))>0 & sum(samp[indsamp[which(is.na(samp[indsamp,"mstat"]))],"mstatc"]%in%4)==sum(is.na(samp[indsamp,"mstat"]))){
    samp[indsamp[which(is.na(samp[indsamp,"mstat"]))],"mstat"] <- data[inddat[which(!is.na(samp[inddat,"mstat"]))][1],"mstat"]
  }
  
  # children
  if(sum(is.na(samp[indsamp,"child"]))>0){
    samp[indsamp[which(is.na(samp[indsamp,"child"]))],"child"] <- data[inddat[which(!is.na(data[inddat,"child"]))][1],"child"]
  }
  
  show(i)
}
return(samp)
}

fivewave <- function(panel){
# first transform data to contain five waves again
w <- c(1,2,4,5,6)
data <- panel[order(panel[,"numberid"]),]
for(i in 1:length(unique(data[,"numberid"]))){
  inddat <- which(data[,"numberid"]==unique(data[,"numberid"])[i])
  waves <- data[inddat,"wave"]
  for(k in 1:5){
    if(!(w[k]%in%waves)){
      if(inddat[1]==1 & k==1){
        datanew <- as.data.frame(matrix(NA,(nrow(data)+1),ncol(data)))
        colnames(datanew) <- colnames(data)
        datanew[2:nrow(datanew),] <- data[1:nrow(data),] 
        datanew[1,c(1:2)] <- data[inddat[1],c(1:2)]
        datanew[1,3] <- w[k]
        data <- datanew
      }
      else{
        indw <- inddat[1] - 2 + k
        datanew <- as.data.frame(matrix(NA,(nrow(data)+1),ncol(data)))
        colnames(datanew) <- colnames(data)
        datanew[1:indw,] <- data[1:indw,]
        datanew[(indw+2):nrow(datanew),] <- data[(indw+1):nrow(data),] 
        datanew[(indw+1),c(1:2)] <- data[inddat[1],c(1:2)]
        datanew[(indw+1),3] <- w[k]
        data <- datanew
      }
    } 
  }
  show(i)
}
return(data)
}