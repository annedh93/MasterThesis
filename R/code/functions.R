#WOOLRIDGE LL ###########################################################################################################
LL <- function(gamma){
  beta <- gamma[1:totreg]
  tau <- gamma[(totreg+1):(totreg+cat)]
  cdfl <- list()
  pdfl <- list()
  tau2 <- c(-Inf,tau,Inf) # set one to 0 for identification
  w <- eval(parse(text=paste(paste0("beta","[",1:totreg,"]*","X","[,,",1:totreg,"]",collapse="+"))))
  
  for(l in 1:length(tau2)){
    cdfl[[l]] <- pnorm(tau2[l] - w)
  }
  
  logl <- matrix(0,N,T)
  for(j in 1:(cat+1)){
    ll <- (y==(j-1))*log(cdfl[[j+1]]-cdfl[[j]])
    ll[(y!=(j-1))] <- 0
    logl <- logl + ll
  }
  return(sum(logl,na.rm=TRUE))
}

#ODERED LOGIT LL ########################################################################################################
LL2 <- function(beta){
  dich <- (y >= 2)
  Pi <- rep(0,N)
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
      
      denom <- 0
      for(l in 1:comb){
        denom <- denom + exp(t(Bi[index,l])%*%x[i,index,]%*%beta)
      }
      Pi[i] <- exp(t(dich[i,index])%*%x[i,index,]%*%beta)/denom  
    }
  }
  return(-sum(log(Pi),na.rm=TRUE))
}

#create dummy matrix ###################################################################################################
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

#create dummy matrix with all dummies ###################################################################################
as.dummy2 <- function(x,name){
  tab <- table(x)
  N <- length(x)
  G <- matrix(NA,N,length(tab))
  for(i in 1:length(tab)){
    G[which(x==names(tab[i])),i] <- 1
    G[which(x!=names(tab[i])),i] <- 0
  }
  colnames(G) <- paste0(name,"d",1:length(tab))
  return(G)
}

#ROUTING ###########################################################################################################
routing <- function(samp,selection=FALSE){
  for(i in 1:(nrow(samp)/5)){
    coup <- FALSE
    wavesp <- NULL
    
    if(selection==TRUE){
    obs <- which(!is.na(indexinc[datsamp[(i*5-4):(i*5)]]))
    indsamp <- c((i*5-5+obs[1]):(i*5-5+obs[length(obs)]))
    inddat <- datsamp[indsamp]
    }else{
      obs <- which(!is.na(samp[(i*5-4):(i*5),4]))
      indsamp <- c((i*5-5+obs[1]):(i*5-5+obs[length(obs)]))
      inddat <- indsamp
    }
    coupleid <- data[inddat,"numbercid"]
    if(sum(coupleid,na.rm=TRUE)>length(inddat)){
      coup <- TRUE
      indcouple <- which(data[,"numbercid"]==coupleid[coupleid>1][1])
      min <- match(inddat,indcouple)
      min <- min[complete.cases(min)]
      indpartner <- indcouple[-min]
      wavesp <- data[indpartner,"wave"]
    }
    waves <- data[inddat,"wave"]
    
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

#MAKE MATRIX WITH FIVE WAVES FOR ALL OBSERVATIONS#############################################################################
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

######## DAS EN VAN SOEST ORDERED LOGIT ######################################################################################################
orderedlogit <- function(y,x,c,regressors){
N <- nrow(y)
T <- ncol(y)
p <- dim(x)[3]
X <- x

# dichotomized ordered logit model
betak <- matrix(0,p,c)
Hbetak <- matrix(0,(c-1)*p,(c-1)*p)
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
          if(length(index)>1){
          sum <- sum + as.vector(exp(t(Bi[index,m])%*%X[i,index,]%*%beta)/denom)*(t(X[i,index,])%*%Bi[index,m] - 
                 num2/as.vector(denom))%*%t(t(X[i,index,])%*%Bi[index,m] - num2/as.vector(denom))
          }else{
            sum <- sum + as.vector(exp(t(Bi[index,m])%*%X[i,index,]%*%beta)/denom)*(as.vector(t(X[i,index,])*Bi[index,m]) - 
                 num2/as.vector(denom))%*%t(as.vector(t(X[i,index,])*Bi[index,m]) - num2/as.vector(denom))
          }
        }
        if(length(index)>1){
          gbeta[,i,k] <- - t(X[i,index,])%*%(as.matrix(dich[i,index]) - num1/as.vector(denom))
        }else{
          gbeta[,i,k] <- - as.vector(t(X[i,index,]))%*%(as.matrix(dich[i,index]) - num1/as.vector(denom))
        }
        Hbeta <- Hbeta + sum
      }
    }
    Hbetak[((k-1)*p-p+1):((k-1)*p),((k-1)*p-p+1):((k-1)*p)] <- Hbeta
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
      mom2[(a*p-p+1):(a*p),(b*p-p+1):(b*p)] <- mom2[(a*p-p+1):(a*p),(b*p-p+1):(b*p)] + gbeta[,i,(a+1)]%*%t(gbeta[,i,(b+1)])
    }
    mom1[(a*p-p+1):(a*p),(a*p-p+1):(a*p)] <- mom1[(a*p-p+1):(a*p),(a*p-p+1):(a*p)] + gbeta[,i,(a+1)]%*%t(gbeta[,i,(a+1)])
  }
}

omega <- solve(mom1)%*%(mom2)%*%solve(mom1)
omega2 <- solve(-Hbetak)%*%(mom2)%*%solve(-Hbetak) # no 1/N because this is variance of beta not sqrt(n)(beta-mean(beta))
I <- matrix(rep(diag(1,p),(c-1)),p*(c-1),p,byrow=TRUE)
betaf <- as.vector(solve(t(I)%*%solve(omega)%*%I)%*%t(I)%*%solve(omega)%*%c(betak[,2:c]))
betaf2 <- as.vector(solve(t(I)%*%solve(omega2)%*%I)%*%t(I)%*%solve(omega2)%*%c(betak[,2:c]))
names(betaf2) <- regressors
var <- solve(t(I)%*%solve(omega)%*%I)
var2 <- solve(t(I)%*%solve(omega2)%*%I)
sd <- sqrt(diag(var2))
z <- betaf2/sd
pvalue <- pnorm(abs(z),lower.tail=FALSE)*2

return(cbind(betaf2,sd,z,pvalue))
}

######## BAETSCHMANN ######################################################################################################
#x <- array(rnorm(100*5*3),dim = c(100,5,3))
#beta <- c(5,-0.5,2)
#ystar <- x[,,1]*beta[1] + x[,,2]*beta[2] + x[,,3]*beta[3] + matrix(rlogis(100*5),100,5)
#y <- (ystar>-1) + (ystar>0) + (ystar>1) + (ystar>2) + 1
orderedlogitB <- function(y,x,cat,regressors){
  N <- nrow(y)
  T <- ncol(y)
  p <- dim(x)[3]
  X <- x
  
  # dichotomized ordered logit model
  gbetak <- array(0,dim=c(p,N,cat))

  beta <- rep(0,p)
  betaold <- NULL
  
  dif <- rep(100,p)
  eps <- rep(10^-3,p)
  
  while(sum(dif>eps)>0){
    Hbeta <- matrix(0,p,p)
    gbeta <- matrix(0,p,1)
    gbetak <- array(0,dim=c(p,N,cat))
    LL <- 0
    for(k in 2:cat){
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
            if(length(index)>1){
              sum <- sum + as.vector(exp(t(Bi[index,m])%*%X[i,index,]%*%beta)/denom)*(t(X[i,index,])%*%Bi[index,m] - 
                     num2/as.vector(denom))%*%t(t(X[i,index,])%*%Bi[index,m] - num2/as.vector(denom))
            }else{
              sum <- sum + as.vector(exp(t(Bi[index,m])%*%X[i,index,]%*%beta)/denom)*(as.vector(t(X[i,index,])*Bi[index,m]) - 
                     num2/as.vector(denom))%*%t(as.vector(t(X[i,index,])*Bi[index,m]) - num2/as.vector(denom))
            }
          }
          if(length(index)>1){
            gbeta <- gbeta - t(X[i,index,])%*%(as.matrix(dich[i,index]) - num1/as.vector(denom))
            gbetak[,i,k] <- - t(X[i,index,])%*%(as.matrix(dich[i,index]) - num1/as.vector(denom))
          }else{
            gbeta <- gbeta - as.vector(t(X[i,index,]))%*%(as.matrix(dich[i,index]) - num1/as.vector(denom))
            gbetak[,i,k] <- - as.vector(t(X[i,index,]))%*%(as.matrix(dich[i,index]) - num1/as.vector(denom))
          }
          Hbeta <- Hbeta + sum
          LL <- LL + log(exp(t(as.matrix(dich[i,index]))%*%X[i,index,]%*%beta)/denom)
        }
      }
    }
    betaold <- cbind(betaold,beta)
    deltabeta <- -solve(Hbeta)%*%gbeta
    beta <- beta + deltabeta
    dif <- abs(betaold[,ncol(betaold)]-beta)
    show(sum(dif>eps))
  }
  
  mom <- matrix(0,p,p)
  for(i in 1:N){
    for(a in 1:(cat-1)){
      for(b in 1:(cat-1)){
        mom <- mom + gbetak[,i,(a+1)]%*%t(gbetak[,i,(b+1)])
      }
    }
  }
  
  betaf <- beta
  names(betaf) <- regressors
  var <- solve(Hbeta)%*%(mom)%*%solve(Hbeta) # no 1/N because this is variance of beta not sqrt(n)(beta-mean(beta))
  sd <- sqrt(diag(var))
  z <- betaf/sd
  pvalue <- pnorm(abs(z),lower.tail=FALSE)*2
  output <- cbind(betaf,sd,z,pvalue)
  rownames(output) <- regressors
  colnames(output) <- c("beta","sd","z","pvalue")
  
  return(output)
}
