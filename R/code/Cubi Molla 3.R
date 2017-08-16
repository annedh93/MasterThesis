############################################################### construct data set #############################################################

load("C:\\Users\\anned\\OneDrive\\Documenten\\master thesis\\R\\data\\samp.Rdata")

# find sample that experiences change in ADL
incidence <- samp[,"IADL"] # "chronic" or "IADL"
indexinc <- (incidence==1) # > 0 for "IADL, == 3 for "chronic
indexinc2 <- matrix(FALSE,nrow(samp),1)
dur <- matrix(0,nrow(samp),1)
sphuslag <- matrix(NA,nrow(samp),1)
first <- matrix(FALSE,nrow(samp),1)
firstsphus <- matrix(NA,nrow(samp),1)
for(i in 1:(nrow(samp)/5)){
  obs <- which(!is.na(indexinc[(i*5-4):(i*5)]))
  if(sum(diff(obs))==(length(obs)-1)){ # check whether observations lie in consecutive waves
    firstobs <- (i*5-5+obs[1])
    lastobs <- (i*5-5+obs[length(obs)])
    incobs <- which(indexinc[firstobs:lastobs]==TRUE)
    waves <- samp[firstobs:lastobs,"wave"][incobs]
    if(length(incobs)!=0 & length(incobs)==length(which(samp[firstobs:lastobs,"IADL"]>0))){ # length(incobs)==length(samp[firstobs:lastobs,"IADL"]>0)
      if(sum(diff(incobs))==(length(incobs)-1) & indexinc[firstobs]==FALSE & indexinc[lastobs]==TRUE & sum(samp[firstobs:lastobs,"chronic"][incobs],na.rm=TRUE)==(length(incobs)*3)){ 
        #length(unique(samp[firstobs:lastobs,"IADL"][incobs]))==1 & length(unique(indexinc[incobs]))==1
        indexinc2[(i*5-4):(i*5)] <- TRUE
        sphuslag[(i*5-3):(i*5)] <- samp[(i*5-4):(i*5-1),"sphus"]
        firstsphus[(i*5-4):(i*5)] <- samp[firstobs,"sphus"]
        first[firstobs] <- TRUE
        for(k in 1:length(incobs)){
          if(k==1){
            dur[(firstobs-1+incobs[1]),1] <- abs(samp[firstobs:lastobs,"age"][incobs[1]]-samp[firstobs:lastobs,"age"][(incobs[1]-1)])/2
          }else{
            dur[(firstobs-1+incobs[k]),1] <- abs(samp[firstobs:lastobs,"age"][incobs[k]]-samp[firstobs:lastobs,"age"][(incobs[k]-1)]) + dur[(firstobs-1+incobs[k-1]),1]
          }
        }
      }
    }
  }
}

samp <- cbind(samp[,1:18],firstsphus,sphuslag,first,samp[,19:23],indexinc,dur,samp[,24:ncol(samp)])
samp <- samp[indexinc2,]
indexsamp <- which(!is.na(samp[,"sphus"]) & !is.na(samp[,"sphuslag"]) & !is.na(samp[,"firstsphus"]))
rsamp <- samp[-which(is.na(samp[,"sphus"]) | is.na(samp[,"sphuslag"]) | is.na(samp[,"firstsphus"])),]
IADLi <- (rsamp[,"IADL"]==1)
rsamp <- cbind(rsamp,IADLi)
rsamp[,"IADLi"] <- as.numeric(rsamp[,"IADLi"])

# form imputations
#install.packages("Amelia")
#library(Amelia)
c <- which(colnames(rsamp)%in%c("numberid","id","wave","country","firstwave","mstat","age","isced",
                                "child","lifesat","sphus","firstsphus","sphuslag","dur","IADL","employment"))
M <- round(100*length(unique(which(is.na(rsamp),arr.ind=TRUE)[,1]))/nrow(rsamp))
imputations <- amelia(rsamp[,c],m=M,ts="wave",cs="id",
                      idvars=c("numberid","country","firstwave"),
                      noms=c("mstat","employment"),
                      ords=c("sphus","firstsphus","sphuslag","lifesat"))
implist <- imputations$imputations

############################################################################################################################
resultslist <- list()
for(m in 1:M){
  
  imp <- implist[[m]]
  # age (16-20,21-30,31-40,41-50,51-60,61-70,71-80,80+)
  imp[(imp[,"age"]<=50),"age"] <- 1
  imp[(imp[,"age"]>50 & imp[,"age"]<=60),"age"] <- 2
  imp[(imp[,"age"]>60 & imp[,"age"]<=70),"age"] <- 3
  imp[(imp[,"age"]>70 & imp[,"age"]<=80),"age"] <- 4
  imp[(imp[,"age"]>80),"age"] <- 5
  # marital status (married, single, widowed, divorced, separated)
  imp[(imp[,"mstat"]==3 | imp[,"mstat"]==4),"mstat"] <- 1 # group Married and Registered partnership together
  imp[(imp[,"mstat"]>=5),"mstat"] <- 2
  #imp[(imp[,"mstat"]>=5),"mstat"] <- imp[(imp[,"mstat"]>=5),"mstat"]-3
  # labour-force status (employed, unemployed, retired, inactive)
  imp[(imp[,"employment"]>=6),"employment"] <- imp[(imp[,"employment"]>=6),"employment"] <- 6
  imp[,"employment"] <- imp[,"employment"]-2
  # education (high school, less than high school, more than high school)
  imp[(imp[,"isced"])<11,"isced"] <- 1
  imp[(imp[,"isced"]>=11 & imp[,"isced"]<=12),"isced"] <- 2
  imp[(imp[,"isced"]>12),"isced"] <- 3
  # sphus from poor to excellent
  imp[,"sphus"] <- abs(imp[,"sphus"]-7)
  imp[,"sphuslag"] <- abs(imp[,"sphuslag"]-7)
  imp[,"firstsphus"] <- abs(imp[,"firstsphus"]-7)
  # life satisfaction into 5 categories
  imp[(imp[,"lifesat"]<=3),"lifesat"] <- 1
  imp[(imp[,"lifesat"]<=4 & imp[,"lifesat"]>3),"lifesat"] <- 2
  imp[(imp[,"lifesat"]<=5 & imp[,"lifesat"]>4),"lifesat"] <- 3
  imp[(imp[,"lifesat"]<=6 & imp[,"lifesat"]>5),"lifesat"] <- 4
  imp[(imp[,"lifesat"]<=7 & imp[,"lifesat"]>6),"lifesat"] <- 5
  imp[(imp[,"lifesat"]<=8 & imp[,"lifesat"]>7),"lifesat"] <- 6
  imp[(imp[,"lifesat"]<=9 & imp[,"lifesat"]>8),"lifesat"] <- 7
  imp[(imp[,"lifesat"]<=10 & imp[,"lifesat"]>9),"lifesat"] <- 8
  
  # dynamic ordered probit
  panel <- as.matrix(cbind(imp[,c("numberid","id","wave","sphus","lifesat","IADL","dur","child")],
                           as.dummy(imp[,"sphuslag"],1,"sphuslag"),as.dummy(imp[,"firstsphus"],1,"firstsphus"),
                           as.dummy(imp[,"mstat"],1,"mstat"),as.dummy(imp[,"employment"],2,"employment"),
                           as.dummy(imp[,"age"],1,"age"),as.dummy(imp[,"isced"],2,"isced")))
  panel2 <- matrix(NA,nrow(samp),ncol(panel))
  panel2[indexsamp,] <- panel
  colnames(panel2) <- colnames(panel)
  panel2 <- as.data.frame(panel2)
  panel2[,c("numberid","id","wave")] <- samp[,c("numberid","id","wave")]
  panel2 <- as.matrix(panel2)
  regressors <- c(paste0("age","d",1:4),paste0("mstat","d",1),paste0("employment","d",1:3),paste0("isced","d",1:2),"child","IADL","dur") 
  #regressors <- c(paste0("age","d",1:4))
  lags <- c(paste0("sphuslag","d",1:3),paste0("firstsphus","d",1:3))
  N <- nrow(panel2)/5 # datch = samp
  T <- 5
  k <- length(regressors) # number of regressors
  l <- length(lags)/2
  cat <- 4 # number of categories - 1
  y <- matrix(as.numeric(panel2[,"sphus"]),N,T,byrow=TRUE) # 0 = Poor, 1 = Fair, 2 = Good, 3 = Excellent
  x <- array(NA,dim=c(N,T,(k+2*l)))
  #x <- array(NA,dim=c(N,T,(k)))
  #for(p in 1:k){
  #   x[,,p] <- matrix(as.numeric(panel2[,regressors[p]]),N,T,byrow=TRUE) 
  #}
  for(p in 1:l){
    x[,,p] <- matrix(as.numeric(panel2[,lags[p]]),N,T,byrow=TRUE) # lagged variable 1:3
  }
  for(p in 1:k){
    q <- l+p
    x[,,q] <- matrix(as.numeric(panel2[,regressors[p]]),N,T,byrow=TRUE) # all regressors 4:17
  }
  for(p in 1:l){
    q <- l+k+p
    x[,,q] <- matrix(as.numeric(panel2[,lags[(l+p)]]),N,T,byrow=TRUE) # first SAH observed in the sample 18:20
  }
  
  #means <- eval(parse(text = paste0("array(c(",paste0("matrix(rowMeans(x[,,",1:k,"],na.rm=TRUE),N,T,byrow=FALSE)",collapse=","),"),dim=c(N,T,k))")))
  #X <- array(c(x,means),dim=c(N,T,(2*k)))
  
  means <- eval(parse(text = paste0("array(c(",paste0("matrix(rowMeans(x[,,",4:(k+3),"],na.rm=TRUE),N,T,byrow=FALSE)",collapse=","),"),dim=c(N,T,k))")))
  X <- array(c(x,means),dim=c(N,T,2*(k+l)))
  totreg <- 2*(k+l)
  
  # optimize log likelihood
  gamma <- c(rep(0,totreg),1,2,3,4)
  par <- optim(gamma,LL,method="BFGS",hessian=TRUE)
  parameters <- par$par
  #names(parameters) <- regressors
  names(parameters) <- c(paste0("sphuslag","d",1:3),regressors,paste0("firstsphus","d",1:3),paste0(regressors,"mean"),"tau1","tau2","tau3","tau4") 
  Hessian <- par$hessian
  cov <- solve(Hessian) 
  sd <- sqrt(diag(cov))
  z <- parameters/sd
  pvalue <- pnorm(abs(z),lower.tail=FALSE)*2
  resultslist[[m]] <- as.matrix(cbind(parameters,sd,z,pvalue))

  print(paste("m is",m))
}

###########################################################################################################################################

# calculate pooled estimate
output <- as.data.frame(matrix(NA,k,5))
rownames(output) <- regressors
colnames(output) <- c("Estimate","Std. Error","z value","Pr(>|t|)","star")

# calculate pooled regression coefficients
coef <- matrix(NA,k,M)
U <- matrix(NA,k,M)
for(m in 1:M){
  result <- resultslist[[m]]
  coef[,m] <- result[,1]
  U[,m] <- (result[,2])^2
}

mipool <- mi.meld(coef,sqrt(U),byrow=FALSE)
output[,"Estimate"] <- t(mipool$q.mi)
output[,"Std. Error"] <- t(mipool$se.mi)
output[,"z value"] <- output[,"Estimate"]/output[,"Std. Error"]
output[,"Pr(>|t|)"] <- pnorm(abs(output[,"z value"]),lower.tail=FALSE)*2
star <- rep("",k)
star[output[,"Pr(>|t|)"] < 0.05] <- "*"
star[output[,"Pr(>|t|)"] < 0.01] <- "**"
star[output[,"Pr(>|t|)"] < 0.001] <- "***"
output[,"star"] <- star

# partial effects
beta <- parameters[1:(2*(k+l))]
tau <- parameters[(2*(k+l)+1):(2*(k+l)+cat)]
cdfl <- list()
pdfl <- list()
tau2 <- c(-Inf,tau,Inf) # set one to 0 for identification
w <- eval(parse(text=paste(paste0("beta","[",1:(2*(k+l)),"]*","x","[,,",1:(2*(k+l)),"]",collapse="+"))))

for(l in 1:length(tau2)){
  pdfl[[l]] <- dnorm(tau2[l] - w)    
}

# scoring Excellent
part <- matrix(NA,8,4)
c <- which(names(parameters)%in%c("lagsphusd1","lagsphusd2","lagsphusd3","LSI","dur","firstsphusd1","firstsphusd2","firstsphusd3"))
colnames(part) <- c("Excellent","Good","Fair","Poor")
rownames(part) <- c("lagsphusd1","lagsphusd2","lagsphusd3","LSI","dur","firstsphusd1","firstsphusd2","firstsphusd3")
for(j in 1:4){
  for(k in 1:8){
    part[k,j] <- -mean((pdfl[[j+1]]-pdfl[[j]])*beta[c[k]],na.rm=TRUE)
  }
}

