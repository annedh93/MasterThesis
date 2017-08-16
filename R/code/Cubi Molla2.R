############################################################### construct data set #############################################################

load("C:\\Users\\anned\\OneDrive\\Documenten\\master thesis\\R\\data\\samp.Rdata")

# find sample that experiences change in ADL
LSI <- (samp[,"chronicn"]==1)
indexinc2 <- matrix(FALSE,nrow(samp),1)
dur <- matrix(0,nrow(samp),1)
sphuslag <- matrix(NA,nrow(samp),1)
first <- matrix(FALSE,nrow(samp),1)
firstsphus <- matrix(NA,nrow(samp),1)
age <- c("age2004","age2007","age2011","age2013","age2015")
durationf <- c(1.5,1.5,2,2,1,1)
duration <- c(4,4,4,4,2,2)
for(i in 1:(nrow(data)/5)){
  indexinc <- (data[(i*5-4):(i*5),"chronicn"]<=1)
  if(sum(indexinc,na.rm=TRUE)==sum(!is.na(data[(i*5-4):(i*5),"chronicn"]))){
    obs <- which(!is.na(indexinc))
    if(sum(diff(obs))==(length(obs)-1)){ # check whether observations lie in consecutive waves
      indexinc <- (data[(i*5-4):(i*5),"chronicn"]==1)
      firstobs <- (i*5-5+obs[1])
      lastobs <- (i*5-5+obs[length(obs)])
      incobs <- which(LSI[firstobs:lastobs]==TRUE)
      waves <- data[firstobs:lastobs,"wave"][incobs]
      if(length(incobs)!=0){
        if(sum(diff(incobs))==(length(incobs)-1) & LSI[firstobs]==FALSE){
          indexinc2[(i*5-4):(i*5)] <- TRUE
          sphuslag[(i*5-3):(i*5)] <- data[(i*5-4):(i*5-1),"sphus"]
          firstsphus[(i*5-4):(i*5)] <- data[firstobs,"sphus"]
          first[firstobs] <- TRUE
          for(k in 1:length(incobs)){
            if(k==1){
              dur[(firstobs-1+incobs[k]),1] <- durationf[waves[k]]
            }else{
              dur[(firstobs-1+incobs[k]),1] <- duration[waves[k]] + dur[(firstobs-1+incobs[k-1]),1]
            }
          }
        }
      }
    }
  }
}

data2 <- cbind(data[,1:16],firstsphus,sphuslag,first,data[,17:18],LSI,dur,data[,19:ncol(data)])
samp <- data2[indexinc2,]
samp <- routing(samp)
rsamp <- samp[-which(is.na(samp[,"sphus"]) | is.na(samp[,"sphuslag"]) | is.na(samp[,"firstsphus"])),]

# remove imputed variables
#rsamp <- samp[,-which(colnames(samp)%in%c("mstatc","mstat2","yedu","child2","maxgrip2","remove"))]

# form imputations
#install.packages("Amelia")
#library(Amelia)
c <- which(colnames(rsamp)%in%c("numberid","id","wave","country","firstwave","mstat","age","isced",
                                "child","sphus","firstsphus","sphuslag","dur","LSI","IADL","employment","tenure"))
imputations <- amelia(rsamp[,c],m=1,ts="wave",cs="id",
                      idvars=c("numberid","country","firstwave"),
                      noms=c("mstat","employment","tenure"),
                      ords=c("sphus","sphuslag"))
imp <- imputations$imputations$imp1

############################################################################################################################
# marital status (married, single, widowed, divorced, separated)
imp[(imp[,"mstat"]==3 | imp[,"mstat"]==4),"mstat"] <- 1 # group Married and Registered partnership together
imp[(imp[,"mstat"]==5 | imp[,"mstat"]==7),"mstat"] <- 2 # group Separated and Divorced together
imp[(imp[,"mstat"]==6),"mstat"] <- 3 # Single
imp[(imp[,"mstat"]==8),"mstat"] <- 4 # Widowed
# labour-force status (employed, unemployed, retired, inactive)
imp[(imp[,"employment"]>=6),"employment"] <- imp[(imp[,"employment"]>=6),"employment"] <- 6
imp[,"employment"] <- imp[,"employment"]-2
# education (high school, less than high school, more than high school)
imp[(imp[,"isced"]<11 & imp[,"isced"]>0),"isced"] <- 2
imp[(imp[,"isced"]>12),"isced"] <- 4
imp[(imp[,"isced"]>=11 & imp[,"isced"]<=12),"isced"] <- 3
imp[(imp[,"isced"]<=0),"isced"] <- 1
# tenure
imp[(imp[,"tenure"]==3 | imp[,"tenure"]==4),"tenure"] <- 1 # owners
imp[(imp[,"tenure"]==5 | imp[,"tenure"]==6),"tenure"] <- 2 # renters
imp[(imp[,"tenure"]==7 | imp[,"tenure"]==8),"tenure"] <- 3 # other
# self-assessed health
imp[(imp[,"sphus"]==3),"sphus"] <- 4
imp[,"sphus"] <- abs(imp[,"sphus"] - 7)
imp[(imp[,"sphuslag"]==3),"sphuslag"] <- 4
imp[,"sphuslag"] <- abs(imp[,"sphuslag"] - 7)
imp[(imp[,"firstsphus"]==3),"firstsphus"] <- 4
imp[,"firstsphus"] <- abs(imp[,"firstsphus"] - 7)

# dynamic ordered probit
panel <- as.data.frame(cbind(imp[,c("numberid","id","wave","sphus","LSI","dur","age","child","isced")],
                             as.dummy(imp[,"sphuslag"],1,"sphuslag"),as.dummy(imp[,"firstsphus"],1,"firstsphus"),
                             as.dummy(imp[,"mstat"],3,"mstat"),as.dummy(imp[,"employment"],2,"employment"),
                             as.dummy(imp[,"tenure"],3,"tenure"),as.dummy(imp[,"isced"],1,"isced")))
panel <- fivewave(panel)
regressors <- c("LSI","dur","child",paste0("mstat","d",1:3),paste0("employment","d",1:3),paste0("tenure","d",1:2),
                paste0("isced","d",1:3)) 
lags <- c(paste0("sphuslag","d",1:3),paste0("firstsphus","d",1:3))
N <- nrow(panel)/5 # datch = samp
T <- 5
k <- length(regressors) # number of regressors
l <- length(lags)/2
cat <- 3 # number of categories - 1
y <- matrix(panel[,"sphus"],N,T,byrow=TRUE) # 0 = Poor, 1 = Fair, 2 = Good, 3 = Excellent
x <- array(NA,dim=c(N,T,(k+2*l)))
for(p in 1:l){
  x[,,p] <- matrix(panel[,lags[p]],N,T,byrow=TRUE) # lagged variable 1:3
}
for(p in 1:k){
  q <- l+p
  x[,,q] <- matrix(panel[,regressors[p]],N,T,byrow=TRUE) # all regressors 4:17
}
for(p in 1:l){
  q <- l+k+p
  x[,,q] <- matrix(panel[,lags[(l+p)]],N,T,byrow=TRUE) # first SAH observed in the sample 18:20
}

means <- eval(parse(text = paste0("array(c(",paste0("matrix(rowMeans(x[,,",4:(k+3),"],na.rm=TRUE),N,T,byrow=FALSE)",collapse=","),"),dim=c(N,T,k))")))
x <- array(c(x,means),dim=c(N,T,2*(k+l)))

# optimize log likelihood
gamma <- c(rep(0,(k+l)*2),1,2,3)
par <- optim(gamma,LL,method="BFGS",hessian=TRUE)
parameters <- par$par
names(parameters) <- c(paste0("sphuslag","d",1:3),regressors,paste0("firstsphus","d",1:3),paste0(regressors,"mean"),"tau1","tau2","tau3") 
Hessian <- par$hessian
cov <- solve(Hessian) 
sd <- sqrt(diag(cov))
z <- parameters/sd
pvalue <- pnorm(abs(z),lower.tail=FALSE)
star <- rep("",length(parameters))
star[pvalue < 0.05] <- "*"
star[pvalue < 0.01] <- "**"
star[pvalue < 0.001] <- "***"
show(as.data.frame(cbind(parameters,sd,z,pvalue,star))) 

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

