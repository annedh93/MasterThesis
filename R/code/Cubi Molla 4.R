############################################################### construct data set #############################################################
load("C:\\Users\\anned\\OneDrive\\Documenten\\master thesis\\R\\data\\samp.Rdata")

# find sample that experiences change in ADL
response <- c("lifesat")
if(response%in%"lifesat"){
  samp <- samp[-which(samp[,"wave"]==1),]
  t <- 4
  cat <- 7
}else{
  t <- 5
  cat <- 4
}
indexinc <- (samp[,"IADL"]>0) # > 0 for "IADL, == 3 for "chronic
indexinc2 <- matrix(FALSE,nrow(samp),1)
dur <- matrix(0,nrow(samp),1)
lag <- matrix(NA,nrow(samp),1)
firsti <- matrix(FALSE,nrow(samp),1)
first <- matrix(NA,nrow(samp),1)

for(i in 1:(nrow(samp)/t)){
  obs <- which(!is.na(indexinc[(i*t-(t-1)):(i*t)]))
  if(sum(diff(obs))==(length(obs)-1)){ # check whether observations lie in consecutive waves
    firstobs <- (i*t-t+obs[1])
    lastobs <- (i*t-t+obs[length(obs)])
    incobs <- which(indexinc[firstobs:lastobs]==TRUE)
    waves <- samp[firstobs:lastobs,"wave"][incobs]
    #if(response%in%"lifesat" & firstobs%in%(i*5-5+1)){
    #  firstobs <- firstobs[2:length(firstobs)]
    #  incobs <- incobs[2:length(incobs)]
    #}
    if(length(incobs)!=0){ # length(incobs)==length(samp[firstobs:lastobs,"IADL"]>0)
      if(sum(diff(incobs))==(length(incobs)-1) & indexinc[firstobs]==FALSE & indexinc[lastobs]==TRUE & sum(samp[firstobs:lastobs,"chronic"][incobs],na.rm=TRUE)==(length(incobs)*3)){ 
        #length(unique(samp[firstobs:lastobs,"IADL"][incobs]))==1 & length(unique(indexinc[incobs]))==1
        indexinc2[(i*t-t+1):(i*t)] <- TRUE
        lag[(i*t-t+2):(i*t)] <- samp[(i*t-t+1):(i*t-1),response]
        first[(i*t-t+1):(i*t)] <- samp[firstobs,response]
        firsti[firstobs] <- TRUE
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

samp <- cbind(samp[,1:25],first,lag,firsti,samp[,26:ncol(samp)],dur) #18/25
samp <- samp[indexinc2,]
indexsamp <- which(!is.na(samp[,response]) & !is.na(samp[,"lag"]) & !is.na(samp[,"first"]))
rsamp <- samp[-which(is.na(samp[,response]) | is.na(samp[,"lag"]) | is.na(samp[,"first"])),]
IADLi <- (rsamp[,"IADL"]>0)
rsamp <- cbind(rsamp,IADLi)
rsamp[,"IADLi"] <- as.numeric(rsamp[,"IADLi"])

# form imputations
#install.packages("Amelia")
#library(Amelia)
c <- which(colnames(rsamp)%in%c("numberid","id","wave","country","firstwave","mstat","age","isced",
                                "child","lifesat","sphus","first","lag","dur","IADL","IADLi","employment"))
M <- round(100*length(unique(which(is.na(rsamp[,c]),arr.ind=TRUE)[,1]))/nrow(rsamp))
if(M<10){M <- 10}
imputations <- amelia(rsamp[,c],m=M,ts="wave",cs="id",
                      idvars=c("numberid","country","firstwave"),
                      noms=c("mstat","employment"),
                      ords=c("sphus","lifesat","first","lag"))
implist <- imputations$imputations

############################################################################################################################
resultslist <- list()
for(m in 1:M){
  
imp <- implist[[m]]
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
# life satisfaction into 5 categories
imp[(imp[,"lifesat"]<=3),"lifesat"] <- 0
imp[(imp[,"lifesat"]>=4),"lifesat"] <- imp[(imp[,"lifesat"]>=4),"lifesat"] - 3

if(response%in%"sphus"){
  imp[,"lag"] <- abs(imp[,"lag"]-7)
  imp[,"first"] <- abs(imp[,"first"]-7)
}else{
  imp[(imp[,"lag"]<=3),"lag"] <- 0
  imp[(imp[,"lag"]>=4),"lag"] <- imp[(imp[,"lag"]>=4),"lag"] -3
  imp[(imp[,"first"]<=3),"first"] <- 0
  imp[(imp[,"first"]>=4),"first"] <- imp[(imp[,"first"]>=4),"first"] -3
}


# dynamic ordered probit
panel <- as.matrix(cbind(imp[,c("numberid","id","wave","sphus","lifesat","IADL","IADLi","dur","child")],
                         as.dummy(imp[,"lag"],1,"lag"),as.dummy(imp[,"first"],1,"first"),
                         as.dummy(imp[,"mstat"],1,"mstat"),as.dummy(imp[,"employment"],2,"employment"),
                         as.dummy(imp[,"age"],1,"age"),as.dummy(imp[,"isced"],2,"isced")))
panel2 <- matrix(NA,nrow(samp),ncol(panel))
panel2[indexsamp,] <- panel
colnames(panel2) <- colnames(panel)
panel2 <- as.data.frame(panel2)
panel2[,c("numberid","id","wave")] <- samp[,c("numberid","id","wave")]
panel2 <- as.matrix(panel2)
regressors <- c(paste0("mstat","d",1),paste0("employment","d",1:3),paste0("isced","d",1:2),"child","IADL","IADLi","dur")
lags <- c(paste0("lag","d",1:cat),paste0("first","d",1:cat))
N <- nrow(panel2)/t 
T <- t
k <- length(regressors) # number of regressors
l <- length(lags)/2
y <- matrix(as.numeric(panel2[,response]),N,T,byrow=TRUE) # 0 = Poor, 1 = Fair, 2 = Good, 3 = Excellent
x <- array(NA,dim=c(N,T,(k+2*l)))
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

means <- eval(parse(text = paste0("array(c(",paste0("matrix(rowMeans(x[,,",(l+1):(l+k),"],na.rm=TRUE),N,T,byrow=FALSE)",collapse=","),"),dim=c(N,T,k))")))
X <- array(c(x,means),dim=c(N,T,2*(k+l)))

# optimize log likelihood
#install.packages("maxLik")
#library(maxLik)
totreg <- 2*(k+l)
gamma <- c(rep(0,totreg),1:cat)
A <- matrix(0,(cat-1),length(gamma))
for(k in 1:(cat-1)){
  A[k,(totreg+k):(totreg+k+1)] <- c(-1,1) 
}
B <- matrix(rep(0,(cat-1)),(cat-1),1)
par <- maxBFGS(LL,start = gamma,constraints = list(ineqA=A,ineqB=B),control = list(reltol = (10^-3)),finalHessian=TRUE) # constraints = list(ineqA=A,ineqB=B)
parameters <- par$estimate
Hessian <- par$hessian
cov <- -solve(Hessian) 
sd <- sqrt(diag(cov))
z <- parameters/sd
pvalue <- pnorm(abs(z),lower.tail=FALSE)
resultslist[[m]] <- as.matrix(cbind(parameters,sd,z,pvalue))

print(paste("m is",m))
}

###########################################################################################################################################

# calculate pooled estimate
output <- as.data.frame(matrix(NA,(totreg+cat),5))
rownames(output) <- c(paste0("lag","d",1:cat),regressors,paste0("first","d",1:cat),paste0(regressors,"mean"),paste0("tau",1:(cat)))
colnames(output) <- c("Estimate","Std. Error","z value","Pr(>|t|)","star")

# calculate pooled regression coefficients
coef <- matrix(NA,(totreg+cat),M)
U <- matrix(NA,(totreg+cat),M)
for(m in 1:M){
  result <- resultslist[[m]]
  coef[,m] <- result[,1]
  U[,m] <- (result[,2])^2
}
Una <- unique(which(is.na(U),arr.ind=TRUE)[,2])
if(length(Una)>0){
  coef <- coef[,-Una]
  U <- U[,-Una] 
}

mipool <- mi.meld(coef,sqrt(U),byrow=FALSE)
output[,"Estimate"] <- t(mipool$q.mi)
output[,"Std. Error"] <- t(mipool$se.mi)
output[,"z value"] <- output[,"Estimate"]/output[,"Std. Error"]
output[,"Pr(>|t|)"] <- pnorm(abs(output[,"z value"]),lower.tail=FALSE)*2
star <- rep("",(totreg+cat))
star[output[,"Pr(>|t|)"] < 0.05] <- "*"
star[output[,"Pr(>|t|)"] < 0.01] <- "**"
star[output[,"Pr(>|t|)"] < 0.001] <- "***"
output[,"star"] <- star

##############################################################################################################################################

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

