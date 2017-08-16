set.seed(2)
#install.packages("plm")
library(plm)

# simulate data
N <- 100
T <- 3
R <- 100
k <- 5
ols_c <- matrix(NA,(k+1),R)
fixed_c <- matrix(NA,k,R)
fixedLL_c <- matrix(NA,k,R)
fixedNR_c <- matrix(NA,k,R)

for(r in 1:R){
  beta <- c(1,-1,0.03,5,2.45) # change this for different number of parameters
  if(length(beta) != k){stop("length beta is different from k")}
  x <- array(rnorm(N*T*k,0,1),dim = c(N,T,k))
  a <- rnorm(N,0,1)
  meanx_i <- rowMeans(x[,,1])
  alpha <- sqrt(T)*meanx_i + a + seq(1/10,(N/10),by=1/10)
  w <- matrix(alpha,nrow=N,ncol=T,byrow=FALSE) + eval(parse(text=paste(paste0("beta","[",1:k,"]*","x","[,,",1:k,"]",collapse="+"))))
  epsilon <- matrix(rnorm(N*T,0,1),nrow=N)
  y <- eval(parse(text=paste(paste0("beta","[",1:k,"]*","x","[,,",1:k,"]",collapse="+")))) + matrix(alpha,nrow=N,ncol=T,byrow=FALSE) +
       matrix(rnorm(N*T,0,1),N,T)
  sample1 <- sample(N,ceiling((3/4)*N))
  sample2 <- sample(N,ceiling((4/5)*N))
  x[-sample1,1:2,] <- NA
  x[-sample2,2,] <- NA
  y[-sample1,1:2] <- NA
  y[-sample2,2] <- NA
  
  # reshape data set
  y2 <- as.vector(t(y))
  x2 <- eval(parse(text=paste(paste0("cbind(",paste0("as.vector(t(x[,,",1:k,"]))",collapse=","),")"))))
  individual <- sort(rep(1:N,T))
  year <- rep(1:T,N)
  panel <- as.data.frame(cbind(individual,year,y2,x2))

  # run analyses
  formula <- as.formula(paste("y2~",paste0("x2[,",1:k,"]",collapse="+")))
  ols <- lm(formula)
  fixed <- plm(y2~x2, data = panel, index = c("individual","year"), model = "within")
  
  
  ols_c[,r] <- ols$coefficients
  fixed_c[,r] <- fixed$coefficients
  fixedLL_c[,r] <- optim(rep(0,k),LLna,"BFGS",x=x,y=y,N=N,T=T,k=k)$par
  fixedNR_c[,r] <- NR(rep(0,k),x,y,N,T,k)
  
  show(r)
}

ols_m <- rowMeans(ols_c)
fixed_m <- rowMeans(fixed_c)
fixedLL_m <- rowMeans(fixedLL_c) # perform different and worse from other two methods
fixedNR_m <- rowMeans(fixedNR_c)

# with actual data, get sample from sample code
# first transform data to contain five waves again
w <- c(1,2,4,5,6)
data <- imp[order(imp[,"numberid"]),]
for(i in 1:length(unique(data[,"numberid"]))){
  inddat <- which(data[,"numberid"]==unique(data[,"numberid"])[i])
  waves <- data[inddat,"wave"]
  for(k in 1:5){
    if(!(w[k]%in%waves)){
      if(inddat[1]==1 & k==1){
        datanew <- as.data.frame(matrix(NA,(nrow(data)+1),ncol(data)))
        colnames(datanew) <- colnames(data)
        datanew[2:nrow(datanew),] <- data[1:nrow(data),] 
        datanew[1,c(1:2,4:5)] <- data[inddat[1],c(1:2,4:5)]
        datanew[1,3] <- w[k]
        data <- datanew
      }
      else{
        indw <- inddat[1] - 2 + k
        datanew <- as.data.frame(matrix(NA,(nrow(data)+1),ncol(data)))
        colnames(datanew) <- colnames(data)
        datanew[1:indw,] <- data[1:indw,]
        datanew[(indw+2):nrow(datanew),] <- data[(indw+1):nrow(data),] 
        datanew[(indw+1),c(1:2,4:5)] <- data[inddat[1],c(1:2,4:5)]
        datanew[(indw+1),3] <- w[k]
        data <- datanew
      }
    } 
  }
  show(i)
}

regressors <- c("age","isced","child","casp","chronic","IADL")
#("gender","age","mstat","isced","child","casp","sphus","chronic","ADL","IADL","employment","maxgrip")
formula <- as.formula(paste(paste("sphus ~"),paste(regressors,collapse="+")))
ols <- lm(formula,data=datach)
fixed <- plm(formula, data = data, index = c("id","wave"), model = "within")

# homemade Newton Rhapson algorithm & maximum likelihood
N <- nrow(data)/5 # datch = samp
T <- 5
k <- length(regressors) # number of regressors
y <- matrix(data[,"sphus"],N,T,byrow=TRUE)
y <- y - 3
x <- array(NA,dim=c(N,T,k))
for(p in 1:k){
  x[,,p] <- matrix(data[,regressors[p]],N,T,byrow=TRUE)
}

fixedLL <- as.vector(t(optim(rep(0,k),LLna,"BFGS",x=xdat,y=ydat,N=N,T=T,k=k)$par))
names(fixedLL) <- regressors
fixedNR <- as.vector(t(NR(rep(0,k),x,y,N,T,k)))
names(fixedNR) <- regressors
