# panel data simulation
#install.packages("plm")
library(mvtnorm)
library(plm)
n <- 100
t <- 2
R <- 100
ols_c <- matrix(NA,3,R)
fixed_c <- matrix(NA,2,R)

for(r in 1:R){
x <- rmvnorm(n*t, c(5,1), matrix(c(1,0.5,0.5,1),nrow = 2))
alpha <- sort(rep(sqrt(t)*mean(x[,1]) + 1:n + rnorm(n,0,1),2))
beta <- c(1,-1)
y <- alpha + x%*%beta + rnorm(n*t,0,1)
individual <- sort(rep(1:n,2))
year <- rep(1:2,n)
data <- cbind(individual,year,y,x)
colnames(data) <- c("individual","year","y","x1","x2")
panel <- as.data.frame(data)

# run analyses
ols <- lm(y~x, data = panel)
fixed <- plm(y~x1 + x2, data = panel, index = c("individual","year"), model = "within")

ols_c[,r] <- ols$coefficients
fixed_c[,r] <- fixed$coefficients
}

ols_m <- rowMeans(ols_c)
fixed_m <- rowMeans(fixed_c)
