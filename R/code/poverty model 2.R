############################################################### construct data set #############################################################
# C_it <- age (16-20,21-30,31-40,41-50,51-60,61-70,71-80,80+), 
#         marital status (married, single, widowed, divorced, separated), 
#         labour-force status (employed, unemployed, retired, inactive), 
#         education (high school, less than high school, more than high school),
#         number of children in the household and wave dummies
# d_0 <- incidence illness
# d_1 <- intensity illness
remove(data,rsamp,samp)
data <-     retention(id,
                      cv_w1,cv_w2,cv_w3,cv_w4,cv_w5,
                      dn_w1,dn_w2,dn_w3,dn_w4,dn_w5,
                      ch_w1,ch_w2,ch_w3,ch_w4,ch_w5,
                      ph_w1,ph_w2,ph_w3,ph_w4,ph_w5,
                      ep_w1,ep_w2,ep_w3,ep_w4,ep_w5,
                      ac_w1,ac_w2,ac_w3,ac_w4,ac_w5,
                      ho_w1,ho_w2,ho_w3,ho_w4,ho_w5,
                      health_w1,health_w2,health_w3,health_w4,health_w5,
                      imp_w1,imp_w2,imp_w3,imp_w4,imp_w5,
                      isced_w1,isced_w2,isced_w3,isced_w4,isced_w5,
                      tech_w1,tech_w2,tech_w3,tech_w4,tech_w5,
                      cv_retention,dn_retention,ch_retention,ph_retention,ep_retention,ac_retention,
                      ho_retention,health_retention,imp_retention,isced_retention,tech_retention,
                      index1,index2,index3,index4,index5)
remove(cv_w1,cv_w2,cv_w3,cv_w4,cv_w5,
       cv_1,cv_2,cv_3,cv_4,cv_5,
       dn_w1,dn_w2,dn_w3,dn_w4,dn_w5,
       ch_w1,ch_w2,ch_w3,ch_w4,ch_w5,
       ph_w1,ph_w2,ph_w3,ph_w4,ph_w5,
       ep_w1,ep_w2,ep_w3,ep_w4,ep_w5,
       ac_w1,ac_w2,ac_w3,ac_w4,ac_w5,
       ho_w1,ho_w2,ho_w3,ho_w4,ho_w5,
       health_w1,health_w2,health_w3,health_w4,health_w5,
       isced_w1,isced_w2,isced_w3,isced_w4,isced_w5,
       tech_w1,tech_w2,tech_w3,tech_w4,tech_w5,
       cv_retention,dn_retention,ch_retention,ph_retention,ep_retention,ac_retention,
       ho_retention,health_retention,isced_retention,tech_retention)

load("C:\\Users\\anned\\OneDrive\\Documenten\\master thesis\\R\\data\\samp.Rdata")

# find sample that experiences change in ADL
incidence <- samp[,"IADL"] # "chronic" or "IADL"
indexinc <- (incidence>0) # > 0 for "IADL, == 3 for "chronic
indexinc2 <- matrix(FALSE,nrow(samp),1)
dur <- matrix(0,nrow(samp),4)
colnames(dur) <- c("chronic01","chronic12","chronic23","chronic34")
for(i in 1:(nrow(samp)/5)){
  obs <- which(!is.na(indexinc[(i*5-4):(i*5)]))
  if(sum(diff(obs))==(length(obs)-1)){ # check whether observations lie in consecutive waves
    firstobs <- (i*5-5+obs[1])
    lastobs <- (i*5-5+obs[length(obs)])
    incobs <- which(indexinc[firstobs:lastobs]==TRUE)
    if(length(incobs)!=0){
      if(sum(diff(incobs))==(length(incobs)-1) & indexinc[firstobs]==FALSE){
        indexinc2[(i*5-4):(i*5)] <- TRUE
        for(k in 1:length(incobs)){
          dur[(firstobs-1+incobs[k]),k] <- 1
        }
      }
    }
  }
}
samp <- cbind(samp[,1:23],dur,samp[,24:ncol(samp)])
samp <- samp[indexinc2,]
indexsamp <- which(!is.na(samp[,"sphus"]))
rsamp <- samp[-which(is.na(samp[,"sphus"])),] # "sphus" or "lifesat"

# find sample that experiences change in ADL
dur <- matrix(0,nrow(samp),1)
for(i in 1:(nrow(samp)/5)){
  obs <- which(!is.na(indexinc[(i*5-4):(i*5)]))
  if(sum(diff(obs))==(length(obs)-1)){ # check whether observations lie in consecutive waves
    firstobs <- (i*5-5+obs[1])
    lastobs <- (i*5-5+obs[length(obs)])
    incobs <- which(indexinc[firstobs:lastobs]==TRUE)
    waves <- samp[firstobs:lastobs,"wave"][incobs]
    if(length(incobs)!=0){
      if(sum(diff(incobs))==(length(incobs)-1) & indexinc[firstobs]==FALSE & indexinc[lastobs]==TRUE & length(unique(samp[firstobs:lastobs,"IADL"][incobs]))==1 &
         sum(samp[firstobs:lastobs,"chronic"][incobs],na.rm=TRUE)==(length(incobs)*3)){ # & length(unique(indexinc[incobs]))==1
        indexinc2[(i*5-4):(i*5)] <- TRUE
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

dur2 <- dur
dur2[which(dur2<=2 & dur2>0)] <- 1
dur2[which(dur2<=5.5 & dur2>2)] <- 2
dur2[which(dur2>5.5)] <- 3
samp <- cbind(samp[,1:23],as.dummy(dur2,2,""),samp[,24:ncol(samp)])
#years <- c("<2","<5.5",">5.5")
years <- c("0","5.5",">5.5")
c <- which(colnames(samp)%in%paste0("d",1:length(years)))
colnames(samp)[c] <- years
samp <- samp[indexinc2,]
datsamp <- which(indexinc2==TRUE)
indexsamp <- which(!is.na(samp[,"lifesat"]))
rsamp <- samp[-which(is.na(samp[,"lifesat"])),]

# form imputations
#install.packages("Amelia")
#library(Amelia)
c <- which(colnames(rsamp)%in%c("numberid","id","wave","country","firstwave","age","mstat","isced",
                                "child","casp","sphus","chronic01","chronic12","chronic23",
                               "chronic34","IADL","employment","lifesat"))
M <- round(100*length(unique(which(is.na(rsamp),arr.ind=TRUE)[,1]))/nrow(rsamp))
imputations <- amelia(rsamp[,c],m=M,ts="wave",cs="id",
                      idvars=c("numberid","country","firstwave"),
                      noms=c("mstat","employment"),
                      ords=c("sphus"))
implist <- imputations$imputations

#################################################### analysis ####################################################################################
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
imp[(imp[,"mstat"]>=5),"mstat"] <- imp[(imp[,"mstat"]>=5),"mstat"]-3
# labour-force status (employed, unemployed, retired, inactive)
imp[(imp[,"employment"]>=6),"employment"] <- imp[(imp[,"employment"]>=6),"employment"] <- 6
imp[,"employment"] <- imp[,"employment"]-2
# education (high school, less than high school, more than high school)
imp[(imp[,"isced"])<11,"isced"] <- 1
imp[(imp[,"isced"]>=11 & imp[,"isced"]<=12),"isced"] <- 2
imp[(imp[,"isced"]>12),"isced"] <- 3

imp[,"sphus"] <- abs(imp[,"sphus"]-8)

# estimate FE model
#install.packages("plm")
#library(plm)
panel <- as.data.frame(cbind(imp[,c("numberid","id","wave","sphus","IADL","lifesat",
                                    "chronic01","chronic12","chronic23","chronic34","child")],
                             as.dummy(imp[,"age"],1,"age"),as.dummy(imp[,"mstat"],1,"mstat"),
                             as.dummy(imp[,"employment"],2,"employment"),as.dummy(imp[,"isced"],2,"isced")))
panel1 <- as.matrix(cbind(imp[,c("numberid","id","wave","sphus","IADL","lifesat",
                                    "chronic01","chronic12","chronic23","chronic34","child")],
                             as.dummy(imp[,"age"],1,"age"),as.dummy(imp[,"mstat"],1,"mstat"),
                             as.dummy(imp[,"employment"],2,"employment"),as.dummy(imp[,"isced"],2,"isced")))
panel2 <- matrix(NA,nrow(samp),ncol(panel))
panel2[indexsamp,] <- panel1
colnames(panel2) <- colnames(panel1)
panel2 <- as.data.frame(panel2)
panel2[,c("numberid","id","wave")] <- samp[,c("numberid","id","wave")]
regressors <- c(paste0("age","d",1:4),paste0("mstat","d",1:4),paste0("employment","d",1:3),paste0("isced","d",1:2),
                "child","chronic01","chronic12","chronic23","chronic34","IADL") 
formula <- as.formula(paste(paste("sphus ~"),paste(regressors,collapse="+")))

resultslist[[m]] <- plm(formula, data = panel, index = c("id","wave"), model = "within")
print(paste("m is",m))
}
N <- nrow(panel2)/5 # datch = samp
T <- 5
k <- length(regressors) # number of regressors
y <- matrix(as.numeric(panel2[,"sphus"]),N,T,byrow=TRUE) - rowMeans(matrix(as.numeric(panel2[,"sphus"]),N,T,byrow=TRUE),na.rm=TRUE) # all regressors  # 0 = Poor, ..., 4 = Excellent
x <- array(NA,dim=c(N,T,k))
for(p in 1:k){
  x[,,p] <- matrix(as.numeric(panel2[,regressors[p]]),N,T,byrow=TRUE) - rowMeans(matrix(as.numeric(panel2[,regressors[p]]),N,T,byrow=TRUE),na.rm=TRUE) # all regressors 
}
w <- eval(parse(text=paste(paste0("result3$coefficients","[",1:k,"]*x[,,",1:k,"]",collapse="+"))))
LL <- sum(log(pnorm(y-w)),na.rm=TRUE)
