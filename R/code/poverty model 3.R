############################################################### construct data set #############################################################

load("C:\\Users\\anned\\OneDrive\\Documenten\\master thesis\\R\\data\\samp.Rdata")

# find sample that experiences change in ADL
incidence <- samp[,"IADL"] # "chronic" or "IADL"
indexinc <- (incidence>0) # > 0 for "IADL, == 3 for "chronic
indexinc2 <- matrix(FALSE,nrow(samp),1)
dur <- matrix(0,nrow(samp),1)
for(i in 1:(nrow(samp)/5)){
  obs <- which(!is.na(indexinc[(i*5-4):(i*5)]))
  if(sum(diff(obs))==(length(obs)-1)){ # check whether observations lie in consecutive waves
    firstobs <- (i*5-5+obs[1])
    lastobs <- (i*5-5+obs[length(obs)])
    incobs <- which(indexinc[firstobs:lastobs]==TRUE)
    if(length(incobs)!=0){
      if(sum(diff(incobs))==(length(incobs)-1) & indexinc[firstobs]==FALSE & indexinc[lastobs]==TRUE & sum(samp[firstobs:lastobs,"chronic"][incobs],na.rm=TRUE)==(length(incobs)*3)){ 
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
years <- c("dur0","dur25","dur5")
#years <- c("dur0","dur2","dur3","dur4")
c <- which(colnames(samp)%in%paste0("d",1:length(years)))
colnames(samp)[c] <- years
samp <- samp[indexinc2,]
datsamp <- which(indexinc2==TRUE)
indexsamp <- which(!is.na(samp[,"lifesat"]))
rsamp <- samp[-which(is.na(samp[,"lifesat"])),]
IADLi <- rsamp[,"IADL"]>0
rsamp <- cbind(rsamp,IADLi)
rsamp[,"IADLi"] <- as.numeric(rsamp[,"IADLi"])
#years <- "IADLi"

# form imputations
#install.packages("Amelia")
#library(Amelia)
c <- which(colnames(rsamp)%in%c("numberid","id","wave","country","firstwave","age","mstat","isced",
                                "child","casp","sphus",years,
                                "IADL","employment","lifesat"))
M <- round(100*length(unique(which(is.na(rsamp[,c]),arr.ind=TRUE)[,1]))/nrow(rsamp))
if(M<10){M <- 10}
imputations <- amelia(rsamp[,c],m=M,ts="wave",cs="id",
                      idvars=c("numberid","country","firstwave"),
                      noms=c("mstat","employment",years),
                      ords=c("sphus","lifesat"))
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
imp[,"sphus"] <- abs(imp[,"sphus"]-8)
# life satisfaction into 5 categories
imp[(imp[,"lifesat"]<=3),"lifesat"] <- 1
imp[(imp[,"lifesat"]<=4 & imp[,"lifesat"]>3),"lifesat"] <- 2
imp[(imp[,"lifesat"]<=5 & imp[,"lifesat"]>4),"lifesat"] <- 3
imp[(imp[,"lifesat"]<=6 & imp[,"lifesat"]>5),"lifesat"] <- 4
imp[(imp[,"lifesat"]<=7 & imp[,"lifesat"]>6),"lifesat"] <- 5
imp[(imp[,"lifesat"]<=8 & imp[,"lifesat"]>7),"lifesat"] <- 6
imp[(imp[,"lifesat"]<=9 & imp[,"lifesat"]>8),"lifesat"] <- 7
imp[(imp[,"lifesat"]<=10 & imp[,"lifesat"]>9),"lifesat"] <- 8

# estimate FE model
#install.packages("plm")
#library(plm)
panel <- as.data.frame(cbind(imp[,c("numberid","id","wave","sphus","IADL","lifesat",
                                years,"child")],
                         as.dummy(imp[,"age"],1,"age"),as.dummy(imp[,"mstat"],1,"mstat"),
                         as.dummy(imp[,"employment"],2,"employment"),as.dummy(imp[,"isced"],2,"isced")))
regressors <- c(paste0("mstat","d",1),paste0("employment","d",1:3),paste0("isced","d",1:2),"child",years,"IADL")
#regressors <- c(paste0("age","d",1:4),paste0("mstat","d",1),paste0("employment","d",1:3),paste0("isced","d",1:2),"child","IADLi","IADL")
formula <- as.formula(paste(paste("sphus ~"),paste(regressors,collapse="+")))
resultslist[[m]] <- plm(formula, data = panel, index = c("id","wave"), model = "within")
print(paste("m is",m))
}

##############################################################################################################################################
# calculate pooled estimate
output <- as.data.frame(matrix(NA,length(regressors),5))
rownames(output) <- regressors
colnames(output) <- c("Estimate","Std. Error","z value","Pr(>|z|)","star")

# calculate pooled regression coefficients
coef <- matrix(NA,length(regressors),M)
U <- matrix(NA,length(regressors),M)
for(m in 1:M){
  result <- summary(resultslist[[m]])$coefficients
  coef[,m] <- result[,1]
  U[,m] <- (result[,2])^2
  print(paste("m is",m))
}

mipool <- mi.meld(coef,sqrt(U),byrow=FALSE)
output[,"Estimate"] <- t(mipool$q.mi)
output[,"Std. Error"] <- t(mipool$se.mi)
output[,"z value"] <- output[,"Estimate"]/output[,"Std. Error"]
output[,"Pr(>|z|)"] <- pnorm(abs(output[,"z value"]),lower.tail=FALSE)*2
star <- rep("",length(regressors))
star[output[,"Pr(>|z|)"] < 0.05] <- "*"
star[output[,"Pr(>|z|)"] < 0.01] <- "**"
star[output[,"Pr(>|z|)"] < 0.001] <- "***"
output[,"star"] <- star
