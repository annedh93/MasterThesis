load("C:\\Users\\anned\\OneDrive\\Documenten\\master thesis\\R\\data\\samp.Rdata")

response <- c("sphus")
if(response%in%"lifesat"){
  cat <- 8
}else{
  cat <- 5
}

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
    waves <- samp[firstobs:lastobs,"wave"][incobs]
    if(length(incobs)!=0){
      if(sum(diff(incobs))==(length(incobs)-1) & indexinc[firstobs]==FALSE & indexinc[lastobs]==TRUE & sum(samp[firstobs:lastobs,"chronic"][incobs],na.rm=TRUE)==(length(incobs)*3)){ 
        #length(unique(samp[firstobs:lastobs,"IADL"][incobs]))==1 & sum(diff(samp[firstobs:lastobs,"IADL"][incobs]))<=-length(incobs-1) & length(waves)>=3
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
dur2[which(dur2<=1 & dur2>0)] <- 1
dur2[which(dur2<=5.5 & dur2>1)] <- 2
dur2[which(dur2>5.5)] <- 3
samp <- cbind(samp[,1:23],as.dummy(dur2,2,""),samp[,24:ncol(samp)])
#samp <- cbind(samp[,1:23],dur2,samp[,24:ncol(samp)])
years <- c("dur0","dur15.5","dur5.5")
#years <- "dur2"
c <- which(colnames(samp)%in%paste0("d",1:length(years)))
colnames(samp)[c] <- years
samp <- samp[indexinc2,]
datsamp <- which(indexinc2==TRUE)
indexsamp <- which(!is.na(samp[,response]))
rsamp <- samp[-which(is.na(samp[,response])),]
rsamp[,"gender"] <- as.numeric(as.factor(rsamp[,"gender"]))

# form imputations
install.packages("Amelia")
library(Amelia)
c <- which(colnames(rsamp)%in%c("numberid","id","wave","country","firstwave","age","mstat","isced","gender",
                                "child","casp",response,years,
                                "IADL","employment"))
M <- round(100*length(unique(which(is.na(rsamp[,c]),arr.ind=TRUE)[,1]))/nrow(rsamp))
if(M<10){M <- 10}
imputations <- amelia(rsamp[,c],m=M,ts="wave",cs="id",
                      idvars=c("numberid","country","firstwave"),
                      noms=c("mstat","employment"),
                      ords=c(response))
implist <- imputations$imputations

#################################################### analysis ####################################################################################
imp <- implist[[1]]  

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
  if(response%in%"lifesat"){
    # life satisfaction into 5 categories
    imp[(imp[,"lifesat"]<=3),"lifesat"] <- 1
    imp[(imp[,"lifesat"]<=4 & imp[,"lifesat"]>3),"lifesat"] <- 2
    imp[(imp[,"lifesat"]<=5 & imp[,"lifesat"]>4),"lifesat"] <- 3
    imp[(imp[,"lifesat"]<=6 & imp[,"lifesat"]>5),"lifesat"] <- 4
    imp[(imp[,"lifesat"]<=7 & imp[,"lifesat"]>6),"lifesat"] <- 5
    imp[(imp[,"lifesat"]<=8 & imp[,"lifesat"]>7),"lifesat"] <- 6
    imp[(imp[,"lifesat"]<=9 & imp[,"lifesat"]>8),"lifesat"] <- 7
    imp[(imp[,"lifesat"]<=10 & imp[,"lifesat"]>9),"lifesat"] <- 8
  }else{
    # sphus from poor to excellent
    imp[,response] <- abs(imp[,response]-8)
  }
  
  # ordered logit
  panel <- as.data.frame(cbind(imp[,c("numberid","id","wave","IADL",response,"gender",
                                  years,"child")],
                           as.dummy(imp[,"mstat"],1,"mstat"),
                           as.dummy(imp[,"employment"],1,"employment"),
                           as.dummy(imp[,"isced"],2,"isced")))

# count changes in number of children and education
childcount <- 0
educcount <- 0
educcount2 <- 0
employcount <- 0
employcount2 <- 0
employcount3 <- 0
gendercount <- 0
IADLcount <- 0
not <- NULL
c <- unique(panel[,"numberid"])
k <- 0
for(i in 1:length(c)){
  if(sum(abs(diff(panel[which(panel[,"numberid"]==c[i]),"child"])))!=0){
    childcount <- childcount + 1
  }
  if(sum(abs(diff(panel[which(panel[,"numberid"]==c[i]),"iscedd1"])))!=0){
    educcount <- educcount + 1
    not <- cbind(not,c[i])
  }
  if(sum(abs(diff(panel[which(panel[,"numberid"]==c[i]),"iscedd2"])))!=0){
    educcount2 <- educcount2 + 1
  }
  if(sum(abs(diff(panel[which(panel[,"numberid"]==c[i]),"employmentd1"])))!=0){
    employcount <- employcount + 1
  }
  if(sum(abs(diff(panel[which(panel[,"numberid"]==c[i]),"employmentd2"])))!=0){
    employcount2 <- employcount2 + 1
  }
  if(sum(abs(diff(panel[which(panel[,"numberid"]==c[i]),"employmentd3"])))!=0){
    employcount3 <- employcount3 + 1
  }
  if(sum(abs(diff(panel[which(panel[,"numberid"]==c[i]),"gender"])))!=0){
    gendercount <- gendercount + 1
  }
  if(sum(abs(diff(panel[which(panel[,"numberid"]==c[i]),"IADL"])))!=0){
    IADLcount <- IADLcount + 1
  }
  k <- k + 1
  show(k)
}

regressors <- c(paste0("mstat","d",1),paste0("employment","d",1:3),paste0("isced","d",1:2),"child",years,"IADL","gender")
formula <- as.formula(paste(paste0(response," ~"),paste(regressors,collapse="+")))
result <- plm(formula, data = panel, index = c("id","wave"), model = "within")

# ordered logit
panel <- as.matrix(cbind(imp[,c("numberid","id","wave","IADL",response,"gender",
                                years,"child")],
                         as.dummy(imp[,"mstat"],1,"mstat"),
                         as.dummy(imp[,"employment"],1,"employment"),
                         as.dummy(imp[,"isced"],2,"isced")))
panel2 <- matrix(NA,nrow(samp),ncol(panel))
panel2[indexsamp,] <- panel
colnames(panel2) <- colnames(panel)
panel2 <- as.data.frame(panel2)
panel2[,c("numberid","id","wave")] <- samp[,c("numberid","id","wave")]
regressors <- c(paste0("mstat","d",1),paste0("employment","d",1:3),paste0("isced","d",1:2),"child",years,"IADL")
N <- nrow(panel2)/5 # datch = samp
T <- 5
k <- length(regressors) # number of regressors
y <- matrix(as.numeric(panel2[,response]),N,T,byrow=TRUE) # 0 = Poor, ..., 4 = Excellent
x <- array(NA,dim=c(N,T,k))
for(p in 1:k){
  x[,,p] <- matrix(as.numeric(panel2[,regressors[p]]),N,T,byrow=TRUE) # all regressors 
}

install.packages("gtools")
library(gtools)
result <- orderedlogitB(y,x,cat,regressors)
