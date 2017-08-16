data <-     retention(id,
                      cv_w1,cv_w2,cv_w3,cv_w4,cv_w5,
                      dn_w1,dn_w2,dn_w3,dn_w4,dn_w5,
                      ch_w1,ch_w2,ch_w3,ch_w4,ch_w5,
                      ph_w1,ph_w2,ph_w3,ph_w4,ph_w5,
                      ep_w1,ep_w2,ep_w3,ep_w4,ep_w5,
                      health_w1,health_w2,health_w3,health_w4,health_w5,
                      imp_w1,imp_w2,imp_w3,imp_w4,imp_w5,
                      isced_w1,isced_w2,isced_w3,isced_w4,isced_w5,
                      cv_retention,dn_retention,ch_retention,ph_retention,ep_retention,health_retention,imp_retention,isced_retention,
                      index1,index2,index3,index4,index5)

# find sample that experiences change in ADL
ADLIADL <- data[,"ADL"]+data[,"IADL"]
data <- cbind(data[,1:23],ADLIADL,data[,24:26])
ADL <- data[,"ADL"]
indexADL <- (ADL>0) 
indexADL2 <- matrix(FALSE,nrow(data),1)
for(i in 1:(nrow(data)/5)){
  obs <- which(!is.na(indexADL[(i*5-4):(i*5)]))
  if(sum(diff(obs))==(length(obs)-1)){ # check whether observations lie in consecutive waves
    firstobs <- (i*5-5+obs[1])
    lastobs <- (i*5-5+obs[length(obs)])
    if(sum(indexADL[firstobs:lastobs])>0 & indexADL[firstobs]==FALSE){
      indexADL2[(i*5-4):(i*5)] <- TRUE
    }
  }
}

samp <- data[indexADL2,]
i1 <- which(is.na(samp[,3])) # remove subjects not present in certain waves
i2 <- which(samp[,"age"]==-9) # remove deceased subjects
samp <- samp[-c(i1,i2),]


# convert the id to numbers
numberid <- as.numeric(as.factor(samp[,"id"]))
samp <- cbind(numberid,samp)
uniq <- unique(numberid)
ind <- NULL
ind2 <- NULL
count <- 0
count2 <- 0
count3 <- 0
count4 <- 0
count5 <- 0
count6 <- 0

# take care of routing
for(i in 1:length(uniq)){
  coup <- FALSE
  wavesp <- NULL
  obs <- samp[which(samp[,"numberid"]%in%uniq[i]),"mergeid"][1]
  indsamp <- which(samp[,"mergeid"]%in%obs)
  inddat <- which(data[,"mergeid"]%in%obs)
  waves <- samp[indsamp,"wave"]
  
  # count how many individuals have constant ADL/IADL
  nul <- which(samp[indsamp,"ADL"]==0)
  notnul <- which(samp[indsamp,"ADL"]!=0)
  lengthnul <- notnul[1]-nul[1]
  na <- which(is.na(samp[indsamp,"ADL"]))
  ADL <- samp[indsamp[-c(nul[1:lengthnul],na)],"ADL"]
  indADL <- indsamp[-c(nul[1:lengthnul],na)]
  
  if(sum(samp[indADL,"ADL"]==mean(samp[indADL,"ADL"]))==length(indADL)){
    count <- count + 1
    if(length(indADL)>1){
      count2 <- count2 + 1
      ind <- c(ind,indsamp)
    }
    else{count5 <- count5 + 1}
  }
  else{
      count3 <- count3 + 1
      if(sum(sort(ADL)==ADL)==length(ADL)){
        count6 <- count6 + 1
        ind <- c(ind,indsamp)
        }
  }
  if(length(indADL)>1){count4 <- count4 + 1}
    
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
      if(samp[indsamp[t],"mstatc"]==4 & !is.na(samp[indsamp[t],"mstatc"]) & is.na(samp[indsamp[t],"mstat"])){
        samp[indsamp[t],"mstat"] <- samp[(indsamp[t]-1),"mstat"]
      }
      if(is.na(samp[indsamp[t],"mstat"]) & waves[t]%in%wavesp & coup){
        tp <- which(wavesp%in%waves[t])
        samp[indsamp[t],"mstat"] <- data[indpartner[tp],"mstat"] 
      }

      # children
      if(is.na(samp[indsamp[t],"child"]) & waves[t]%in%wavesp & coup){
        tp <- which(wavesp%in%waves[t])
        samp[indsamp[t],"child"] <- data[indpartner[tp],"child"]
      }
      if(is.na(samp[indsamp[t],"child"])){
        samp[indsamp[t],"child"] <- samp[indsamp[which(!is.na(samp[indsamp,"child"]))][1],"child"]
      }
    
      # education
      if(is.na(samp[indsamp[t],"isced"])){
        if(t > 1){
          samp[indsamp[t],"isced"] <- samp[(indsamp[t]-1),"isced"]
        }
        if(is.na(samp[indsamp[t],"isced"]) & sum(!is.na(samp[indsamp,"isced"]))>0){
          samp[indsamp[t],"isced"] <- samp[indsamp[which(!is.na(samp[indsamp,"isced"]))][1],"isced"]
        }
      }
    
  }
  
  # children
  if(sum(is.na(samp[indsamp,"child"]))>0){
    samp[indsamp[which(is.na(samp[indsamp,"child"]))],"child"] <- samp[indsamp[which(!is.na(samp[indsamp,"child"]))][1],"child"]
  }
  
  show(i)
}

# remove imputed variables
c <- which(colnames(samp)%in%c("mstatc","mstat2","yedu","child2","maxgrip2"))
samp <- samp[,-c]
samp[,"gender"] <- as.numeric(as.factor(samp[,"gender"]))

# what missings are where
missings <- matrix(NA,1,ncol(samp))
rownames(missings) <- c("number missing")
colnames(missings) <- colnames(samp)
for(p in 1:ncol(samp)){
  missings[1,p] <- sum(is.na(samp[,p]))
}
show(missings)

# form imputations
install.packages("Amelia")
library(Amelia)
c <- which(colnames(samp)%in%c("numberid","id","wave","country","firstwave","gender","age","mstat","isced","child","casp","sphus","chronic",
                               "ADL","IADL","employment","maxgrip"))
imputations <- amelia(samp[,c],m=1,ts="wave",cs="id",
                      idvars=c("numberid","country","firstwave"),
                      noms=c("gender","mstat","chronic","employment"),
                      ords=c("sphus"))
imp <- imputations$imputations$imp1

# do some panel tests
install.packages("plm")
library(plm)
regressors <- c("gender","age","isced","child","casp","chronic","ADL","IADL","maxgrip")
formula <- as.formula(paste(paste("sphus ~"),paste(regressors,collapse="+")))
fixed <- plm(formula, data = imp, index = c("id","wave"), model = "within")
