data <-     retention(waves,
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

# routing marital status
w <- c(1,2,4,5,6)
for(k in 2:5){
  ind <- which(data[,"mstatc"]==4 & data[,"wave"]==w[k]) 
  data[ind,"mstat"] <- data[(ind-1),"mstat"]
}

# find sample that experiences change in ADL
ADL <- data[,"ADL"]
indexADL <- (ADL>0) 
indexADL2 <- matrix(FALSE,nrow(data),1)
for(i in 1:(nrow(data)/5)){
  obs <- which(!is.na(indexADL[(i*5-4):(i*5)]))
  if(sum(diff(obs))==(length(obs)-1)){ # check whether observations lie in consecutive waves
    firstobs <- (i*5-5+obs[1])
    lastobs <- (i*5-5+obs[length(obs)])
    if(sum(indexADL[firstobs:lastobs])>0 & sum(indexADL[firstobs:lastobs])<length(obs) & indexADL[firstobs]==FALSE){
      indexADL2[(i*5-4):(i*5)] <- TRUE
    }
  }
}

samp <- data[indexADL2,]

i <- which(is.na(samp[,3])) # remove subjects not present in certain waves
i2 <- which(samp[,"age"]==-9) # remove deceased subjects
samp <- samp[-c(i,i2),]

# routing children
nas <- which(is.na(samp[,"child"]))
ind <- unique(samp[nas,"mergeid"])
for(l in 1:length(ind)){
  obs <- ind[l]
  indsamp <- which(samp[,"mergeid"]%in%obs)
  indobs <- which(data[,"mergeid"]%in%obs)
  couple <- data[indobs[1],"coupleid"]
  if(!(couple%in%"")){
    indcouple <- which(data[,"coupleid"]%in%couple)
    min <- match(indobs,indcouple)
    min <- min[complete.cases(min)]
    indpart <- indcouple[-min]
    replace <- which(!is.na(data[indpart,"child"]))
    if(!is.na(indpart[replace[1]])){
      samp[indsamp,"child"] <- data[indpart[replace[1]],"child"] 
    }
  } 
  if(is.na(sum(samp[indsamp,"child"]))){
    nach <- which(is.na(samp[indsamp,"child"]))
    samp[indsamp[nach],"child"] <- samp[indsamp[-nach][1],"child"]
  }
} 

# routing marital status
for(k in 2:5){
  ind <- which(samp[,"mstatc"]==4 & samp[,"wave"]==w[k]) 
  samp[ind,"mstat"] <- samp[(ind-1),"mstat"]
}
nas <- which(is.na(samp[,"mstat"]))
ind <- unique(samp[nas,"mergeid"])
for(l in 1:length(ind)){
  obs <- ind[l]
  indsamp <- which(samp[,"mergeid"]%in%obs)
  indobs <- which(data[,"mergeid"]%in%obs)
  indnam <- which(is.na(samp[indsamp,"mstat"]))
  for(k in 1:length(indnam)){
    if(samp[indsamp[indnam],"mstatc"][k]==4 & !is.na(samp[indsamp[indnam],"mstatc"][k])){
      samp[indsamp[indnam],"mstat"][k] <- samp[(indsamp[indnam]-1),"mstat"][k]
    }
    if(is.na(samp[indsamp[indnam],"mstatc"][k])){
      couple <- samp[indsamp[indnam],"coupleid"][k]
      if(!(couple%in%"")){
        indcouple <- which(data[,"coupleid"]%in%couple)
        min <- match(indobs,indcouple)
        min <- min[complete.cases(min)]
        indpart <- indcouple[-min]
        replace <- which(!is.na(data[indpart,"child"]))
        if(!is.na(indpart[replace[1]])){
          samp[indsamp[indnam],"mstat"][k] <- data[indpart[replace[1]],"mstat"] 
        }
      }
    }
  }
}  

# remove imputed variables
samp <- samp[,-c(12,15,24,26)]
samp[,"gender"] <- as.numeric(as.factor(samp[,"gender"]))

# what missings are where
missings <- matrix(NA,1,ncol(samp))
rownames(missings) <- c("number missing")
colnames(missings) <- colnames(samp)
for(p in 1:ncol(samp)){
  missings[1,p] <- sum(is.na(samp[,p]))
}

# try form imputations
install.packages("Amelia")
library(Amelia)
imputations <- amelia(samp[,c(1:2,5:8,11,13,15:17,19:20)],m=1,ts="wave",cs="id",
                      idvars=c("country","firstwave"),
                      noms=c("gender","mstat","chronic"),
                      ords=c("sphus"))
imp <- imputations$imputations$imp1
