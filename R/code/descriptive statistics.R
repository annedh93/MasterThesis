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
    waves <- samp[firstobs:lastobs,"wave"][incobs]
    if(length(incobs)!=0){
      if(sum(diff(incobs))==(length(incobs)-1) & indexinc[firstobs]==FALSE & indexinc[lastobs]==TRUE & sum(samp[firstobs:lastobs,"chronic"][incobs],na.rm=TRUE)==(length(incobs)*3)){ 
        #length(unique(samp[firstobs:lastobs,"IADL"][incobs]))==1 & length(unique(indexinc[incobs]))==1
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
samp <- cbind(samp[,1:23],as.dummy(dur2,2,""),dur,samp[,24:ncol(samp)])
#years <- c("<2","<5.5",">5.5")
years <- c("0","25.5",">5.5")
c <- which(colnames(samp)%in%paste0("d",1:length(years)))
colnames(samp)[c] <- years
rsamp <- samp[indexinc2,]

# marital status (married, single, widowed, divorced, separated)
rsamp[which(!is.na(rsamp[,"mstat"]) & rsamp[,"mstat"]==3 | rsamp[,"mstat"]==4),"mstat"] <- 1 # group Married and Registered partnership together
rsamp[which(!is.na(rsamp[,"mstat"]) & rsamp[,"mstat"]>=5),"mstat"] <- 2
# labour-force status (employed, unemployed, retired, inactive)
rsamp[which(!is.na(rsamp[,"employment"]) & (rsamp[,"employment"]>=6)),"employment"] <- 6
rsamp[,"employment"] <- rsamp[,"employment"]-2
# education (high school, less than high school, more than high school)
rsamp[which(!is.na(rsamp[,"isced"]) & rsamp[,"isced"]<11),"isced"] <- 1
rsamp[which(!is.na(rsamp[,"isced"]) & rsamp[,"isced"]>=11 & rsamp[,"isced"]<=12),"isced"] <- 2
rsamp[which(!is.na(rsamp[,"isced"]) & rsamp[,"isced"]>12),"isced"] <- 3
# sphus from poor to excellent
rsamp[,"sphus"] <- abs(rsamp[,"sphus"]-8)
# life satisfaction into 5 categories
rsamp[which(!is.na(rsamp[,"lifesat"]) & rsamp[,"lifesat"]<=3),"lifesat"] <- 1
rsamp[which(!is.na(rsamp[,"lifesat"]) & rsamp[,"lifesat"]<=4 & rsamp[,"lifesat"]>3),"lifesat"] <- 2
rsamp[which(!is.na(rsamp[,"lifesat"]) & rsamp[,"lifesat"]<=5 & rsamp[,"lifesat"]>4),"lifesat"] <- 3
rsamp[which(!is.na(rsamp[,"lifesat"]) & rsamp[,"lifesat"]<=6 & rsamp[,"lifesat"]>5),"lifesat"] <- 4
rsamp[which(!is.na(rsamp[,"lifesat"]) & rsamp[,"lifesat"]<=7 & rsamp[,"lifesat"]>6),"lifesat"] <- 5
rsamp[which(!is.na(rsamp[,"lifesat"]) & rsamp[,"lifesat"]<=8 & rsamp[,"lifesat"]>7),"lifesat"] <- 6
rsamp[which(!is.na(rsamp[,"lifesat"]) & rsamp[,"lifesat"]<=9 & rsamp[,"lifesat"]>8),"lifesat"] <- 7
rsamp[which(!is.na(rsamp[,"lifesat"]) & rsamp[,"lifesat"]<=10 & rsamp[,"lifesat"]>9),"lifesat"] <- 8

panel <- as.matrix(cbind(rsamp[,c("IADL","dur",years,"age","child","chronic")],as.numeric(as.factor(rsamp[,"gender"])),as.numeric(rsamp[,"IADL"]>0),
                         as.dummy2(rsamp[,"lifesat"],"lifesat"),as.dummy2(rsamp[,"sphus"],"sphus"),
                         as.dummy2(rsamp[,"mstat"],"mstat"),
                         as.dummy2(rsamp[,"employment"],"employment"),as.dummy2(rsamp[,"isced"],"isced")))
panel[,"chronic"] <- abs(panel[,"chronic"]-4)
panel[,"as.numeric(as.factor(rsamp[, \"gender\"]))"] <- abs(panel[,"as.numeric(as.factor(rsamp[, \"gender\"]))"]-2)
descriptives <- cbind(apply(panel,2,mean,na.rm=TRUE),apply(panel,2,sd,na.rm=TRUE))
colnames(descriptives) <- c("mean","sd")

dissphus <- matrix(0,5,5)
wave <- c(1,2,4,5,6)
k <- 1
for(i in wave){
  dissphus[,k] <- table(rsamp[which(rsamp[,"wave"]==i),"sphus"])/sum(table(rsamp[which(rsamp[,"wave"]==i),"sphus"]))  
  k <- k + 1
}

dislifesat <- matrix(0,8,5)
k <- 2
for(i in wave[2:length(wave)]){
  dislifesat[,k] <- table(rsamp[which(rsamp[,"wave"]==i),"lifesat"])/sum(table(rsamp[which(rsamp[,"wave"]==i),"lifesat"]))
  k <- k + 1
}

cumsat <- cbind(descriptives[10:17,],rep(0,8))
cumprev <- 0
for(i in 1:8){
  cumsat[i,3] <- cumsat[i,1] + cumprev 
  cumprev <- cumprev + cumsat[i,1]
}
cumsat2 <- cbind(descriptives[18:22,],rep(0,5))
cumprev <- 0
for(i in 1:5){
  cumsat2[i,3] <- cumsat2[i,1] + cumprev 
  cumprev <- cumprev + cumsat2[i,1]
}