data <-     retention(ph_waves,
                      cv_w1,cv_w2,cv_w3,cv_w4,cv_w5,
                      dn_w1,dn_w2,dn_w3,dn_w4,dn_w5,
                      ch_w1,ch_w2,ch_w3,ch_w4,ch_w5,
                      ph_w1,ph_w2,ph_w3,ph_w4,ph_w5,
                      ep_w1,ep_w2,ep_w3,ep_w4,ep_w5,
                      ac_w1,ac_w2,ac_w3,ac_w4,ac_w5,
                      ex_w1,ex_w2,ex_w3,ex_w4,ex_w5,
                      health_w1,health_w2,health_w3,health_w4,health_w5,
                      cv_retention,dn_retention,ch_retention,ph_retention,ep_retention,ac_retention,ex_retention,health_retention,
                      index1,index2,index3,index4,index5)

# make plot captureing relationship between ADL and self-reported health
cor(data[,"self perceived health"],data[,"ADL"],use="pairwise.complete.obs")
cor(data[,"self perceived health"],data[,"IADL"],use="pairwise.complete.obs")
cor(data[,"casp"],data[,"ADL"],use="pairwise.complete.obs")
cor(data[,"casp"],data[,"IADL"],use="pairwise.complete.obs")

# find sample that becomes ill during waves
chronic <- data[,"chronic"]
indexill <- (chronic==3) 
chronicindex <- matrix(FALSE,nrow(data),1)
for(i in 1:(nrow(data)/5)){
  obs <- which(!is.na(indexill[(i*5-4):(i*5)]))
  if(sum(diff(obs))==(length(obs)-1)){ # check whether observations lie in consecutive waves
    firstobs <- (i*5-5+obs[1])
    lastobs <- (i*5-5+obs[length(obs)])
    if(sum(indexill[firstobs:lastobs])>0 & sum(indexill[firstobs:lastobs])<length(obs) & indexill[firstobs]==FALSE){
      chronicindex[(i*5-4):(i*5)] <- TRUE
    }
  }
}

# form sample
chronicdata <- data[chronicindex,][,c(1,2,12,13,14,16,17)]
# make plots
ADLc <- table(chronicdata[(chronicdata[,"chronic"]==3),"ADL"])
ADLnc <- table(chronicdata[(chronicdata[,"chronic"]==4),"ADL"])
IADLc <- table(chronicdata[(chronicdata[,"chronic"]==3),"IADL"])
IADLnc <- table(chronicdata[(chronicdata[,"chronic"]==4),"IADL"])
tab <- matrix(0,4,10)
tab[c(1,3),1:7] <- t(cbind(ADLc,ADLnc))[,3:9]
tab[c(2,4),] <- t(cbind(IADLc,IADLnc))[,3:12]
colour <- c("palegreen","seagreen3","palevioletred1","orchid")
barplot(tab,col=colour,beside=TRUE)
legend("topright",c("chronic ADL","chronic IADL","not chornic ADL","not chronic IADL"),cex = 1,fill=colour)

wave <- chronicdata[,"wave"]
i <- 1
obs <- which(!is.na(chronicdata[(i*5-4):(i*5),"self perceived health"]))
firstobs <- (i*5-5+obs[1])
lastobs <- (i*5-5+obs[length(obs)])
chr <- which(chronicdata[(i*5-4):(i*5),"chronic"]==3)
firstchr <- (i*5-5+chr[1])
plot(1:5,chronicdata[1:5,"self perceived health"],type = "l",col="navy")
lines((chr[1]-1):chr[length(chr)],chronicdata[(firstchr-1):(firstchr),"self perceived health"],col="goldenrod1")

for(i in 2:(nrow(chronicdata)/5)){
  if(i==6){
    bla <- 0
  }
  obs <- which(!is.na(chronicdata[(i*5-4):(i*5),"self perceived health"]))
  firstobs <- (i*5-5+obs[1])
  lastobs <- (i*5-5+obs[length(obs)])
  chr <- which(chronicdata[(i*5-4):(i*5),"chronic"]==3)
  firstchr <- (i*5-5+chr[1])
  lines(obs[1]:(chr[1]-1),chronicdata[firstobs:(firstchr-1),"self perceived health"],col="navy")
  lines((chr[1]-1):obs[length(obs)],chronicdata[(firstchr-1):lastobs,"self perceived health"],col="goldenrod1")
}

# biplot 
#install.packages("ggplot2")
#library(devtools)
#install_github("vqv/ggbiplot")
library(ggbiplot)
illnessdata <- retention_diseases(ph_waves,ph_w1,ph_w2,ph_w4,ph_w5,ph_w6,ph_retention,index1,index2,index3,index4,index5)
ci <- illnessdata[chronicindex,]
ci[is.na(ci)] <- 2
res <- princomp(ci[,4:25],cor=TRUE)
biplot(res, choices=1:2, cex = c(1.5,0.5), xlabs = rep(".",nrow(ci[,4:25])), ylabs = c(colnames(ci[,4:10]),rep(".",15)), asp=1)
ggbiplot(res,choices = 1:2)
