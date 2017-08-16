data <-     retention(ph_waves,
                      cv_w1,cv_w2,cv_w3,cv_w4,cv_w5,
                      dn_w1,dn_w2,dn_w3,dn_w4,dn_w5,
                      ch_w1,ch_w2,ch_w3,ch_w4,ch_w5,
                      ph_w1,ph_w2,ph_w3,ph_w4,ph_w5,
                      ep_w1,ep_w2,ep_w3,ep_w4,ep_w5,
                      health_w1,health_w2,health_w3,health_w4,health_w5,
                      imp_w1,imp_w2,imp_w3,imp_w4,imp_w5,
                      cv_retention,dn_retention,ch_retention,ph_retention,ep_retention,health_retention,imp_retention,
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

# biplot 
#install.packages("ggplot2")
#library(devtools)
#install_github("vqv/ggbiplot")
library(ggbiplot)
illnessdata <- retention_diseases(ph_waves,ph_w1,ph_w2,ph_w4,ph_w5,ph_w6,ph_retention,index1,index2,index3,index4,index5)
ci <- illnessdata[chronicindex,]
# individuals 11659 30006 33984, only ones who filled out 1 (refuse)
ind <- which(ci[,4]==1)
ind2 <- which(ci[,4]==2)
ci <- ci[-c(ind,ind2),]
ci <- ci[-which(rowSums(ci[,4:24],na.rm=TRUE)==0),] 
ci[is.na(ci)] <- 3
res <- princomp(ci[,4:24],cor=TRUE)
#biplot(res, choices=1:2, cex = c(1.5,0.5), xlabs = rep(".",nrow(ci[,4:24])), ylabs = c(colnames(ci[,4:10]),rep(".",15)), asp=1)
ggbiplot(res,choices = 1:2,alpha=0.5,varname.adjust=1.5)

#nonlinear Pca
#ci <- ci[,4:ncol(ci)]
ci <- ci[,c(4,5,6,7,8,9,13,14,15,16,17)]
library(homals)                   # Load the homals package for nonlinear PCA
source("C:\\Users\\anned\\Downloads\\plot.homals.R") # Overrides default homals plot command
source("C:\\Users\\anned\\Downloads\\rescale.R") 
res <- homals(data = ci, ndim = 2, rank = 1, level = "ordinal")
res <- rescale(res)
plot(res, plot.type = "biplot", asp = 1, cex.lab = 1, main = NULL)     # biplot of component loadings and object scores: 
# only available in my own plot.homals!
plot(res, plot.type = "objplot", asp = 1)    # Plot of the object scores
plot(res, plot.type = "loadplot", asp = 1)   # Loadings plot
plot(res, plot.type = "trfplot",  col = "blue") # All transformation plots
plot(res, plot.type = "jointplot", asp = 1)  # biplot of objects and categories  
plot(res, plot.type = "jointplot", var.subset = 1, asp = 1)  # biplot of objects and categories  
plot(res, plot.type = "labplot", var.subset = 21, asp = 1) # object score plot labeled by a12 (doctors are murderers)
