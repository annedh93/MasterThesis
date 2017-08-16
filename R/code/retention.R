retention_cv <- function(ph_unique,cv_w1,cv_w2,cv_w4,cv_w5,cv_w6,cv_retention){
  output <- cbind.data.frame(ph_unique,
                                    as.character(cv_w6[cv_retention[,5],"mergeid"]),
                                    cv_w6[cv_retention[,4],"country"],
                                    cv_w6[cv_retention[,4],"firstwave"],
                                    cv_w1[cv_retention[,1],"gender"],
                                    cv_w2[cv_retention[,2],"gender"],
                                    cv_w4[cv_retention[,3],"gender"],
                                    cv_w5[cv_retention[,4],"gender"],
                                    cv_w6[cv_retention[,5],"gender"],
                                    cv_w1[cv_retention[,1],"age2004"],
                                    cv_w2[cv_retention[,2],"age2007"],
                                    cv_w4[cv_retention[,3],"age2011"],
                                    cv_w5[cv_retention[,4],"age2013"],
                                    cv_w6[cv_retention[,5],"age2015"])

colnames(output) <- cbind("id", "mergeid","country","firstwave","gender2004","gender2007","gender2011","gender2013","gender2015",
                                   "age2004","age2007","age2011","age2013","age2015")
K <- ncol(output) - 1
output[,"mergeid"] <- as.character(output[,"mergeid"])

for(k in 2:4){
  for(i in 1:4){
    nam <- colnames(output)[k]
    index <- which(is.na(output[,nam])) # indices of missing values in variable "nam"
    if(i > 2){
      var <- eval(parse(text=paste("cv_w",(i+1),sep="")))
      output[index,nam] <- as.character(var[cv_retention[index,i],nam])
    }else{
      var <- eval(parse(text=paste("cv_w",i,sep="")))
      output[index,nam] <- as.character(var[cv_retention[index,i],nam])
    }
  }
}

return(output)
}

retention_diseases <- function(ph_unique,ph_w1,ph_w2,ph_w4,ph_w5,ph_w6,ph_retention){
output <- cbind.data.frame(ph_unique,
                           ph_w1[ph_retention[,1],"ph004_"],
                           ph_w2[ph_retention[,2],"ph004_"],
                           ph_w4[ph_retention[,3],"ph004_"],
                           ph_w5[ph_retention[,4],"ph004_"],
                           ph_w6[ph_retention[,5],"ph004_"],
                                       ph_w1[ph_retention[,1],"ph006d1"],
                                       ph_w2[ph_retention[,2],"ph006d1"],
                                       ph_w4[ph_retention[,3],"ph006d1"],
                                       ph_w5[ph_retention[,4],"ph006d1"],
                                       ph_w6[ph_retention[,5],"ph006d1"],
                                       ph_w1[ph_retention[,1],"ph006d2"],
                                       ph_w2[ph_retention[,2],"ph006d2"],
                                       ph_w4[ph_retention[,3],"ph006d2"],
                                       ph_w5[ph_retention[,4],"ph006d2"],
                                       ph_w6[ph_retention[,5],"ph006d2"],
                                       ph_w1[ph_retention[,1],"ph006d3"],
                                       ph_w2[ph_retention[,2],"ph006d3"],
                                       ph_w4[ph_retention[,3],"ph006d3"],
                                       ph_w5[ph_retention[,4],"ph006d3"],
                                       ph_w6[ph_retention[,5],"ph006d3"],
                                       ph_w1[ph_retention[,1],"ph006d4"],
                                       ph_w2[ph_retention[,2],"ph006d4"],
                                       ph_w4[ph_retention[,3],"ph006d4"],
                                       ph_w5[ph_retention[,4],"ph006d4"],
                                       ph_w6[ph_retention[,5],"ph006d4"],
                                       ph_w1[ph_retention[,1],"ph006d5"],
                                       ph_w2[ph_retention[,2],"ph006d5"],
                                       ph_w4[ph_retention[,3],"ph006d5"],
                                       ph_w5[ph_retention[,4],"ph006d5"],
                                       ph_w6[ph_retention[,5],"ph006d5"],
                                       ph_w1[ph_retention[,1],"ph006d6"],
                                       ph_w2[ph_retention[,2],"ph006d6"],
                                       ph_w4[ph_retention[,3],"ph006d6"],
                                       ph_w5[ph_retention[,4],"ph006d6"],
                                       ph_w6[ph_retention[,5],"ph006d6"],
                           ph_w1[ph_retention[,1],"ph006d7"],
                           ph_w2[ph_retention[,2],"ph006d7"],
                           ph_w4[ph_retention[,3],"ph006d7"],
                           ph_w1[ph_retention[,1],"ph006d8"],
                           ph_w2[ph_retention[,2],"ph006d8"],
                           ph_w4[ph_retention[,3],"ph006d8"],
                           ph_w1[ph_retention[,1],"ph006d9"],
                           ph_w2[ph_retention[,2],"ph006d9"],
                           ph_w4[ph_retention[,3],"ph006d9"],
                                       ph_w1[ph_retention[,1],"ph006d10"],
                                       ph_w2[ph_retention[,2],"ph006d10"],
                                       ph_w4[ph_retention[,3],"ph006d10"],
                                       ph_w5[ph_retention[,4],"ph006d10"],
                                       ph_w6[ph_retention[,5],"ph006d10"],
                                       ph_w1[ph_retention[,1],"ph006d11"],
                                       ph_w2[ph_retention[,2],"ph006d11"],
                                       ph_w4[ph_retention[,3],"ph006d11"],
                                       ph_w5[ph_retention[,4],"ph006d11"],
                                       ph_w6[ph_retention[,5],"ph006d11"],
                                       ph_w1[ph_retention[,1],"ph006d12"],
                                       ph_w2[ph_retention[,2],"ph006d12"],
                                       ph_w4[ph_retention[,3],"ph006d12"],
                                       ph_w5[ph_retention[,4],"ph006d12"],
                                       ph_w6[ph_retention[,5],"ph006d12"],
                                       ph_w1[ph_retention[,1],"ph006d13"],
                                       ph_w2[ph_retention[,2],"ph006d13"],
                                       ph_w4[ph_retention[,3],"ph006d13"],
                                       ph_w5[ph_retention[,4],"ph006d13"],
                                       ph_w6[ph_retention[,5],"ph006d13"],
                                       ph_w1[ph_retention[,1],"ph006d14"],
                                       ph_w2[ph_retention[,2],"ph006d14"],
                                       ph_w4[ph_retention[,3],"ph006d14"],
                                       ph_w5[ph_retention[,4],"ph006d14"],
                                       ph_w6[ph_retention[,5],"ph006d14"])
colnames(output) <- cbind("id",
                          "chronically ill 1","2","4","5","6",
                                   "heart attack 1","2","4","5","6",
                                   "high blood pressure or hypertension 1","2","4","5","6",
                                   "high blood cholesterol 1","2","4","5","6",
                                   "stroke or cerebral vascular disease 1","2","4","5","6",
                                   "diabetes or high blood sugar 1","2","4","5","6",
                                   "chronic lung disease 1","2","4","5","6",
                          "Asthma 1","2","4",
                          "Arthritis 1","2","4",
                          "Osteoporosis 1", "2","4",
                                   "cancer or malignant tumor 1","2","4","5","6",
                                   "stomach or duodenal ulcer","2","4","5","6",
                                   "Parkinson disease","2","4","5","6",
                                   "cataracts","2","4","5","6",
                                   "hip fracture","2","4","5","6")
return(output)
}

retention_health <- function(ph_unique,health_w1,health_w2,health_w4,health_w5,health_w6,health_retention){
          output <- cbind.data.frame(ph_unique,
                                    as.character(health_w6[health_retention[,5],"mergeid"]),
                                    health_w1[health_retention[,1],"casp"],
                                    health_w2[health_retention[,2],"casp"],
                                    health_w4[health_retention[,3],"casp"],
                                    health_w5[health_retention[,4],"casp"],
                                    health_w6[health_retention[,5],"casp"],
                                    health_w1[health_retention[,1],"sphus"],
                                    health_w2[health_retention[,2],"sphus"],
                                    health_w4[health_retention[,3],"sphus"],
                                    health_w5[health_retention[,4],"sphus"],
                                    health_w6[health_retention[,5],"sphus"])

colnames(output) <- cbind("id","mergeid","casp 1","casp 2", "casp 4", "casp 5","casp 6",
                          "self-perceived health 1","self-perceived health 2","self-perceived health 4","self-perceived health 5","self-perceived health 6")
output[,"mergeid"] <- as.character(output[,"mergeid"])
for(i in 1:4){
  nam <- "mergeid"
  index <- which(is.na(output[,nam])) # indices of missing values in variable "nam"
  if(i > 2){
    var <- eval(parse(text=paste("health_w",(i+1),sep="")))
    output[index,nam] <- as.character(var[health_retention[index,i],nam])
  }else{
    var <- eval(parse(text=paste("health_w",i,sep="")))
    output[index,nam] <- as.character(var[health_retention[index,i],nam])
  }
}

return(output)
}