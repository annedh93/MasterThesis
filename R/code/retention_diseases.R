retention_diseases <- function(ph_unique,ph_w1,ph_w2,ph_w4,ph_w5,ph_retention){
output <- cbind.data.frame(ph_unique,
                                       ph_w1[ph_retention[,1],"ph006d1"],
                                       ph_w2[ph_retention[,2],"ph006d1"],
                                       ph_w4[ph_retention[,3],"ph006d1"],
                                       ph_w5[ph_retention[,4],"ph006d1"],
                                       ph_w1[ph_retention[,1],"ph006d2"],
                                       ph_w2[ph_retention[,2],"ph006d2"],
                                       ph_w4[ph_retention[,3],"ph006d2"],
                                       ph_w5[ph_retention[,4],"ph006d2"],
                                       ph_w1[ph_retention[,1],"ph006d3"],
                                       ph_w2[ph_retention[,2],"ph006d3"],
                                       ph_w4[ph_retention[,3],"ph006d3"],
                                       ph_w5[ph_retention[,4],"ph006d3"],
                                       ph_w1[ph_retention[,1],"ph006d4"],
                                       ph_w2[ph_retention[,2],"ph006d4"],
                                       ph_w4[ph_retention[,3],"ph006d4"],
                                       ph_w5[ph_retention[,4],"ph006d4"],
                                       ph_w1[ph_retention[,1],"ph006d5"],
                                       ph_w2[ph_retention[,2],"ph006d5"],
                                       ph_w4[ph_retention[,3],"ph006d5"],
                                       ph_w5[ph_retention[,4],"ph006d5"],
                                       ph_w1[ph_retention[,1],"ph006d6"],
                                       ph_w2[ph_retention[,2],"ph006d6"],
                                       ph_w4[ph_retention[,3],"ph006d6"],
                                       ph_w5[ph_retention[,4],"ph006d6"],
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
                                       ph_w1[ph_retention[,1],"ph006d11"],
                                       ph_w2[ph_retention[,2],"ph006d11"],
                                       ph_w4[ph_retention[,3],"ph006d11"],
                                       ph_w5[ph_retention[,4],"ph006d11"],
                                       ph_w1[ph_retention[,1],"ph006d12"],
                                       ph_w2[ph_retention[,2],"ph006d12"],
                                       ph_w4[ph_retention[,3],"ph006d12"],
                                       ph_w5[ph_retention[,4],"ph006d12"],
                                       ph_w1[ph_retention[,1],"ph006d13"],
                                       ph_w2[ph_retention[,2],"ph006d13"],
                                       ph_w4[ph_retention[,3],"ph006d13"],
                                       ph_w5[ph_retention[,4],"ph006d13"],
                                       ph_w1[ph_retention[,1],"ph006d14"],
                                       ph_w2[ph_retention[,2],"ph006d14"],
                                       ph_w4[ph_retention[,3],"ph006d14"],
                                       ph_w5[ph_retention[,4],"ph006d14"])
colnames(output) <- cbind("id","heart attack 1","2","4","5",
                                   "high blood pressure or hypertension 1","2","4","5",
                                   "high blood cholesterol 1","2","4","5",
                                   "stroke or cerebral vascular disease 1","2","4","5",
                                   "diabetes or high blood sugar 1","2","4","5",
                                   "chronic lung disease 1","2","4","5",
                          "Asthma 1","2","4",
                          "Arthritis 1","2","4",
                          "Osteoporosis 1", "2","4",
                                   "cancer or malignant tumor 1","2","4","5",
                                   "stomach or duodenal ulcer","2","4","5",
                                   "Parkinson disease","2","4","5",
                                   "cataracts","2","4","5",
                                   "hip fracture","2","4","5")
return(output)
}

retention_health <- function(ph_unique,health_w1,health_w2,health_w4,health_w5,health_retention){
          output <- cbind.data.frame(ph_unique,
                                    as.character(health_w5[health_retention[,4],"mergeid"]),
                                    health_w1[health_retention[,1],"casp"],
                                    health_w2[health_retention[,2],"casp"],
                                    health_w4[health_retention[,3],"casp"],
                                    health_w5[health_retention[,4],"casp"],
                                    health_w1[health_retention[,1],"sphus"],
                                    health_w2[health_retention[,2],"sphus"],
                                    health_w4[health_retention[,3],"sphus"],
                                    health_w5[health_retention[,4],"sphus"])

colnames(output) <- cbind("id","mergeid","casp 1","casp 2", "casp 4", "casp 5",
                          "self-perceived health 1","self-perceived health 2","self-perceived health 4","self-perceived health 5")
output[,"mergeid"] <- as.character(output[,"mergeid"])
for(i in 1:3){
  nam <- "mergeid"
  index <- which(is.na(output[,nam])) # indices of missing values in variable "nam"
  if(i == 3){
    var <- eval(parse(text=paste("health_w",(i+1),sep="")))
    output[index,nam] <- as.character(var[health_retention[index,i],nam])
  }else{
    var <- eval(parse(text=paste("health_w",i,sep="")))
    output[index,nam] <- as.character(var[health_retention[index,i],nam])
  }
}

return(output)
}